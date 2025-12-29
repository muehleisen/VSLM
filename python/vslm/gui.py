import sys
import os
import numpy as np
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                             QHBoxLayout, QPushButton, QLabel, QGroupBox, 
                             QRadioButton, QButtonGroup, QFileDialog, QMessageBox,
                             QProgressDialog, QFrame, QSizePolicy)
from PyQt6.QtCore import Qt, QThread, pyqtSignal, pyqtSlot

# Matplotlib Integration
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

# Import our backend modules
from .core import VSLMCore
from .playback import AudioPlayer
from . import analysis

class AnalysisWorker(QThread):
    """
    Runs analysis tasks in the background to prevent GUI freezing.
    Emits 'finished' signal with data when done.
    """
    result_ready = pyqtSignal(object)
    error_occurred = pyqtSignal(str)

    def __init__(self, function, *args, **kwargs):
        super().__init__()
        self.function = function
        self.args = args
        self.kwargs = kwargs

    def run(self):
        try:
            result = self.function(*self.args, **self.kwargs)
            self.result_ready.emit(result)
        except Exception as e:
            # traceback could be added here for debug
            self.error_occurred.emit(str(e))

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        
        self.setWindowTitle("VSLM - Virtual Sound Level Meter (Python Port)")
        self.resize(1000, 700)
        
        # Initialize Core Engines
        self.core = VSLMCore()
        self.player = AudioPlayer()
        self.worker = None # Placeholder for background thread
        
        # --- UI Setup ---
        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)
        self.main_layout = QHBoxLayout(self.central_widget)
        
        self.create_left_panel()
        self.create_right_panel()
        
        # Status Bar
        self.status_bar = self.statusBar()
        self.status_bar.showMessage("Ready. Please load a measurement file.")

    def create_left_panel(self):
        """Creates the control dashboard (Left side)"""
        panel = QWidget()
        layout = QVBoxLayout(panel)
        panel.setFixedWidth(280) # Fixed width like a tool palette
        
        # 1. File I/O Group
        io_group = QGroupBox("File & Calibration")
        io_layout = QVBoxLayout()
        
        self.btn_load_meas = QPushButton("Load Measurement (.wav)")
        self.btn_load_meas.clicked.connect(self.load_measurement)
        
        self.btn_load_cal = QPushButton("Set Calibration")
        self.btn_load_cal.clicked.connect(self.set_calibration)
        
        self.btn_play = QPushButton("Play Audio")
        self.btn_play.clicked.connect(self.toggle_playback)
        self.btn_play.setEnabled(False)
        
        io_layout.addWidget(self.btn_load_meas)
        io_layout.addWidget(self.btn_load_cal)
        io_layout.addWidget(self.btn_play)
        io_group.setLayout(io_layout)
        
        # 2. File Info Display
        self.info_label = QLabel("File: None\nLength: 0s\nFs: 0 Hz\nCal Factor: 1.0")
        self.info_label.setFrameStyle(QFrame.Shape.StyledPanel | QFrame.Shadow.Sunken)
        self.info_label.setStyleSheet("background-color: #f0f0f0; padding: 5px;")
        
        # 3. Weighting Group
        wtg_group = QGroupBox("Frequency Weighting")
        wtg_layout = QVBoxLayout()
        self.wtg_bg = QButtonGroup(self)
        
        rb_a = QRadioButton("A-Weighting")
        rb_c = QRadioButton("C-Weighting")
        rb_z = QRadioButton("Flat (Z)")
        rb_a.setChecked(True) # Default
        
        self.wtg_bg.addButton(rb_a, 1)
        self.wtg_bg.addButton(rb_c, 2)
        self.wtg_bg.addButton(rb_z, 3)
        
        wtg_layout.addWidget(rb_a)
        wtg_layout.addWidget(rb_c)
        wtg_layout.addWidget(rb_z)
        wtg_group.setLayout(wtg_layout)
        
        # 4. Speed Group
        spd_group = QGroupBox("Meter Speed")
        spd_layout = QVBoxLayout()
        self.spd_bg = QButtonGroup(self)
        
        rb_slow = QRadioButton("Slow (1.0s)")
        rb_fast = QRadioButton("Fast (125ms)")
        rb_imp = QRadioButton("Impulse (35ms)")
        rb_slow.setChecked(True)
        
        self.spd_bg.addButton(rb_slow, 1)
        self.spd_bg.addButton(rb_fast, 2)
        self.spd_bg.addButton(rb_imp, 3)
        
        spd_layout.addWidget(rb_slow)
        spd_layout.addWidget(rb_fast)
        spd_layout.addWidget(rb_imp)
        spd_group.setLayout(spd_layout)
        
        # 5. Analysis Mode Group
        mode_group = QGroupBox("Analysis Mode")
        mode_layout = QVBoxLayout()
        self.mode_bg = QButtonGroup(self)
        
        modes = [
            ("Sound Level (Lp)", "lp"),
            ("Leq / Dose", "leq"),
            ("Octave Band (ANSI)", "octave"),
            ("1/3 Octave (ANSI)", "third"),
            ("PSD (Power Spectral Density)", "psd"),
            ("Spectrogram", "spec")
        ]
        
        # Create radio buttons dynamically
        for i, (name, tag) in enumerate(modes):
            rb = QRadioButton(name)
            rb.setProperty("tag", tag) # Store string tag
            if i == 0: rb.setChecked(True)
            self.mode_bg.addButton(rb, i)
            mode_layout.addWidget(rb)
            
        mode_group.setLayout(mode_layout)
        
        # 6. Analyze Button
        self.btn_analyze = QPushButton("ANALYZE")
        self.btn_analyze.setFixedHeight(50)
        self.btn_analyze.setStyleSheet("font-weight: bold; font-size: 14px; background-color: #d0e0ff;")
        self.btn_analyze.clicked.connect(self.start_analysis)
        self.btn_analyze.setEnabled(False)
        
        # Add all to main layout
        layout.addWidget(io_group)
        layout.addWidget(self.info_label)
        layout.addWidget(wtg_group)
        layout.addWidget(spd_group)
        layout.addWidget(mode_group)
        layout.addStretch() # Push everything up
        layout.addWidget(self.btn_analyze)
        
        self.main_layout.addWidget(panel)

    def create_right_panel(self):
        """Creates the plotting area (Right side)"""
        right_widget = QWidget()
        layout = QVBoxLayout(right_widget)
        
        # Matplotlib Figure
        self.figure = Figure(figsize=(5, 4), dpi=100)
        self.canvas = FigureCanvas(self.figure)
        self.ax = self.figure.add_subplot(111)
        self.ax.text(0.5, 0.5, "Load a file to begin", ha='center', va='center')
        self.ax.axis('off')
        
        # Toolbar
        self.toolbar = NavigationToolbar(self.canvas, self)
        
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        
        self.main_layout.addWidget(right_widget, stretch=1)

    # --- Slots & Logic ---

    def load_measurement(self):
        fname, _ = QFileDialog.getOpenFileName(self, "Open Audio File", "", "Audio Files (*.wav)")
        if fname:
            try:
                self.core.load_file(fname)
                self.update_info_display()
                self.btn_analyze.setEnabled(True)
                self.btn_play.setEnabled(True)
                self.status_bar.showMessage(f"Loaded: {os.path.basename(fname)}")
            except Exception as e:
                QMessageBox.critical(self, "Error Loading File", str(e))

    def update_info_display(self):
        if self.core.audio_data is not None:
            duration = len(self.core.audio_data) / self.core.fs
            txt = (f"File: {os.path.basename(self.core.filename)}\n"
                   f"Length: {duration:.2f} s\n"
                   f"Fs: {self.core.fs} Hz\n"
                   f"Cal Factor: {self.core.cal_factor:.4f}")
            self.info_label.setText(txt)

    def set_calibration(self):
        # In VSLM, you can load a file OR type a number. Simple dialog for now:
        if self.core.audio_data is None:
            QMessageBox.warning(self, "Warning", "Please load a measurement file first (or a cal tone).")
            return
            
        val, ok = QFileDialog.getOpenFileName(self, "Load Calibration Tone (Optional)", "", "Audio Files (*.wav)")
        
        # If user picked a file, load it temporarily to get RMS
        ref_rms = None
        if ok and val:
             # Logic to read cal file RMS would go here (using soundfile)
             # For simplicity in this example, we ask for manual dB input
             pass

        from PyQt6.QtWidgets import QInputDialog
        db_val, ok = QInputDialog.getDouble(self, "Calibration", 
                                          "Enter Calibrator Level (dB):", 94.0, 0, 150, 1)
        if ok:
            try:
                # If a separate cal file wasn't loaded, we assume the CURRENT file is the tone
                # In a full app, you'd handle the separate cal file logic specifically
                factor = self.core.set_calibration(db_val) 
                self.update_info_display()
                QMessageBox.information(self, "Calibration", f"Factor set to {factor:.4f}")
            except Exception as e:
                QMessageBox.critical(self, "Error", str(e))

    def toggle_playback(self):
        if self.player.is_playing:
            self.player.stop()
            self.btn_play.setText("Play Audio")
        else:
            if self.core.audio_data is not None:
                self.player.play(self.core.audio_data, self.core.fs)
                self.btn_play.setText("Stop Audio")

    def get_selected_weighting(self):
        if self.wtg_bg.checkedId() == 1: return 'A'
        if self.wtg_bg.checkedId() == 2: return 'C'
        return 'Z'

    def get_selected_speed(self):
        if self.spd_bg.checkedId() == 1: return 'Slow'
        if self.spd_bg.checkedId() == 2: return 'Fast'
        return 'Impulse'

    def start_analysis(self):
        # Get settings
        mode_btn = self.mode_bg.checkedButton()
        if not mode_btn: return
        
        mode_tag = mode_btn.property("tag")
        weighting = self.get_selected_weighting()
        speed = self.get_selected_speed()
        
        self.status_bar.showMessage(f"Analyzing {mode_tag.upper()}... please wait.")
        self.btn_analyze.setEnabled(False)
        self.central_widget.setEnabled(False) # Lock UI
        
        # Define the task based on mode
        if mode_tag == 'lp':
            self.worker = AnalysisWorker(self.core.calculate_lp, weighting=weighting, speed=speed)
        elif mode_tag == 'leq':
            self.worker = AnalysisWorker(self.core.calculate_leq, weighting=weighting)
        elif mode_tag == 'octave':
            self.worker = AnalysisWorker(analysis.calculate_ansi_bands, self.core, weighting=weighting, resolution='octave')
        elif mode_tag == 'third':
            self.worker = AnalysisWorker(analysis.calculate_ansi_bands, self.core, weighting=weighting, resolution='third')
        elif mode_tag == 'psd':
            self.worker = AnalysisWorker(analysis.calculate_psd, self.core)
        elif mode_tag == 'spec':
            # Spectrogram logic usually returns raw matrix, we handle plotting in update_plot
            # For simplicity, we'll pass a lambda or specialized function if needed
            # But Matplotlib specgram computes it internally fast enough usually, 
            # OR we compute PSDs over time.
            # Let's just trigger the plot update directly for Spectrogram since Matplotlib does the math 
            # internally in `specgram`.
            self.handle_analysis_result(('spec', None)) 
            return 

        self.worker.result_ready.connect(lambda res: self.handle_analysis_result((mode_tag, res)))
        self.worker.error_occurred.connect(self.handle_analysis_error)
        self.worker.start()

    def handle_analysis_error(self, msg):
        self.central_widget.setEnabled(True)
        self.btn_analyze.setEnabled(True)
        self.status_bar.showMessage("Analysis Failed.")
        QMessageBox.critical(self, "Analysis Error", msg)

    def handle_analysis_result(self, payload):
        mode, result = payload
        
        # Unlock UI
        self.central_widget.setEnabled(True)
        self.btn_analyze.setEnabled(True)
        self.status_bar.showMessage("Analysis Complete.")
        
        # Clear Plot
        self.figure.clear()
        self.ax = self.figure.add_subplot(111)
        
        # --- Plotting Logic based on Mode ---
        
        if mode == 'lp':
            t, lp = result
            self.ax.plot(t, lp)
            self.ax.set_title(f"Sound Pressure Level ({self.get_selected_weighting()}-Weighted, {self.get_selected_speed()})")
            self.ax.set_xlabel("Time (s)")
            self.ax.set_ylabel("Lp (dB)")
            self.ax.grid(True)
            
        elif mode == 'leq':
            # Result is a single float number
            leq_val = result
            self.ax.axis('off')
            self.ax.text(0.5, 0.6, f"Leq ({self.get_selected_weighting()})", ha='center', fontsize=16)
            self.ax.text(0.5, 0.4, f"{leq_val:.2f} dB", ha='center', fontsize=30, fontweight='bold', color='blue')
            
        elif mode in ['octave', 'third']:
            freqs, levels = result
            # Bar chart
            x_pos = np.arange(len(freqs))
            self.ax.bar(x_pos, levels, width=0.8, color='green', alpha=0.7)
            self.ax.set_xticks(x_pos)
            
            # Smart tick labeling for 1/3 octave (skip some if too crowded)
            labels = [f"{int(f)}" if f < 1000 else f"{f/1000:.1f}k" for f in freqs]
            if mode == 'third':
                # Show every 3rd label
                for i in range(len(labels)):
                    if i % 3 != 0: labels[i] = ""
                    
            self.ax.set_xticklabels(labels, rotation=45)
            self.ax.set_title(f"{'Octave' if mode=='octave' else '1/3 Octave'} Band Levels")
            self.ax.set_ylabel("dB")
            self.ax.grid(axis='y')

        elif mode == 'psd':
            f, lpxx = result
            self.ax.semilogx(f, lpxx)
            self.ax.set_title("Power Spectral Density")
            self.ax.set_xlabel("Frequency (Hz)")
            self.ax.set_ylabel("dB/Hz")
            self.ax.grid(True, which="both")
            self.ax.set_xlim(20, self.core.fs/2)
            
        elif mode == 'spec':
            # Compute Spectrogram using Matplotlib's native optimized routine
            data = self.core.audio_data * self.core.cal_factor
            if self.get_selected_weighting() != 'Z':
                # Apply filter before specgram
                data = self.core.apply_weighting_filter(data, self.get_selected_weighting())
                
            Pxx, freqs, bins, im = self.ax.specgram(data, Fs=self.core.fs, NFFT=4096, noverlap=2048, cmap='inferno')
            self.ax.set_title(f"Spectrogram ({self.get_selected_weighting()}-Weighted)")
            self.ax.set_xlabel("Time (s)")
            self.ax.set_ylabel("Frequency (Hz)")
            self.figure.colorbar(im, ax=self.ax).set_label('Intensity (dB)')

        self.canvas.draw()

# --- Entry Point ---
if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())