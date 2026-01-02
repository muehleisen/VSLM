# python2/vslm/gui/main.py
import sys
import os
from PySide6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                               QHBoxLayout, QPushButton, QLabel, QGroupBox, 
                               QFileDialog, QMessageBox, QFrame, QButtonGroup, 
                               QRadioButton, QProgressBar)
from PySide6.QtCore import Qt

# Imports
from ..analysis_engine import StreamProcessor
from .waveform import WaveformDialog
from .widgets import MatplotlibWidget

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("VSLM 2.0 (Python)")
        self.resize(1100, 750)
        
        # State
        self.filepath = None
        self.start_time = 0.0
        self.end_time = None # None means end of file
        self.cal_factor = 1.0
        self.block_size_ms = 100 # Default Analysis Block Size
        
        # UI Setup
        self._init_ui()
        self.status_bar = self.statusBar()
        self.status_bar.showMessage("Ready. Load a file to begin.")

    def _init_ui(self):
        # Main Widget
        central = QWidget()
        self.setCentralWidget(central)
        main_layout = QHBoxLayout(central)
        
        # --- LEFT PANEL (Controls) ---
        left_panel = QWidget()
        left_panel.setFixedWidth(300)
        left_layout = QVBoxLayout(left_panel)
        
        # 1. File Group
        grp_file = QGroupBox("File & Selection")
        layout_file = QVBoxLayout()
        
        self.btn_load = QPushButton("Load WAV File")
        self.btn_load.clicked.connect(self.on_load_file)
        
        self.btn_select = QPushButton("Select File Section")
        self.btn_select.clicked.connect(self.on_select_section)
        self.btn_select.setEnabled(False) # Disabled until file loaded
        
        self.lbl_info = QLabel("No File Loaded")
        self.lbl_info.setFrameStyle(QFrame.StyledPanel | QFrame.Sunken)
        self.lbl_info.setWordWrap(True)
        
        layout_file.addWidget(self.btn_load)
        layout_file.addWidget(self.btn_select)
        layout_file.addWidget(self.lbl_info)
        grp_file.setLayout(layout_file)
        left_layout.addWidget(grp_file)
        
        # 2. Settings Group (Weighting)
        grp_weight = QGroupBox("Weighting")
        layout_weight = QHBoxLayout()
        self.bg_weight = QButtonGroup()
        for i, text in enumerate(['A', 'C', 'Z']):
            rb = QRadioButton(text)
            if text == 'A': rb.setChecked(True)
            self.bg_weight.addButton(rb, i)
            layout_weight.addWidget(rb)
        grp_weight.setLayout(layout_weight)
        left_layout.addWidget(grp_weight)
        
        # 3. Settings Group (Speed)
        grp_speed = QGroupBox("Speed")
        layout_speed = QHBoxLayout()
        self.bg_speed = QButtonGroup()
        for i, text in enumerate(['Slow', 'Fast']):
            rb = QRadioButton(text)
            if text == 'Fast': rb.setChecked(True)
            self.bg_speed.addButton(rb, i)
            layout_speed.addWidget(rb)
        grp_speed.setLayout(layout_speed)
        left_layout.addWidget(grp_speed)
        
        # 4. Mode Group
        grp_mode = QGroupBox("Analysis Mode")
        layout_mode = QVBoxLayout()
        self.bg_mode = QButtonGroup()
        modes = ["Level vs Time", "Leq", "Octave Bands", "1/3 Octave Bands"]
        for i, m in enumerate(modes):
            rb = QRadioButton(m)
            if i == 0: rb.setChecked(True)
            self.bg_mode.addButton(rb, i)
            layout_mode.addWidget(rb)
        grp_mode.setLayout(layout_mode)
        left_layout.addWidget(grp_mode)
        
        left_layout.addStretch()
        
        # --- NEW: Progress Bar ---
        self.progress = QProgressBar()
        self.progress.setRange(0, 100)
        self.progress.setValue(0)
        self.progress.setTextVisible(False) # Minimal look
        # Optional styling
        self.progress.setStyleSheet("""
            QProgressBar { height: 10px; border: 1px solid grey; border-radius: 2px; } 
            QProgressBar::chunk { background-color: #3b82f6; }
        """)
        left_layout.addWidget(self.progress)
        
        # 5. Analyze Button
        self.btn_analyze = QPushButton("ANALYZE")
        self.btn_analyze.setStyleSheet("font-weight: bold; font-size: 14px; height: 40px; background-color: #dbeafe;")
        self.btn_analyze.clicked.connect(self.on_analyze)
        self.btn_analyze.setEnabled(False)
        left_layout.addWidget(self.btn_analyze)
        
        main_layout.addWidget(left_panel)
        
        # --- RIGHT PANEL (Results) ---
        self.plot_panel = MatplotlibWidget()
        main_layout.addWidget(self.plot_panel, stretch=1)

    # --- Actions ---

    def on_load_file(self):
        fname, _ = QFileDialog.getOpenFileName(self, "Open WAV", "", "WAV Files (*.wav)")
        if fname:
            self.filepath = fname
            self.start_time = 0.0
            
            # Get info
            try:
                from soundfile import info
                inf = info(fname)
                self.end_time = inf.duration
                
                # Update Info Label with Block Size
                self.lbl_info.setText(f"File: {os.path.basename(fname)}\n"
                                      f"Fs: {inf.samplerate} Hz\n"
                                      f"Dur: {inf.duration:.2f} s\n"
                                      f"Block: {self.block_size_ms} ms")
                
                self.btn_select.setEnabled(True)
                self.btn_analyze.setEnabled(True)
                self.status_bar.showMessage("File loaded.")
            except Exception as e:
                QMessageBox.critical(self, "Error", str(e))

    def on_select_section(self):
        if not self.filepath: return
        
        dlg = WaveformDialog(self.filepath, self)
        
        if self.end_time:
            dlg.viewer.region.setRegion([self.start_time, self.end_time])
            
        if dlg.exec():
            s, e = dlg.get_selection()
            self.start_time = s
            self.end_time = e
            
            current_text = self.lbl_info.text().split("\nSelection:")[0]
            self.lbl_info.setText(f"{current_text}\nSelection: {s:.2f}s - {e:.2f}s")

    def on_analyze(self):
        if not self.filepath: return
        
        # Gather Settings
        w_btn = self.bg_weight.checkedButton()
        weighting = w_btn.text() if w_btn else 'A'
        mode_id = self.bg_mode.checkedId() 
        
        self.status_bar.showMessage("Analyzing...")
        self.btn_analyze.setEnabled(False) 
        self.progress.setValue(0)
        
        # Force UI repaint before heavy loop
        QApplication.processEvents() 
        
        try:
            processor = StreamProcessor(self.filepath, cal_factor=self.cal_factor)
            
            do_bands = (mode_id >= 2)
            res = 'octave' if mode_id == 2 else 'third'
            
            # 1. Calculate Total Blocks for Progress
            # We assume processing the whole file (engine default) for now
            total_blocks = int(processor.duration * 1000 / self.block_size_ms)
            self.progress.setRange(0, total_blocks)
            
            # 2. Initialize Generator
            gen = processor.run_analysis(
                block_size_ms=self.block_size_ms, 
                weighting=weighting,
                do_band_analysis=do_bands,
                band_resolution=res
            )
            
            results = []
            
            # 3. Iterate and Update Progress
            for i, block in enumerate(gen):
                results.append(block)
                
                # Update progress bar
                self.progress.setValue(i + 1)
                
                # Allow GUI to update (prevent freezing)
                QApplication.processEvents()
            
            # Filter results by time selection
            filtered_results = [
                r for r in results 
                if self.start_time <= r['time'] <= self.end_time
            ]
            
            self._plot_results(filtered_results, mode_id, weighting)
            self.status_bar.showMessage("Analysis Complete.")
            
        except Exception as e:
            QMessageBox.critical(self, "Analysis Failed", str(e))
            self.status_bar.showMessage("Error.")
        finally:
            self.btn_analyze.setEnabled(True)

    def _plot_results(self, results, mode_id, weighting):
        self.plot_panel.reset()
        ax = self.plot_panel.ax
        
        if not results: return

        if mode_id == 0: # Level vs Time
            times = [r['time'] for r in results]
            levels = [r['leq'] for r in results]
            ax.plot(times, levels)
            ax.set_title(f"Level vs Time ({weighting}-Weighted)")
            ax.set_xlabel("Time (s)")
            ax.set_ylabel("Level (dB)")
            ax.grid(True)
            
        elif mode_id == 2 or mode_id == 3: # Bands
            import numpy as np
            
            freqs = results[0]['band_freqs']
            n_bands = len(freqs)
            energy_sums = np.zeros(n_bands)
            
            for r in results:
                levels = r['bands']
                pressures = (10**(levels/10.0)) * (20e-6**2)
                energy_sums += pressures
                
            mean_pressure = energy_sums / len(results)
            mean_db = 10 * np.log10(mean_pressure / (20e-6**2) + 1e-30)
            
            x = np.arange(len(freqs))
            ax.bar(x, mean_db)
            ax.set_xticks(x)
            
            lbls = []
            for f in freqs:
                if f >= 1000: lbls.append(f"{f/1000:.0f}k")
                else: lbls.append(f"{f:.0f}")
            
            if mode_id == 3:
                lbls = [l if i%3==0 else "" for i, l in enumerate(lbls)]
                
            ax.set_xticklabels(lbls, rotation=90)
            ax.set_title("Average Spectrum")
            ax.set_ylabel("Level (dB)")
            ax.grid(axis='y')

        self.plot_panel.draw()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    w = MainWindow()
    w.show()
    sys.exit(app.exec())