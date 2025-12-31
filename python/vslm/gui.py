# vslm/gui.py
import sys
import os
import inspect
from typing import Tuple, List, Optional, Any

import numpy as np
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, 
    QPushButton, QLabel, QGroupBox, QRadioButton, QButtonGroup, 
    QFileDialog, QMessageBox, QFrame, QProgressBar, QInputDialog, 
    QSizePolicy, QLayout
)
from PySide6.QtCore import Qt, QThread, Signal, QSize

# Matplotlib Integration
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

# Backend Imports
from .core import VSLMCore
from .playback import AudioPlayer
from . import analysis

# --- Constants & Styles ---

BTN_SIZE = QSize(50, 30)

STYLE_TOGGLE_BTN = """
    QPushButton {
        border: 2px solid #aaa;
        border-radius: 8px;
        background-color: #f5f5f5;
    }
    QPushButton:checked {
        background-color: #3b82f6;
        border-color: #1d4ed8;
    }
    QPushButton:hover {
        border-color: #3b82f6;
    }
"""

STYLE_INFO_LABEL = """
    QLabel {
        background-color: #f0f0f0; 
        padding: 5px; 
        border: 1px inset #ccc;
    }
"""

STYLE_ANALYZE_BTN = """
    QPushButton {
        font-weight: bold; 
        font-size: 14px; 
        background-color: #d0e0ff;
    }
"""

STYLE_PROGRESS_BAR = """
    QProgressBar { 
        height: 10px; 
        border: 1px solid grey; 
        border-radius: 2px; 
    } 
    QProgressBar::chunk { 
        background-color: #3b82f6; 
    }
"""


class AnalysisWorker(QThread):
    """
    Executes analysis functions in a background thread to keep the GUI responsive.
    """
    result_ready = Signal(object)
    error_occurred = Signal(str)
    progress_updated = Signal(int)

    def __init__(self, function: Any, *args: Any, **kwargs: Any) -> None:
        super().__init__()
        self.function = function
        self.args = args
        self.kwargs = kwargs

    def run(self) -> None:
        try:
            # Inspect function to see if it accepts a progress callback
            sig = inspect.signature(self.function)
            if 'progress_callback' in sig.parameters:
                self.kwargs['progress_callback'] = self.emit_progress
            
            result = self.function(*self.args, **self.kwargs)
            self.result_ready.emit(result)
        except Exception as e:
            self.error_occurred.emit(str(e))

    def emit_progress(self, value: int) -> None:
        self.progress_updated.emit(value)


class MainWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        
        self.setWindowTitle("VSLM - Virtual Sound Level Meter (Python Port)")
        self.resize(1050, 700)
        
        # Initialize Backend
        self.core = VSLMCore()
        self.player = AudioPlayer()
        self.worker: Optional[AnalysisWorker] = None 
        
        # Main UI Container
        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)
        self.main_layout = QHBoxLayout(self.central_widget)
        
        # UI Elements (initialized in create methods)
        self.info_label: QLabel
        self.btn_analyze: QPushButton
        self.btn_play: QPushButton
        self.wtg_bg: QButtonGroup
        self.spd_bg: QButtonGroup
        self.mode_bg: QButtonGroup
        self.progress_bar: QProgressBar
        
        self._create_left_panel()
        self._create_right_panel()
        
        self.status_bar = self.statusBar()
        self.status_bar.showMessage("Ready. Please load a measurement file.")

    # --- UI Construction ---

    def _create_left_panel(self) -> None:
        """Constructs the left-hand settings sidebar."""
        panel = QWidget()
        panel.setFixedWidth(280)
        layout = QVBoxLayout(panel)
        
        layout.addWidget(self._create_io_group())
        layout.addWidget(self._create_info_display())
        layout.addWidget(self._create_weighting_group())
        layout.addWidget(self._create_speed_group())
        layout.addWidget(self._create_mode_group())
        
        layout.addStretch() 
        layout.addWidget(self._create_progress_bar())
        layout.addWidget(self._create_analyze_button())
        
        self.main_layout.addWidget(panel)

    def _create_right_panel(self) -> None:
        """Constructs the right-hand plotting area."""
        right_widget = QWidget()
        layout = QVBoxLayout(right_widget)
        
        self.figure = Figure(figsize=(5, 4), dpi=100)
        self.canvas = FigureCanvas(self.figure)
        self.ax = self.figure.add_subplot(111)
        self.ax.text(0.5, 0.5, "Load a file to begin", ha='center', va='center')
        self.ax.axis('off')
        
        self.toolbar = NavigationToolbar(self.canvas, self)
        
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        
        self.main_layout.addWidget(right_widget, stretch=1)

    # --- Group Creators ---

    def _create_io_group(self) -> QGroupBox:
        group = QGroupBox("File & Calibration")
        layout = QVBoxLayout()
        
        btn_load = QPushButton("Load Measurement (.wav)")
        btn_load.clicked.connect(self.load_measurement)
        
        btn_cal = QPushButton("Set Calibration")
        btn_cal.clicked.connect(self.set_calibration)
        
        self.btn_play = QPushButton("Play Audio")
        self.btn_play.clicked.connect(self.toggle_playback)
        self.btn_play.setEnabled(False)
        
        layout.addWidget(btn_load)
        layout.addWidget(btn_cal)
        layout.addWidget(self.btn_play)
        group.setLayout(layout)
        return group

    def _create_info_display(self) -> QLabel:
        self.info_label = QLabel("File: None\nLength: 0s\nFs: 0 Hz\nCal Factor: 1.0")
        self.info_label.setFrameStyle(QFrame.Shape.StyledPanel | QFrame.Shadow.Sunken)
        self.info_label.setStyleSheet(STYLE_INFO_LABEL)
        return self.info_label

    def _create_weighting_group(self) -> QGroupBox:
        group = QGroupBox("Frequency Weighting")
        layout = QHBoxLayout()
        layout.setSpacing(10)
        
        self.wtg_bg = QButtonGroup(self)
        options = [("A", 1), ("C", 2), ("Flat (Z)", 3)]
        
        for text, uid in options:
            pair = self._create_vertical_pair(text, uid, self.wtg_bg)
            layout.addWidget(pair)
            
        group.setLayout(layout)
        return group

    def _create_speed_group(self) -> QGroupBox:
        group = QGroupBox("Meter Speed")
        layout = QHBoxLayout()
        layout.setSpacing(10)
        
        self.spd_bg = QButtonGroup(self)
        options = [
            ("Slow\n(1.0s)", 1), 
            ("Fast\n(125ms)", 2), 
            ("Impulse\n(35ms/1.5s)", 3)
        ]
        
        for text, uid in options:
            pair = self._create_vertical_pair(text, uid, self.spd_bg)
            layout.addWidget(pair)
            
        group.setLayout(layout)
        return group

    def _create_mode_group(self) -> QGroupBox:
        group = QGroupBox("Analysis Mode")
        layout = QVBoxLayout()
        self.mode_bg = QButtonGroup(self)
        
        modes = [
            ("Sound Level (Lp)", "lp"),
            ("Leq / Dose", "leq"),
            ("Octave Band (ANSI)", "octave"),
            ("1/3 Octave (ANSI)", "third"),
            ("PSD", "psd"),
            ("Spectrogram", "spec")
        ]
        
        for i, (name, tag) in enumerate(modes):
            rb = QRadioButton(name)
            rb.setProperty("tag", tag)
            if i == 0: rb.setChecked(True)
            self.mode_bg.addButton(rb, i)
            layout.addWidget(rb)
            
        group.setLayout(layout)
        return group

    def _create_progress_bar(self) -> QProgressBar:
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        self.progress_bar.setTextVisible(False)
        self.progress_bar.setStyleSheet(STYLE_PROGRESS_BAR)
        self.progress_bar.setVisible(False)
        return self.progress_bar

    def _create_analyze_button(self) -> QPushButton:
        self.btn_analyze = QPushButton("ANALYZE")
        self.btn_analyze.setFixedHeight(40)
        self.btn_analyze.setStyleSheet(STYLE_ANALYZE_BTN)
        self.btn_analyze.clicked.connect(self.start_analysis)
        self.btn_analyze.setEnabled(False)
        return self.btn_analyze

    def _create_vertical_pair(self, text: str, btn_id: int, 
                              group: QButtonGroup) -> QWidget:
        """Helper to create a 'Label above Button' widget pair."""
        container = QWidget()
        layout = QVBoxLayout(container)
        layout.setContentsMargins(0, 5, 0, 5)
        layout.setSpacing(4)
        
        lbl = QLabel(text)
        lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
        lbl.setStyleSheet("font-weight: bold; color: #444; font-size: 11px;")
        
        btn = QPushButton("")
        btn.setCheckable(True)
        btn.setFixedSize(BTN_SIZE)
        btn.setCursor(Qt.CursorShape.PointingHandCursor)
        btn.setStyleSheet(STYLE_TOGGLE_BTN)
        
        if btn_id == 1:
            btn.setChecked(True)
            
        group.addButton(btn, btn_id)
        
        layout.addWidget(lbl)
        layout.addWidget(btn, alignment=Qt.AlignmentFlag.AlignCenter)
        return container

    # --- Logic ---

    def load_measurement(self) -> None:
        fname, _ = QFileDialog.getOpenFileName(self, "Open Audio File", "", "Audio Files (*.wav)")
        if fname:
            try:
                self.core.load_file(fname)
                self._update_info_text()
                self.btn_analyze.setEnabled(True)
                self.btn_play.setEnabled(True)
                self.status_bar.showMessage(f"Loaded: {os.path.basename(fname)}")
            except Exception as e:
                QMessageBox.critical(self, "Error Loading File", str(e))

    def _update_info_text(self) -> None:
        if self.core.audio_data is not None:
            duration = len(self.core.audio_data) / self.core.fs
            txt = (f"File: {os.path.basename(self.core.filename)}\n"
                   f"Length: {duration:.2f} s\n"
                   f"Fs: {self.core.fs} Hz\n"
                   f"Cal Factor: {self.core.cal_factor:.4f}")
            self.info_label.setText(txt)

    def set_calibration(self) -> None:
        if self.core.audio_data is None:
            QMessageBox.warning(self, "Warning", "Please load a measurement file first.")
            return
            
        db_val, ok = QInputDialog.getDouble(self, "Calibration", 
                                          "Enter Calibrator Level (dB):", 94.0, 0, 150, 1)
        if ok:
            try:
                factor = self.core.set_calibration(db_val) 
                self._update_info_text()
                QMessageBox.information(self, "Calibration", f"Factor set to {factor:.4f}")
            except Exception as e:
                QMessageBox.critical(self, "Error", str(e))

    def toggle_playback(self) -> None:
        if self.player.is_playing:
            self.player.stop()
            self.btn_play.setText("Play Audio")
        else:
            if self.core.audio_data is not None:
                self.player.play(self.core.audio_data, self.core.fs)
                self.btn_play.setText("Stop Audio")

    def get_selected_weighting(self) -> str:
        uid = self.wtg_bg.checkedId()
        return {1: 'A', 2: 'C', 3: 'Z'}.get(uid, 'A')

    def get_selected_speed(self) -> str:
        uid = self.spd_bg.checkedId()
        return {1: 'Slow', 2: 'Fast', 3: 'Impulse'}.get(uid, 'Slow')

    def start_analysis(self) -> None:
        mode_btn = self.mode_bg.checkedButton()
        if not mode_btn: return
        
        mode_tag = mode_btn.property("tag")
        weighting = self.get_selected_weighting()
        speed = self.get_selected_speed()
        
        self.status_bar.showMessage(f"Analyzing {mode_tag.upper()}... please wait.")
        self._set_ui_busy(True)
        
        # Determine Task
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
            # Spectrogram is fast enough or handled differently; simple pass-through here
            self.handle_analysis_result(('spec', None)) 
            return 

        if self.worker:
            self.worker.result_ready.connect(lambda res: self.handle_analysis_result((mode_tag, res)))
            self.worker.error_occurred.connect(self.handle_analysis_error)
            self.worker.progress_updated.connect(self.progress_bar.setValue)
            self.worker.start()

    def _set_ui_busy(self, busy: bool) -> None:
        self.btn_analyze.setEnabled(not busy)
        self.central_widget.setEnabled(not busy)
        self.progress_bar.setVisible(busy)
        if busy:
            self.progress_bar.setValue(0)

    def handle_analysis_error(self, msg: str) -> None:
        self._set_ui_busy(False)
        self.status_bar.showMessage("Analysis Failed.")
        QMessageBox.critical(self, "Analysis Error", msg)

    def handle_analysis_result(self, payload: Tuple[str, Any]) -> None:
        mode, result = payload
        self._set_ui_busy(False)
        self.status_bar.showMessage("Analysis Complete.")
        
        self.figure.clear()
        self.ax = self.figure.add_subplot(111)
        
        weighting = self.get_selected_weighting()

        if mode == 'lp':
            t, lp = result
            self.ax.plot(t, lp)
            self.ax.set_title(f"Sound Pressure Level ({weighting}-Weighted, {self.get_selected_speed()})")
            self.ax.set_xlabel("Time (s)")
            self.ax.set_ylabel("Lp (dB)")
            self.ax.grid(True)
            
        elif mode == 'leq':
            leq_val = result
            self.ax.axis('off')
            self.ax.text(0.5, 0.6, f"Leq ({weighting})", ha='center', fontsize=16)
            self.ax.text(0.5, 0.4, f"{leq_val:.2f} dB", ha='center', fontsize=30, fontweight='bold', color='blue')
            
        elif mode in ['octave', 'third']:
            freqs, levels = result
            x_pos = np.arange(len(freqs))
            self.ax.bar(x_pos, levels, width=0.8, color='green', alpha=0.7)
            self.ax.set_xticks(x_pos)
            
            # Smart labeling
            labels = []
            for f in freqs:
                lbl = f"{int(f)}" if f < 1000 else f"{f/1000:.1f}k"
                labels.append(lbl)
                
            if mode == 'third':
                # Skip labels to avoid crowding
                labels = [lbl if i % 3 == 0 else "" for i, lbl in enumerate(labels)]
                    
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
            data = self.core.audio_data * self.core.cal_factor
            if weighting != 'Z':
                data = self.core.apply_weighting_filter(data, weighting)
                
            Pxx, freqs, bins, im = self.ax.specgram(
                data, Fs=self.core.fs, NFFT=4096, noverlap=2048, cmap='inferno'
            )
            self.ax.set_title(f"Spectrogram ({weighting}-Weighted)")
            self.ax.set_xlabel("Time (s)")
            self.ax.set_ylabel("Frequency (Hz)")
            self.figure.colorbar(im, ax=self.ax).set_label('Intensity (dB)')

        self.canvas.draw()


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())