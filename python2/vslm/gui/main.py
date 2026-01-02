# python2/vslm/gui/main.py
import sys
from pathlib import Path
from PySide6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                               QHBoxLayout, QPushButton, QLabel, QGroupBox, 
                               QFileDialog, QMessageBox, QFrame, QButtonGroup, 
                               QRadioButton, QProgressBar, QComboBox)

# Imports
from .waveform import WaveformDialog
from .calibration_dialog import CalibrationDialog
from .widgets import MatplotlibWidget
from .workers import AnalysisWorker
from .plotter import ResultPlotter # New Import

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("VSLM 2.0 (Python)")
        self.resize(1150, 800)
        
        # State
        self.filepath: Path | None = None
        self.start_time: float = 0.0
        self.end_time: float | None = None
        self.cal_factor: float = 1.0 
        self.block_size_ms: float = 100.0
        
        self.worker: AnalysisWorker | None = None
        
        self._init_ui()
        self.status_bar = self.statusBar()
        self.status_bar.showMessage("Ready. Load a file to begin.")

    def _init_ui(self):
        # ... (UI Layout Code - Unchanged from previous step) ...
        central = QWidget()
        self.setCentralWidget(central)
        main_layout = QHBoxLayout(central)
        
        # --- LEFT PANEL ---
        left_panel = QWidget()
        left_panel.setFixedWidth(320)
        left_layout = QVBoxLayout(left_panel)
        self.left_panel = left_panel 
        
        # 1. File Group
        grp_file = QGroupBox("File & Selection")
        layout_file = QVBoxLayout()
        self.btn_load = QPushButton("Load WAV File")
        self.btn_load.clicked.connect(self.on_load_file)
        self.btn_select = QPushButton("Select File Section")
        self.btn_select.clicked.connect(self.on_select_section)
        self.btn_select.setEnabled(False)
        self.btn_cal = QPushButton("Calibrate...")
        self.btn_cal.clicked.connect(self.on_calibrate)
        self.btn_cal.setStyleSheet("background-color: #f3f4f6;")
        
        self.lbl_info_header = QLabel("File Info")
        self.lbl_info_header.setStyleSheet("font-weight: bold; font-size: 10px; margin-top: 5px;")
        self.lbl_info = QLabel("No File Loaded")
        self.lbl_info.setFrameStyle(QFrame.StyledPanel | QFrame.Sunken)
        self.lbl_info.setStyleSheet("padding: 3px;")
        
        layout_file.addWidget(self.btn_load)
        layout_file.addWidget(self.btn_select)
        layout_file.addWidget(self.btn_cal) 
        layout_file.addWidget(self.lbl_info_header)
        layout_file.addWidget(self.lbl_info)
        grp_file.setLayout(layout_file)
        left_layout.addWidget(grp_file)
        
        # 2. Weighting
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
        
        # 3. Speed
        grp_speed = QGroupBox("Speed (Lp Mode)")
        layout_speed = QHBoxLayout()
        self.bg_speed = QButtonGroup()
        for i, text in enumerate(['Slow', 'Fast', 'Impulse']):
            rb = QRadioButton(text)
            if text == 'Fast': rb.setChecked(True)
            self.bg_speed.addButton(rb, i)
            layout_speed.addWidget(rb)
        grp_speed.setLayout(layout_speed)
        left_layout.addWidget(grp_speed)
        
        # 4. Mode
        grp_mode = QGroupBox("Analysis Mode")
        layout_mode = QVBoxLayout()
        self.bg_mode = QButtonGroup()
        modes = ["Level vs Time (Lp)", "LEQ Analysis", "Octave Bands", "1/3 Octave Bands"]
        for i, m in enumerate(modes):
            rb = QRadioButton(m)
            if i == 0: rb.setChecked(True)
            self.bg_mode.addButton(rb, i)
            layout_mode.addWidget(rb)
        grp_mode.setLayout(layout_mode)
        left_layout.addWidget(grp_mode)
        
        # 5. LEQ Settings
        grp_leq = QGroupBox("LEQ Settings")
        layout_leq = QHBoxLayout()
        layout_leq.addWidget(QLabel("Plot Interval:"))
        self.combo_leq_int = QComboBox()
        self.combo_leq_int.addItems(["100 ms", "1 sec", "10 sec", "1 min", "15 min", "1 hour"])
        self.combo_leq_int.setCurrentIndex(1)
        layout_leq.addWidget(self.combo_leq_int)
        grp_leq.setLayout(layout_leq)
        left_layout.addWidget(grp_leq)
        
        left_layout.addStretch()
        
        # Progress & Analyze
        self.progress = QProgressBar()
        self.progress.setRange(0, 100)
        self.progress.setValue(0)
        self.progress.setTextVisible(False)
        self.progress.setStyleSheet("QProgressBar { height: 10px; border: 1px solid grey; } QProgressBar::chunk { background-color: #3b82f6; }")
        left_layout.addWidget(self.progress)
        
        self.btn_analyze = QPushButton("ANALYZE")
        self.btn_analyze.setStyleSheet("font-weight: bold; font-size: 14px; height: 40px; background-color: #dbeafe;")
        self.btn_analyze.clicked.connect(self.on_analyze_click)
        self.btn_analyze.setEnabled(False)
        left_layout.addWidget(self.btn_analyze)
        
        main_layout.addWidget(left_panel)
        
        # --- RIGHT PANEL ---
        self.plot_panel = MatplotlibWidget()
        main_layout.addWidget(self.plot_panel, stretch=1)

    # --- Actions ---

    def _update_file_info_label(self, inf=None):
        if not self.filepath:
            self.lbl_info.setText(f"No File Loaded\nCal Factor: {self.cal_factor:.4f}")
            return
            
        if inf is None:
            from soundfile import info
            inf = info(str(self.filepath))
            
        self.lbl_info.setText(f"File: {self.filepath.name}\n"
                              f"Fs: {inf.samplerate} Hz\n"
                              f"Dur: {inf.duration:.1f} s   Block Size: {self.block_size_ms} ms\n"
                              f"Cal Factor: {self.cal_factor:.4f}")

    def on_load_file(self):
        fname, _ = QFileDialog.getOpenFileName(self, "Open WAV", "", "WAV Files (*.wav)")
        if fname:
            self.filepath = Path(fname)
            self.start_time = 0.0
            try:
                from soundfile import info
                inf = info(str(self.filepath))
                self.end_time = inf.duration
                self._update_file_info_label(inf)
                self.btn_select.setEnabled(True)
                self.btn_analyze.setEnabled(True)
                self.status_bar.showMessage("File loaded.")
            except Exception as e:
                QMessageBox.critical(self, "Error", str(e))

    def on_select_section(self):
        if not self.filepath: return
        dlg = WaveformDialog(str(self.filepath), self)
        if self.end_time:
            dlg.viewer.region.setRegion([self.start_time, self.end_time])
        if dlg.exec():
            s, e = dlg.get_selection()
            self.start_time = s
            self.end_time = e
            curr = self.lbl_info.text().split("\nSelection:")[0]
            self.lbl_info.setText(f"{curr}\nSelection: {s:.2f}s - {e:.2f}s")

    def on_calibrate(self):
        start = self.start_time
        end = self.end_time if self.end_time else 0.0
        dlg = CalibrationDialog(self.cal_factor, self.filepath, start, end, self)
        if dlg.exec():
            self.cal_factor = dlg.get_factor()
            self._update_file_info_label()
            self.status_bar.showMessage(f"Calibration updated: {self.cal_factor:.4f}")

    def toggle_inputs(self, enabled: bool):
        self.btn_load.setEnabled(enabled)
        self.btn_select.setEnabled(enabled)
        self.btn_cal.setEnabled(enabled) 
        self.combo_leq_int.setEnabled(enabled)
        for child in self.left_panel.findChildren(QGroupBox):
             if child.title() != "File & Selection": 
                 child.setEnabled(enabled)

    def on_analyze_click(self):
        if self.worker is not None:
            self.status_bar.showMessage("Stopping...")
            self.worker.stop()
            self.btn_analyze.setEnabled(False)
            return

        if not self.filepath: return
        
        w_btn = self.bg_weight.checkedButton()
        weighting = w_btn.text() if w_btn else 'A'
        
        s_btn = self.bg_speed.checkedButton()
        speed = s_btn.text() if s_btn else 'Fast'
        
        mode_id = self.bg_mode.checkedId()
        match mode_id:
            case 2:
                do_bands = True
                res = 'octave'
            case 3:
                do_bands = True
                res = 'third'
            case _:
                do_bands = False
                res = 'octave'
        
        self.toggle_inputs(False)
        self.btn_analyze.setText("STOP")
        self.btn_analyze.setStyleSheet("font-weight: bold; font-size: 14px; height: 40px; background-color: #fca5a5;")
        self.progress.setValue(0)
        self.status_bar.showMessage(f"Analyzing ({speed})...")
        
        self.worker = AnalysisWorker(
            self.filepath, 
            self.cal_factor, 
            self.block_size_ms,
            weighting,
            do_bands,
            res,
            speed
        )
        
        self.worker.sig_total_blocks.connect(self.progress.setMaximum)
        self.worker.sig_progress.connect(self.progress.setValue)
        self.worker.sig_finished.connect(lambda res: self.on_analysis_finished(res, mode_id, weighting, speed))
        self.worker.sig_error.connect(self.on_analysis_error)
        self.worker.finished.connect(self.on_worker_stopped)
        self.worker.finished.connect(self.worker.deleteLater)
        self.worker.start()

    def on_worker_stopped(self):
        self.worker = None
        self.toggle_inputs(True)
        self.btn_analyze.setText("ANALYZE")
        self.btn_analyze.setStyleSheet("font-weight: bold; font-size: 14px; height: 40px; background-color: #dbeafe;")
        self.btn_analyze.setEnabled(True)
        self.progress.setValue(0)

    def on_analysis_error(self, msg):
        QMessageBox.critical(self, "Analysis Error", msg)
        self.status_bar.showMessage("Error occurred.")

    def on_analysis_finished(self, results, mode_id, weighting, speed):
        self.status_bar.showMessage("Processing Results...")
        
        if self.end_time:
            filtered = [r for r in results if self.start_time <= r['time'] <= self.end_time]
        else:
            filtered = results

        self._plot_results(filtered, mode_id, weighting, speed)
        self.status_bar.showMessage("Analysis Complete.")

    def _plot_results(self, results: list, mode_id: int, weighting: str, speed: str):
        # Delegate to new ResultPlotter class
        leq_int_txt = self.combo_leq_int.currentText()
        ResultPlotter.plot(
            self.plot_panel.figure,
            results,
            mode_id,
            weighting,
            speed,
            leq_int_txt,
            self.block_size_ms
        )
        # Redraw
        self.plot_panel.draw()

    def closeEvent(self, event):
        if self.worker is not None and self.worker.isRunning():
            reply = QMessageBox.question(
                self, 'Analysis Running',
                "An analysis is currently running.\nDo you want to stop it and exit?",
                QMessageBox.Yes | QMessageBox.No, QMessageBox.No
            )
            if reply == QMessageBox.Yes:
                self.worker.stop()
                self.worker.wait()
                event.accept()
            else:
                event.ignore()
        else:
            event.accept()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    w = MainWindow()
    w.show()
    sys.exit(app.exec())