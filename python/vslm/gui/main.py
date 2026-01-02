import sys
from pathlib import Path
from PySide6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                               QHBoxLayout, QPushButton, QLabel, QGroupBox, 
                               QFileDialog, QMessageBox, QFrame, QButtonGroup, 
                               QRadioButton, QProgressBar, QComboBox)
from PySide6.QtGui import QAction, QDesktopServices
from PySide6.QtCore import QUrl

# --- VSLM Imports ---
from .waveform import WaveformDialog
from .calibration_dialog import CalibrationDialog
from .about_dialog import AboutDialog 
from .widgets import MatplotlibWidget
from .workers import AnalysisWorker
from .plotter import ResultPlotter
from ..export import ResultsExporter
from ..settings import SettingsManager, AppSettings

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("VSLM 2.0 (Python)")
        self.resize(1150, 800)
        
        # 1. Load Default Settings
        self.settings_mgr = SettingsManager()
        self.settings = self.settings_mgr.load() 
        
        # 2. Application State
        self.filepath: Path | None = None
        self.start_time: float = 0.0
        self.end_time: float | None = None
        self.cal_factor: float = self.settings.calibration_factor
        self.block_size_ms: float = self.settings.block_size_ms
        
        self.last_results: list = [] 
        self.has_unsaved_data: bool = False # DIRTY FLAG
        self.worker: AnalysisWorker | None = None
        
        # 3. Setup UI
        self._init_menu_bar()
        self._init_ui()
        self._apply_settings_to_ui()
        
        self.status_bar = self.statusBar()
        self.status_bar.showMessage("Ready. Load a file to begin.")

    def _init_menu_bar(self):
        menu_bar = self.menuBar()
        
        # --- File Menu ---
        menu_file = menu_bar.addMenu("File")
        
        # Settings Submenu (Load/Save)
        menu_settings = menu_file.addMenu("Settings")
        
        act_load_sets = QAction("Load Settings...", self)
        act_load_sets.triggered.connect(self.on_action_load_settings)
        menu_settings.addAction(act_load_sets)
        
        act_save_sets = QAction("Save Settings...", self)
        act_save_sets.triggered.connect(self.on_action_save_settings)
        menu_settings.addAction(act_save_sets)
        
        menu_file.addSeparator()
        
        act_quit = QAction("Quit", self)
        act_quit.setShortcut("Ctrl+Q")
        act_quit.triggered.connect(self.close) # Triggers closeEvent
        menu_file.addAction(act_quit)
        
        # --- Export Menu ---
        self.menu_export = menu_bar.addMenu("Export")
        self.menu_export.setEnabled(False) 
        
        act_export_csv = QAction("Save Results (CSV)...", self)
        act_export_csv.triggered.connect(self.on_export_csv)
        self.menu_export.addAction(act_export_csv)
        
        act_save_fig = QAction("Save Plot Figure...", self)
        act_save_fig.triggered.connect(self.on_action_save_figure)
        self.menu_export.addAction(act_save_fig)
        
        # --- Help Menu ---
        menu_help = menu_bar.addMenu("Help")
        
        act_docs = QAction("Documentation", self)
        act_docs.triggered.connect(lambda: self.on_open_url("https://example.com/docs"))
        menu_help.addAction(act_docs)
        
        act_tuts = QAction("Tutorials", self)
        act_tuts.triggered.connect(lambda: self.on_open_url("https://example.com/tutorials"))
        menu_help.addAction(act_tuts)
        
        menu_help.addSeparator()
        
        act_about = QAction("About VSLM", self)
        act_about.triggered.connect(self.on_about)
        menu_help.addAction(act_about)

    def _init_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        main_layout = QHBoxLayout(central)
        
        # --- LEFT PANEL ---
        left_panel = QWidget()
        left_panel.setFixedWidth(320)
        left_layout = QVBoxLayout(left_panel)
        self.left_panel = left_panel 
        
        # Group 1: File & Calibration
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
        
        # Group 2: Weighting
        grp_weight = QGroupBox("Weighting")
        layout_weight = QHBoxLayout()
        self.bg_weight = QButtonGroup()
        for i, text in enumerate(['A', 'C', 'Z']):
            rb = QRadioButton(text)
            self.bg_weight.addButton(rb, i)
            layout_weight.addWidget(rb)
        grp_weight.setLayout(layout_weight)
        left_layout.addWidget(grp_weight)
        
        # Group 3: Speed
        grp_speed = QGroupBox("Speed (Lp Mode)")
        layout_speed = QHBoxLayout()
        self.bg_speed = QButtonGroup()
        for i, text in enumerate(['Slow', 'Fast', 'Impulse']):
            rb = QRadioButton(text)
            self.bg_speed.addButton(rb, i)
            layout_speed.addWidget(rb)
        grp_speed.setLayout(layout_speed)
        left_layout.addWidget(grp_speed)
        
        # Group 4: Mode
        grp_mode = QGroupBox("Analysis Mode")
        layout_mode = QVBoxLayout()
        self.bg_mode = QButtonGroup()
        modes = ["Level vs Time (Lp)", "LEQ Analysis", "Octave Bands", "1/3 Octave Bands"]
        for i, m in enumerate(modes):
            rb = QRadioButton(m)
            self.bg_mode.addButton(rb, i)
            layout_mode.addWidget(rb)
        grp_mode.setLayout(layout_mode)
        left_layout.addWidget(grp_mode)
        
        # Group 5: LEQ Settings
        grp_leq = QGroupBox("LEQ Settings")
        layout_leq = QHBoxLayout()
        layout_leq.addWidget(QLabel("Plot Interval:"))
        self.combo_leq_int = QComboBox()
        self.combo_leq_int.addItems(["100 ms", "1 sec", "10 sec", "1 min", "15 min", "1 hour"])
        layout_leq.addWidget(self.combo_leq_int)
        grp_leq.setLayout(layout_leq)
        left_layout.addWidget(grp_leq)
        
        left_layout.addStretch()
        
        # Action Buttons
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

    def _apply_settings_to_ui(self):
        """Applies loaded settings to UI widgets."""
        for btn in self.bg_weight.buttons():
            if btn.text() == self.settings.weighting:
                btn.setChecked(True)
                break
        else:
            self.bg_weight.button(0).setChecked(True)

        for btn in self.bg_speed.buttons():
            if btn.text() == self.settings.speed:
                btn.setChecked(True)
                break
        else:
            self.bg_speed.button(1).setChecked(True) 

        self.bg_mode.button(self.settings.analysis_mode_index).setChecked(True)
        self.combo_leq_int.setCurrentIndex(self.settings.leq_interval_index)
        self.cal_factor = self.settings.calibration_factor
        self._update_file_info_label()

    def _scrape_ui_to_settings(self):
        """Updates settings object from current UI state."""
        w_btn = self.bg_weight.checkedButton()
        if w_btn: self.settings.weighting = w_btn.text()
        
        s_btn = self.bg_speed.checkedButton()
        if s_btn: self.settings.speed = s_btn.text()
        
        self.settings.analysis_mode_index = self.bg_mode.checkedId()
        self.settings.leq_interval_index = self.combo_leq_int.currentIndex()
        self.settings.calibration_factor = self.cal_factor

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

    # --- Menu Actions ---

    def on_action_load_settings(self):
        fname, _ = QFileDialog.getOpenFileName(self, "Load Settings", self.settings.last_directory, "YAML Files (*.yaml);;All Files (*)")
        if fname:
            self.settings = self.settings_mgr.load(Path(fname))
            self._apply_settings_to_ui()
            self.status_bar.showMessage(f"Settings loaded from {Path(fname).name}")

    def on_action_save_settings(self):
        self._scrape_ui_to_settings()
        fname, _ = QFileDialog.getSaveFileName(self, "Save Settings", self.settings.last_directory, "YAML Files (*.yaml)")
        if fname:
            self.settings_mgr.save(self.settings, Path(fname))
            self.status_bar.showMessage(f"Settings saved to {Path(fname).name}")

    def on_action_save_figure(self):
        QMessageBox.information(self, "Not Implemented", "Save Plot functionality coming soon.")

    def on_open_url(self, url):
        QDesktopServices.openUrl(QUrl(url))

    def on_about(self):
        dlg = AboutDialog(self)
        dlg.exec()

    # --- Main Actions ---

    def on_load_file(self):
        start_dir = self.settings.last_directory
        fname, _ = QFileDialog.getOpenFileName(self, "Open WAV", start_dir, "WAV Files (*.wav)")
        
        if fname:
            path = Path(fname)
            self.filepath = path
            self.start_time = 0.0
            
            self.settings.last_directory = str(path.parent)
            
            try:
                from soundfile import info
                inf = info(str(self.filepath))
                self.end_time = inf.duration
                self._update_file_info_label(inf)
                
                self.btn_select.setEnabled(True)
                self.btn_analyze.setEnabled(True)
                self.menu_export.setEnabled(False)
                # Reset Unsaved Data flag when loading new file
                self.has_unsaved_data = False 
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
            self.settings.calibration_factor = self.cal_factor
            self._update_file_info_label()
            self.status_bar.showMessage(f"Calibration updated: {self.cal_factor:.4f}")

    def toggle_inputs(self, enabled: bool):
        self.btn_load.setEnabled(enabled)
        self.btn_select.setEnabled(enabled)
        self.btn_cal.setEnabled(enabled) 
        self.combo_leq_int.setEnabled(enabled)
        # Enable Export Menu if we have results and not running
        can_export = (len(self.last_results) > 0)
        self.menu_export.setEnabled(enabled and can_export)
        
        for child in self.left_panel.findChildren(QGroupBox):
             if child.title() != "File & Selection": 
                 child.setEnabled(enabled)

    def on_analyze_click(self):
        # 1. STOP Logic
        if self.worker is not None:
            self.status_bar.showMessage("Stopping...")
            self.worker.stop()
            self.btn_analyze.setEnabled(False)
            return

        # 2. START Logic
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
        
        self.toggle_inputs(False) # This handles disabling export menu
        self.btn_analyze.setText("STOP")
        self.btn_analyze.setStyleSheet("font-weight: bold; font-size: 14px; height: 40px; background-color: #fca5a5;")
        self.progress.setValue(0)
        self.status_bar.showMessage(f"Analyzing ({speed})...")
        
        self.worker = AnalysisWorker(
            filepath=self.filepath, 
            cal_factor=self.cal_factor, 
            block_size_ms=self.block_size_ms,
            weighting=weighting,
            do_bands=do_bands,
            band_res=res,
            speed=speed,
            band_order=self.settings.band_filter_order,
            ref_pressure=self.settings.ref_pressure
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
            
        self.last_results = filtered
        # Mark as UNSAVED since we have fresh results
        self.has_unsaved_data = True
        
        self.menu_export.setEnabled(True)

        self._plot_results(filtered, mode_id, weighting, speed)
        self.status_bar.showMessage("Analysis Complete.")

    def _plot_results(self, results: list, mode_id: int, weighting: str, speed: str):
        leq_int_txt = self.combo_leq_int.currentText()
        cur_std = self.settings.current_dose_standard
        dose_params = self.settings.dose_standards.get(cur_std)
        
        ResultPlotter.plot(
            self.plot_panel.figure,
            results,
            mode_id,
            weighting,
            speed,
            leq_int_txt,
            self.block_size_ms,
            dose_params=dose_params,
            ref_pressure=self.settings.ref_pressure
        )
        self.plot_panel.draw()

    def on_export_csv(self):
            if not self.last_results: return

            default_name = self.filepath.stem + "_results.csv" if self.filepath else "results.csv"
            path_str, _ = QFileDialog.getSaveFileName(self, "Export CSV", default_name, "CSV Files (*.csv)")
            if not path_str: return
            out_path = Path(path_str)
            
            w_btn = self.bg_weight.checkedButton()
            weighting = w_btn.text() if w_btn else 'A'
            s_btn = self.bg_speed.checkedButton()
            speed = s_btn.text() if s_btn else 'Fast'
            mode_id = self.bg_mode.checkedId()
            
            # Retrieve Settings
            cur_std = self.settings.current_dose_standard
            dose_params = self.settings.dose_standards.get(cur_std)
            ref_pressure = self.settings.ref_pressure

            try:
                match mode_id:
                    case 1: # LEQ
                        interval_txt = self.combo_leq_int.currentText()
                        ResultsExporter.export_leq(
                            out_path, 
                            self.last_results, 
                            self.block_size_ms, 
                            interval_txt, 
                            weighting,
                            dose_params,   # Passed
                            ref_pressure   # Passed
                        )
                    case 0: # Lp
                        ResultsExporter.export_lp(out_path, self.last_results, weighting, speed)
                    case 2 | 3: # Spectrum
                        ResultsExporter.export_spectrum(
                            out_path, 
                            self.last_results, 
                            weighting,
                            ref_pressure   # Passed
                        )
                
                # Export Successful: Mark as SAVED
                self.has_unsaved_data = False
                
                self.status_bar.showMessage(f"Exported to {out_path.name}")
                QMessageBox.information(self, "Export Successful", f"Data saved to:\n{out_path}")
            except Exception as e:
                QMessageBox.critical(self, "Export Failed", str(e))

    def closeEvent(self, event):
        """Save settings and check for unsaved data on exit."""
        # 1. Thread Safety Check
        if self.worker is not None and self.worker.isRunning():
            reply = QMessageBox.question(
                self, 'Analysis Running',
                "An analysis is currently running.\nDo you want to stop it and exit?",
                QMessageBox.Yes | QMessageBox.No, QMessageBox.No
            )
            if reply == QMessageBox.Yes:
                self.worker.stop()
                self.worker.wait()
            else:
                event.ignore()
                return

        # 2. Unsaved Data Check
        if self.has_unsaved_data:
            reply = QMessageBox.question(
                self, 'Unsaved Results',
                "There may be unsaved analysis results.\nAre you sure you want to quit?",
                QMessageBox.Yes | QMessageBox.No, QMessageBox.No
            )
            if reply == QMessageBox.No:
                event.ignore()
                return

        # 3. Save Settings
        self._scrape_ui_to_settings()
        self.settings_mgr.save(self.settings) 
        
        event.accept()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    w = MainWindow()
    w.show()
    sys.exit(app.exec())