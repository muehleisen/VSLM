import sys
from pathlib import Path
from PySide6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                               QHBoxLayout, QPushButton, QLabel, QGroupBox, 
                               QFileDialog, QMessageBox, QFrame, QButtonGroup, 
                               QRadioButton, QProgressBar, QComboBox)
from PySide6.QtGui import QAction, QDesktopServices
from PySide6.QtCore import QUrl

# --- VSLM Imports ---
from .waveform_dialog import WaveformDialog
from .calibration_dialog import CalibrationDialog
from .about_dialog import AboutDialog 
from .plot_widget import MatplotlibWidget
from .analysis_worker import AnalysisWorker
from .plot_manager import ResultPlotter
from ..result_exporter import ResultsExporter
from ..settings_manager import SettingsManager, AppSettings
from ..constants import LEQ_INTERVAL_MAP # New Import


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("VSLM 2.0 (Python)")
        self.resize(1150, 800)
        
        self.settings_mgr = SettingsManager()
        self.settings = self.settings_mgr.load() 
        
        self.filepath: Path | None = None
        self.start_time: float = 0.0
        self.end_time: float | None = None
        self.cal_factor: float = self.settings.calibration_factor
        self.block_size_ms: float = self.settings.block_size_ms
        
        self.last_results: list = [] 
        self.has_unsaved_data: bool = False 
        self.worker: AnalysisWorker | None = None
        
        self._init_menu_bar()
        self._init_ui()
        self._apply_settings_to_ui()
        
        self.status_bar = self.statusBar()
        self.status_bar.showMessage("Ready. Load a file to begin.")

    def _init_menu_bar(self):
        menu_bar = self.menuBar()
        menu_file = menu_bar.addMenu("File")
        
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
        act_quit.triggered.connect(self.close) 
        menu_file.addAction(act_quit)
        
        self.menu_export = menu_bar.addMenu("Export")
        self.menu_export.setEnabled(False) 
        act_export_csv = QAction("Save Results (CSV)...", self)
        act_export_csv.triggered.connect(self.on_export_csv)
        self.menu_export.addAction(act_export_csv)
        
        act_save_fig = QAction("Save Plot Figure...", self)
        act_save_fig.triggered.connect(self.on_action_save_figure)
        self.menu_export.addAction(act_save_fig)
        
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
        
        left_panel = QWidget()
        left_panel.setFixedWidth(320)
        left_layout = QVBoxLayout(left_panel)
        self.left_panel = left_panel 
        
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
        
        self.lbl_info = QLabel("No File Loaded")
        self.lbl_info.setFrameStyle(QFrame.StyledPanel | QFrame.Sunken)
        
        layout_file.addWidget(self.btn_load)
        layout_file.addWidget(self.btn_select)
        layout_file.addWidget(self.btn_cal) 
        layout_file.addWidget(QLabel("File Info"))
        layout_file.addWidget(self.lbl_info)
        grp_file.setLayout(layout_file)
        left_layout.addWidget(grp_file)
        
        grp_weight = QGroupBox("Weighting")
        layout_weight = QHBoxLayout()
        self.bg_weight = QButtonGroup()
        for i, text in enumerate(['A', 'C', 'Z']):
            rb = QRadioButton(text)
            self.bg_weight.addButton(rb, i)
            layout_weight.addWidget(rb)
        grp_weight.setLayout(layout_weight)
        left_layout.addWidget(grp_weight)
        
        grp_speed = QGroupBox("Speed (Lp Mode)")
        layout_speed = QHBoxLayout()
        self.bg_speed = QButtonGroup()
        for i, text in enumerate(['Slow', 'Fast', 'Impulse']):
            rb = QRadioButton(text)
            self.bg_speed.addButton(rb, i)
            layout_speed.addWidget(rb)
        grp_speed.setLayout(layout_speed)
        left_layout.addWidget(grp_speed)
        
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
        
        grp_leq = QGroupBox("LEQ Settings")
        layout_leq = QHBoxLayout()
        layout_leq.addWidget(QLabel("Plot Interval:"))
        self.combo_leq_int = QComboBox()
        
        # --- REFACTOR: Use Map to Populate ---
        for key, (label, _) in LEQ_INTERVAL_MAP.items():
            self.combo_leq_int.addItem(label, key) # Store Enum Key as UserData
            
        layout_leq.addWidget(self.combo_leq_int)
        grp_leq.setLayout(layout_leq)
        left_layout.addWidget(grp_leq)
        
        left_layout.addStretch()
        
        self.progress = QProgressBar()
        self.progress.setTextVisible(False)
        left_layout.addWidget(self.progress)
        
        self.btn_analyze = QPushButton("ANALYZE")
        self.btn_analyze.setStyleSheet("font-weight: bold; font-size: 14px; height: 40px; background-color: #dbeafe;")
        self.btn_analyze.clicked.connect(self.on_analyze_click)
        self.btn_analyze.setEnabled(False)
        left_layout.addWidget(self.btn_analyze)
        
        main_layout.addWidget(left_panel)
        
        self.plot_panel = MatplotlibWidget()
        self.plot_panel.sig_scaling_changed.connect(self.on_scaling_changed)
        main_layout.addWidget(self.plot_panel, stretch=1)

    def _apply_settings_to_ui(self):
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
        
        self.plot_panel.set_plot_settings(
            self.settings.plot_autoscale,
            self.settings.plot_ymin,
            self.settings.plot_ymax
        )
        self._update_file_info_label()

    def _scrape_ui_to_settings(self):
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
                              f"Dur: {inf.duration:.1f} s\n"
                              f"Cal Factor: {self.cal_factor:.4f}")

    def on_scaling_changed(self, auto, ymin, ymax):
        self.settings.plot_autoscale = auto
        self.settings.plot_ymin = ymin
        self.settings.plot_ymax = ymax
        if self.last_results:
            self._redraw_plot()

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
        if self.plot_panel.toolbar:
            self.plot_panel.toolbar.save_figure()
        else:
            QMessageBox.information(self, "Info", "Use the floppy disk icon on the plot to save.")

    def on_open_url(self, url):
        QDesktopServices.openUrl(QUrl(url))

    def on_about(self):
        dlg = AboutDialog(self)
        dlg.exec()

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
            self._update_file_info_label()

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
        can_export = (len(self.last_results) > 0)
        self.menu_export.setEnabled(enabled and can_export)
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
        self.btn_analyze.setStyleSheet("background-color: #fca5a5;")
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
        self.btn_analyze.setStyleSheet("background-color: #dbeafe;")
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
        self.has_unsaved_data = True
        self.menu_export.setEnabled(True)

        self._plot_results(filtered, mode_id, weighting, speed)
        self.status_bar.showMessage("Analysis Complete.")

    def _plot_results(self, results: list, mode_id: int, weighting: str, speed: str):
        # --- REFACTOR: Retrieve Enum Key ---
        leq_int_key = self.combo_leq_int.currentData()
        
        cur_std = self.settings.current_dose_standard
        dose_params = self.settings.dose_standards.get(cur_std)
        
        ResultPlotter.plot(
            self.plot_panel.figure,
            results,
            mode_id,
            weighting,
            speed,
            leq_int_key, # Pass Key
            self.block_size_ms,
            dose_params=dose_params,
            ref_pressure=self.settings.ref_pressure,
            autoscale=self.settings.plot_autoscale,
            ymin=self.settings.plot_ymin,
            ymax=self.settings.plot_ymax
        )
        self.plot_panel.draw()

    def _redraw_plot(self):
        if not self.last_results: return
        w_btn = self.bg_weight.checkedButton()
        weighting = w_btn.text() if w_btn else 'A'
        s_btn = self.bg_speed.checkedButton()
        speed = s_btn.text() if s_btn else 'Fast'
        mode_id = self.bg_mode.checkedId()
        self._plot_results(self.last_results, mode_id, weighting, speed)

    def on_export_csv(self):
        if not self.last_results: return
        path_str, _ = QFileDialog.getSaveFileName(self, "Export CSV", "results.csv", "CSV Files (*.csv)")
        if not path_str: return
        
        out_path = Path(path_str)
        w_btn = self.bg_weight.checkedButton()
        weighting = w_btn.text() if w_btn else 'A'
        s_btn = self.bg_speed.checkedButton()
        speed = s_btn.text() if s_btn else 'Fast'
        mode_id = self.bg_mode.checkedId()
        
        cur_std = self.settings.current_dose_standard
        dose_params = self.settings.dose_standards.get(cur_std)
        ref_pressure = self.settings.ref_pressure

        try:
            match mode_id:
                case 1: # LEQ
                    # --- REFACTOR: Retrieve Enum Key ---
                    interval_key = self.combo_leq_int.currentData()
                    
                    ResultsExporter.export_leq(
                        out_path, 
                        self.last_results, 
                        self.block_size_ms, 
                        interval_key, # Pass Key
                        weighting,
                        dose_params,
                        ref_pressure
                    )
                case 0: # Lp
                    ResultsExporter.export_lp(out_path, self.last_results, weighting, speed)
                case 2 | 3: # Spectrum
                    ResultsExporter.export_spectrum(
                        out_path, 
                        self.last_results, 
                        weighting,
                        ref_pressure
                    )
            
            self.has_unsaved_data = False
            self.status_bar.showMessage(f"Exported to {out_path.name}")
            QMessageBox.information(self, "Export Successful", f"Data saved to:\n{out_path}")
        except Exception as e:
            QMessageBox.critical(self, "Export Failed", str(e))

    def closeEvent(self, event):
        if self.worker is not None and self.worker.isRunning():
            self.worker.stop()
            self.worker.wait()
        self._scrape_ui_to_settings()
        self.settings_mgr.save(self.settings) 
        event.accept()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    w = MainWindow()
    w.show()
    sys.exit(app.exec())