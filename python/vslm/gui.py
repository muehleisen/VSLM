# vslm/gui.py
import sys
import os
from typing import Tuple, Any

import numpy as np
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, 
    QPushButton, QLabel, QGroupBox, QButtonGroup, 
    QFileDialog, QMessageBox, QFrame, QProgressBar
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QAction

# Backend Imports
from .core import VSLMCore
from .playback import AudioPlayer
from . import analysis

# Refactored Modules
from . import dialogs
from . import widgets
from . import workers
from .settings import VSLMSettings

# --- Constants ---

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

class MainWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        
        self.setWindowTitle("VSLM - Virtual Sound Level Meter (Python Port)")
        self.resize(1050, 700)
        
        # Initialize Backend
        self.core = VSLMCore()
        self.player = AudioPlayer()
        self.worker = None 
        
        # Cache for dynamic updates
        self.cached_result: Tuple[str, Any] = None 
        
        # Initialize Settings Dataclass
        self.settings = VSLMSettings()

        # Main UI Container
        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)
        self.main_layout = QHBoxLayout(self.central_widget)
        
        # UI Elements
        self.info_label: QLabel
        self.btn_analyze: QPushButton
        self.btn_play: QPushButton
        self.wtg_bg: QButtonGroup
        self.spd_bg: QButtonGroup
        self.mode_bg: QButtonGroup
        self.progress_bar: QProgressBar
        self.plot_widget: widgets.VSLMPlotWidget
        
        self._create_menubar()
        self._create_left_panel()
        self._create_right_panel()
        
        self.status_bar = self.statusBar()
        self.status_bar.showMessage("Ready. Please load a measurement file.")

    # --- UI Construction ---

    def _create_menubar(self) -> None:
        menubar = self.menuBar()

        # File Menu
        file_menu = menubar.addMenu("File")
        act_load_settings = QAction("Load Settings", self)
        act_load_settings.triggered.connect(lambda: self._placeholder_action("Load Settings"))
        file_menu.addAction(act_load_settings)
        act_save_settings = QAction("Save Settings", self)
        act_save_settings.triggered.connect(lambda: self._placeholder_action("Save Settings"))
        file_menu.addAction(act_save_settings)
        file_menu.addSeparator()
        act_quit = QAction("Quit", self)
        act_quit.triggered.connect(self.confirm_quit) 
        file_menu.addAction(act_quit)

        # Lpplot Menu
        lpplot_menu = menubar.addMenu("Lpplot")
        act_plot_spacing = QAction("Set Plot Time Spacing", self)
        act_plot_spacing.triggered.connect(self.dlg_plot_spacing)
        lpplot_menu.addAction(act_plot_spacing)
        act_plot_scales = QAction("Plot Scales", self)
        act_plot_scales.triggered.connect(self.dlg_lpplot_scales)
        lpplot_menu.addAction(act_plot_scales)

        # Leq Menu
        leq_menu = menubar.addMenu("Leq")
        act_leq_int = QAction("LEQ Integration Time", self)
        act_leq_int.triggered.connect(self.dlg_leq_integration)
        leq_menu.addAction(act_leq_int)
        act_leq_perc = QAction("Set LEQ Percentile", self)
        act_leq_perc.triggered.connect(self.dlg_leq_percentile)
        leq_menu.addAction(act_leq_perc)
        act_leq_scales = QAction("Plot Scales", self)
        act_leq_scales.triggered.connect(self.dlg_leq_scales)
        leq_menu.addAction(act_leq_scales)
        act_leq_dose = QAction("Noise Dose Criterion", self)
        act_leq_dose.triggered.connect(self.dlg_noise_dose)
        leq_menu.addAction(act_leq_dose)

        # Band Menu
        band_menu = menubar.addMenu("Band")
        act_resolution = QAction("Resolution", self)
        act_resolution.triggered.connect(self.dlg_band_resolution)
        band_menu.addAction(act_resolution)
        act_method = QAction("Method", self)
        act_method.triggered.connect(self.dlg_band_method)
        band_menu.addAction(act_method)

        # PSD Menu
        psd_menu = menubar.addMenu("PSD")
        act_psd_fft = QAction("FFT Size", self)
        act_psd_fft.triggered.connect(self.dlg_psd_fft)
        psd_menu.addAction(act_psd_fft)
        act_psd_overlap = QAction("Overlap", self)
        act_psd_overlap.triggered.connect(self.dlg_psd_overlap)
        psd_menu.addAction(act_psd_overlap)
        act_psd_window = QAction("Window", self)
        act_psd_window.triggered.connect(self.dlg_psd_window)
        psd_menu.addAction(act_psd_window)
        
        # New PSD Scales option
        act_psd_scales = QAction("Plot Scales", self)
        act_psd_scales.triggered.connect(self.dlg_psd_scales)
        psd_menu.addAction(act_psd_scales)

        # Spectrogram Menu
        spec_menu = menubar.addMenu("Spectrogram")
        act_spec_fft = QAction("FFT Size", self)
        act_spec_fft.triggered.connect(self.dlg_spec_fft)
        spec_menu.addAction(act_spec_fft)
        
        act_spec_overlap = QAction("Overlap", self)
        act_spec_overlap.triggered.connect(self.dlg_spec_overlap)
        spec_menu.addAction(act_spec_overlap)
        
        act_spec_slice = QAction("Slice Length", self)
        act_spec_slice.triggered.connect(self.dlg_spec_slice)
        spec_menu.addAction(act_spec_slice)
        
        act_spec_scale = QAction("Plot Scale", self)
        act_spec_scale.triggered.connect(self.dlg_spec_scales)
        spec_menu.addAction(act_spec_scale)
        act_spec_view = QAction("View", self)
        act_spec_view.triggered.connect(self.dlg_spec_view)
        spec_menu.addAction(act_spec_view)

        # Help Menu
        help_menu = menubar.addMenu("Help")
        act_about = QAction("About VSLM", self)
        act_about.triggered.connect(lambda: dialogs.AboutDialog.show(self))
        help_menu.addAction(act_about)

    def _create_left_panel(self) -> None:
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
        right_widget = QWidget()
        layout = QVBoxLayout(right_widget)
        
        self.plot_widget = widgets.VSLMPlotWidget()
        self.plot_widget.show_placeholder()
        
        self.figure = self.plot_widget.figure
        self.canvas = self.plot_widget.canvas
        self.ax = self.plot_widget.ax
        
        layout.addWidget(self.plot_widget)
        self.main_layout.addWidget(right_widget, stretch=1)

    # --- Group Creators ---

    def _create_io_group(self) -> QGroupBox:
        group = QGroupBox("File & Calibration")
        layout = QVBoxLayout()
        
        btn_load = QPushButton("Load Measurement (.wav)")
        btn_load.clicked.connect(self.load_measurement)
        
        btn_cal = QPushButton("Set Calibration")
        btn_cal.clicked.connect(self.dlg_calibration)
        
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
        
        map_w = {'A': 1, 'C': 2, 'Z': 3}
        default_id = map_w.get(self.settings.frequency_weighting, 1)

        for text, uid in options:
            pair = widgets.LabelledToggleButton(text, uid, self.wtg_bg)
            if uid == default_id: 
                pair.setChecked(True)
            layout.addWidget(pair)
            
        group.setLayout(layout)
        return group

    def _create_speed_group(self) -> QGroupBox:
        group = QGroupBox("Meter Speed")
        layout = QHBoxLayout()
        layout.setSpacing(10)
        self.spd_bg = QButtonGroup(self)
        options = [("Slow\n(1.0s)", 1), ("Fast\n(125ms)", 2), ("Impulse\n(35ms/1.5s)", 3)]
        
        map_s = {'Slow': 1, 'Fast': 2, 'Impulse': 3}
        default_id = map_s.get(self.settings.meter_speed, 1)

        for text, uid in options:
            pair = widgets.LabelledToggleButton(text, uid, self.spd_bg)
            if uid == default_id:
                pair.setChecked(True)
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
            from PySide6.QtWidgets import QRadioButton
            rb = QRadioButton(name)
            rb.setProperty("tag", tag)
            if tag == self.settings.analysis_mode:
                rb.setChecked(True)
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

    # --- Slots (Dialogs) ---

    def dlg_calibration(self):
        if self.core.audio_data is None:
            QMessageBox.warning(self, "Warning", "Please load a measurement file first.")
            return
        dlg = dialogs.CalibrationDialog(current_val=94.0, parent=self)
        if dlg.exec():
            try:
                factor = self.core.set_calibration(dlg.get_value())
                self._update_info_text()
                QMessageBox.information(self, "Calibration", f"Factor set to {factor:.4f}")
            except Exception as e:
                QMessageBox.critical(self, "Error", str(e))

    def dlg_plot_spacing(self):
        dlg = dialogs.PlotTimeSpacingDialog(self.settings.plot_time_spacing, self)
        if dlg.exec(): 
            self.settings.plot_time_spacing = dlg.get_value()

    def dlg_lpplot_scales(self):
        cmin, cmax = self.settings.lpplot_y_min, self.settings.lpplot_y_max
        auto = self.settings.lpplot_autoscale
        
        dlg = dialogs.PlotScalesDialog(cmin, cmax, auto, self)
        if dlg.exec(): 
            self.settings.lpplot_y_min, self.settings.lpplot_y_max, self.settings.lpplot_autoscale = dlg.get_values()
            # Dynamic Update for Lp
            if self.cached_result and self.cached_result[0] == 'lp':
                self.handle_analysis_result(self.cached_result)

    def dlg_leq_integration(self):
        dlg = dialogs.LeqIntegrationDialog(self.settings.leq_integration_time, self)
        if dlg.exec(): 
            self.settings.leq_integration_time = dlg.get_value()

    def dlg_leq_percentile(self):
        dlg = dialogs.LeqPercentileDialog(self.settings.leq_percentile, self)
        if dlg.exec(): 
            self.settings.leq_percentile = dlg.get_value()

    def dlg_leq_scales(self):
        cmin, cmax = self.settings.leq_y_min, self.settings.leq_y_max
        auto = self.settings.leq_autoscale
        
        dlg = dialogs.PlotScalesDialog(cmin, cmax, auto, self)
        if dlg.exec(): 
            self.settings.leq_y_min, self.settings.leq_y_max, self.settings.leq_autoscale = dlg.get_values()
            if self.cached_result and self.cached_result[0] == 'leq':
                self.handle_analysis_result(self.cached_result)

    def dlg_noise_dose(self):
        e = self.settings.dose_exchange_rate
        t = self.settings.dose_threshold
        c = self.settings.dose_criterion
        dlg = dialogs.NoiseDoseDialog(e, t, c, self)
        if dlg.exec(): 
            self.settings.dose_exchange_rate, self.settings.dose_threshold, self.settings.dose_criterion = dlg.get_values()

    def dlg_band_resolution(self):
        dlg = dialogs.BandResolutionDialog(self.settings.band_resolution_index, self)
        if dlg.exec(): pass 

    def dlg_band_method(self):
        dlg = dialogs.BandMethodDialog(self.settings.band_method_index, self)
        if dlg.exec(): pass

    def dlg_psd_fft(self):
        dlg = dialogs.FFTSizeDialog(self.settings.psd_fft_size, self)
        if dlg.exec(): 
            self.settings.psd_fft_size = dlg.get_value()

    def dlg_psd_overlap(self):
        dlg = dialogs.OverlapDialog(self.settings.psd_overlap_percent, self)
        if dlg.exec(): 
            self.settings.psd_overlap_percent = dlg.get_value()

    def dlg_psd_window(self):
        dlg = dialogs.WindowDialog(self.settings.psd_window, self)
        if dlg.exec(): 
            self.settings.psd_window = dlg.get_value()

    def dlg_psd_scales(self):
        cmin, cmax = self.settings.psd_y_min, self.settings.psd_y_max
        auto = self.settings.psd_autoscale
        
        dlg = dialogs.PlotScalesDialog(cmin, cmax, auto, self)
        if dlg.exec(): 
            self.settings.psd_y_min, self.settings.psd_y_max, self.settings.psd_autoscale = dlg.get_values()
            # Dynamic Update
            if self.cached_result and self.cached_result[0] == 'psd':
                self.handle_analysis_result(self.cached_result)

    def dlg_spec_fft(self):
        dlg = dialogs.FFTSizeDialog(self.settings.spec_fft_size, self)
        if dlg.exec(): 
            self.settings.spec_fft_size = dlg.get_value()

    def dlg_spec_overlap(self):
        """Dialog for Spectrogram Overlap %."""
        dlg = dialogs.OverlapDialog(self.settings.spec_overlap_percent, self)
        if dlg.exec():
            self.settings.spec_overlap_percent = dlg.get_value()
            self.settings.spec_use_overlap = True

    def dlg_spec_slice(self):
        dlg = dialogs.SliceLengthDialog(self.settings.spec_slice_length, self)
        if dlg.exec(): 
            self.settings.spec_slice_length = dlg.get_value()
            self.settings.spec_use_overlap = False

    def dlg_spec_scales(self):
        cmin, cmax = self.settings.spec_y_min, self.settings.spec_y_max
        auto = self.settings.spec_autoscale
        
        dlg = dialogs.PlotScalesDialog(cmin, cmax, auto, self)
        if dlg.exec(): 
            self.settings.spec_y_min, self.settings.spec_y_max, self.settings.spec_autoscale = dlg.get_values()
            # Dynamic Update
            if self.cached_result and self.cached_result[0] == 'spec':
                self.handle_analysis_result(self.cached_result)

    def dlg_spec_view(self):
        c = self.settings.spec_colormap
        dlg = dialogs.SpectrogramViewDialog(c, self)
        if dlg.exec(): 
            self.settings.spec_colormap = dlg.get_value()
            if self.cached_result and self.cached_result[0] == 'spec':
                self.status_bar.showMessage(f"Updating colormap to {self.settings.spec_colormap}...")
                self.handle_analysis_result(self.cached_result)

    # --- Main Logic ---

    def _placeholder_action(self, name: str) -> None:
        self.status_bar.showMessage(f"Menu action '{name}' triggered (Not implemented).")

    def confirm_quit(self) -> None:
        reply = QMessageBox.question(self, "Confirm Quit", "Are you sure you want to quit?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No, QMessageBox.StandardButton.No)
        if reply == QMessageBox.StandardButton.Yes: QApplication.instance().quit()

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

    def toggle_playback(self) -> None:
        if self.player.is_playing:
            self.player.stop()
            self.btn_play.setText("Play Audio")
        else:
            if self.core.audio_data is not None:
                self.player.play(self.core.audio_data, self.core.fs)
                self.btn_play.setText("Stop Audio")

    def get_selected_weighting(self) -> str:
        return {1: 'A', 2: 'C', 3: 'Z'}.get(self.wtg_bg.checkedId(), 'A')

    def get_selected_speed(self) -> str:
        return {1: 'Slow', 2: 'Fast', 3: 'Impulse'}.get(self.spd_bg.checkedId(), 'Slow')

    def start_analysis(self) -> None:
        mode_btn = self.mode_bg.checkedButton()
        if not mode_btn: return
        
        mode_tag = mode_btn.property("tag")
        weighting = self.get_selected_weighting()
        speed = self.get_selected_speed()
        
        self.settings.analysis_mode = mode_tag
        self.settings.frequency_weighting = weighting
        self.settings.meter_speed = speed
        
        self.status_bar.showMessage(f"Analyzing {mode_tag.upper()}... please wait.")
        self._set_ui_busy(True)
        
        if mode_tag == 'lp':
            self.worker = workers.AnalysisWorker(self.core.calculate_lp, weighting=weighting, speed=speed)
        elif mode_tag == 'leq':
            self.worker = workers.AnalysisWorker(self.core.calculate_leq, weighting=weighting)
        elif mode_tag == 'octave':
            self.worker = workers.AnalysisWorker(analysis.calculate_ansi_bands, self.core, weighting=weighting, resolution='octave')
        elif mode_tag == 'third':
            self.worker = workers.AnalysisWorker(analysis.calculate_ansi_bands, self.core, weighting=weighting, resolution='third')
        elif mode_tag == 'psd':
            self.worker = workers.AnalysisWorker(
                analysis.calculate_psd, 
                self.core,
                nfft=self.settings.psd_fft_size,
                overlap_percent=self.settings.psd_overlap_percent,
                window=self.settings.psd_window
            )
        elif mode_tag == 'spec':
            # Logic: If using overlap, force slice_len to None so backend uses overlap ratio
            s_len = self.settings.spec_slice_length
            if self.settings.spec_use_overlap:
                s_len = None
            
            self.worker = workers.AnalysisWorker(
                analysis.calculate_spectrogram, 
                self.core, 
                weighting=weighting,
                nfft=self.settings.spec_fft_size,
                slice_len=s_len,
                overlap_ratio=self.settings.spec_overlap_percent / 100.0
            )

        if self.worker:
            self.worker.result_ready.connect(lambda res: self.handle_analysis_result((mode_tag, res)))
            self.worker.error_occurred.connect(self.handle_analysis_error)
            self.worker.progress_updated.connect(self.progress_bar.setValue)
            self.worker.start()

    def _set_ui_busy(self, busy: bool) -> None:
        self.btn_analyze.setEnabled(not busy)
        self.central_widget.setEnabled(not busy)
        self.progress_bar.setVisible(busy)
        if busy: self.progress_bar.setValue(0)

    def handle_analysis_error(self, msg: str) -> None:
        self._set_ui_busy(False)
        self.status_bar.showMessage("Analysis Failed.")
        QMessageBox.critical(self, "Analysis Error", msg)

    def handle_analysis_result(self, payload: Tuple[str, Any]) -> None:
        mode, result = payload
        
        self.cached_result = payload
        
        self._set_ui_busy(False)
        self.status_bar.showMessage("Analysis Complete.")
        
        weighting = self.settings.frequency_weighting
        speed = self.settings.meter_speed

        if mode == 'lp':
            self.ax = self.plot_widget.prepare_plot()
            t, lp = result
            self.ax.plot(t, lp)
            self.ax.set_title(f"Sound Pressure Level ({weighting}-Weighted, {speed})")
            self.ax.set_xlabel("Time (s)")
            self.ax.set_ylabel("Lp (dB)")
            self.ax.grid(True)
            
            if self.settings.lpplot_autoscale and len(lp) > 0:
                ymin, ymax = np.min(lp) - 5, np.max(lp) + 5
                self.ax.set_ylim(ymin, ymax)
            else:
                self.ax.set_ylim(self.settings.lpplot_y_min, self.settings.lpplot_y_max)
            
        elif mode == 'leq':
            self.ax = self.plot_widget.prepare_plot()
            leq_val = result
            self.ax.axis('off')
            self.ax.text(0.5, 0.6, f"Leq ({weighting})", ha='center', fontsize=16)
            self.ax.text(0.5, 0.4, f"{leq_val:.2f} dB", ha='center', fontsize=30, fontweight='bold', color='blue')
            
        elif mode in ['octave', 'third']:
            self.ax = self.plot_widget.prepare_plot()
            freqs, levels = result
            x_pos = np.arange(len(freqs))
            self.ax.bar(x_pos, levels, width=0.8, color='green', alpha=0.7)
            self.ax.set_xticks(x_pos)
            labels = [f"{int(f)}" if f < 1000 else f"{f/1000:.1f}k" for f in freqs]
            if mode == 'third':
                labels = [lbl if i % 3 == 0 else "" for i, lbl in enumerate(labels)]
            self.ax.set_xticklabels(labels, rotation=45)
            
            # --- Updated Title with Weighting ---
            mode_text = 'Octave' if mode == 'octave' else '1/3 Octave'
            self.ax.set_title(f"{mode_text} Band Levels ({weighting}-Weighted)")
            
            self.ax.set_ylabel("dB")
            self.ax.grid(axis='y')

        elif mode == 'psd':
            self.ax = self.plot_widget.prepare_plot()
            f, lpxx = result
            self.ax.semilogx(f, lpxx)
            self.ax.set_title("Power Spectral Density")
            self.ax.set_xlabel("Frequency (Hz)")
            self.ax.set_ylabel("dB/Hz")
            self.ax.grid(True, which="both")
            self.ax.set_xlim(20, self.core.fs/2)
            
            if self.settings.psd_autoscale and len(lpxx) > 0:
                ymin, ymax = np.min(lpxx) - 5, np.max(lpxx) + 5
                self.ax.set_ylim(ymin, ymax)
            else:
                self.ax.set_ylim(self.settings.psd_y_min, self.settings.psd_y_max)
            
        elif mode == 'spec':
            f, t, Sxx_db = result
            
            cmap_name = self.settings.spec_colormap
            
            if self.settings.spec_autoscale:
                vmin, vmax = np.min(Sxx_db), np.max(Sxx_db)
            else:
                vmin = self.settings.spec_y_min
                vmax = self.settings.spec_y_max

            self.ax = self.plot_widget.prepare_plot()
            im = self.ax.pcolormesh(
                t, f, Sxx_db, 
                shading='auto', 
                cmap=cmap_name,
                vmin=vmin, 
                vmax=vmax
            )
            self.figure.colorbar(im, ax=self.ax).set_label('Intensity (dB)')
            
            self.ax.set_title(f"Spectrogram ({weighting}-Weighted)")
            self.ax.set_xlabel("Time (s)")
            self.ax.set_ylabel("Frequency (Hz)")

        self.plot_widget.draw()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())