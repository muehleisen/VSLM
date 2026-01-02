import sys
import numpy as np
from pathlib import Path
from PySide6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                               QHBoxLayout, QPushButton, QLabel, QGroupBox, 
                               QFileDialog, QMessageBox, QFrame, QButtonGroup, 
                               QRadioButton, QProgressBar, QComboBox)
from PySide6.QtCore import Qt

# Imports
from .. import leq
from .waveform import WaveformDialog
from .widgets import MatplotlibWidget
from .workers import AnalysisWorker

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
        
        # Worker Reference
        self.worker: AnalysisWorker | None = None
        
        self._init_ui()
        self.status_bar = self.statusBar()
        self.status_bar.showMessage("Ready. Load a file to begin.")

    def _init_ui(self):
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
        self.lbl_info_header = QLabel("File Info")
        self.lbl_info_header.setStyleSheet("font-weight: bold; font-size: 10px; margin-top: 5px;")
        self.lbl_info = QLabel("No File Loaded")
        self.lbl_info.setFrameStyle(QFrame.StyledPanel | QFrame.Sunken)
        self.lbl_info.setStyleSheet("padding: 3px;")
        layout_file.addWidget(self.btn_load)
        layout_file.addWidget(self.btn_select)
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
        self.combo_leq_int.addItems(["1 sec", "10 sec", "1 min", "15 min", "1 hour"])
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

    def on_load_file(self):
        fname, _ = QFileDialog.getOpenFileName(self, "Open WAV", "", "WAV Files (*.wav)")
        if fname:
            self.filepath = Path(fname)
            self.start_time = 0.0
            try:
                from soundfile import info
                inf = info(str(self.filepath))
                self.end_time = inf.duration
                self.lbl_info.setText(f"File: {self.filepath.name}\n"
                                      f"Fs: {inf.samplerate} Hz\n"
                                      f"Dur: {inf.duration:.1f} s   Block Size: {self.block_size_ms} ms")
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

    def toggle_inputs(self, enabled: bool):
        self.btn_load.setEnabled(enabled)
        self.btn_select.setEnabled(enabled)
        self.combo_leq_int.setEnabled(enabled)
        for child in self.left_panel.findChildren(QGroupBox):
             if child.title() != "File & Selection": 
                 child.setEnabled(enabled)

    def on_analyze_click(self):
        # 1. Handle Stop Request
        if self.worker is not None:
            self.status_bar.showMessage("Stopping...")
            self.worker.stop()
            self.btn_analyze.setEnabled(False)
            return

        # 2. Handle Start Request
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
        
        # Clean up C++ resources immediately upon thread completion
        self.worker.finished.connect(self.worker.deleteLater)
        
        self.worker.start()

    def on_worker_stopped(self):
        """Cleanup after thread exit."""
        self.worker = None
        self.toggle_inputs(True)
        self.btn_analyze.setText("ANALYZE")
        self.btn_analyze.setStyleSheet("font-weight: bold; font-size: 14px; height: 40px; background-color: #dbeafe;")
        self.btn_analyze.setEnabled(True)
        self.progress.setValue(0) # Clear progress bar

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
        fig = self.plot_panel.figure
        fig.clear()
        
        if not results: 
            self.plot_panel.draw()
            return

        match mode_id:
            case 1: # LEQ MODE
                int_txt = self.combo_leq_int.currentText()
                match int_txt:
                    case "1 sec": interval = 1.0
                    case "10 sec": interval = 10.0
                    case "1 min": interval = 60.0
                    case "15 min": interval = 900.0
                    case "1 hour": interval = 3600.0
                    case _: interval = 1.0
                
                stats = leq.calculate_leq_analysis(results, self.block_size_ms, interval)
                
                ax1 = fig.add_subplot(2, 1, 1)
                if len(stats.history['time']) > 0:
                    t_plot = list(stats.history['time'])
                    t_plot.append(t_plot[-1] + interval)
                    l_plot = list(stats.history['leq'])
                    l_plot.append(l_plot[-1])
                    ax1.step(t_plot, l_plot, where='post', color='b', linewidth=1.5)
                
                ax1.set_title(f"LEQ History ({int_txt} interval, {weighting}-weighted)")
                ax1.set_ylabel("LEQ (dB)")
                ax1.grid(True)
                
                ax2 = fig.add_subplot(2, 1, 2)
                ax2.axis('off')
                col1, col2, col3 = 0.05, 0.35, 0.65
                
                ax2.text(0.5, 0.95, f"Overall LEQ: {stats.overall:.1f} dB", 
                         ha='center', fontsize=14, fontweight='bold', color='blue')
                ax2.text(col1, 0.80, f"Lmax: {stats.max:.1f} dB")
                ax2.text(col1, 0.65, f"Lmin: {stats.min:.1f} dB")
                ax2.text(col1, 0.50, f"L10: {stats.ln[10]:.1f} dB")
                ax2.text(col1, 0.35, f"L50: {stats.ln[50]:.1f} dB")
                ax2.text(col1, 0.20, f"L90: {stats.ln[90]:.1f} dB")
                ax2.text(col2, 0.80, f"L20: {stats.ln[20]:.1f} dB")
                ax2.text(col2, 0.65, f"L30: {stats.ln[30]:.1f} dB")
                ax2.text(col2, 0.50, f"L40: {stats.ln[40]:.1f} dB")
                ax2.text(col2, 0.35, f"L60: {stats.ln[60]:.1f} dB")
                ax2.text(col2, 0.20, f"L80: {stats.ln[80]:.1f} dB")
                ax2.text(col3, 0.80, f"Dose ({stats.dose['standard']})", fontweight='bold')
                ax2.text(col3, 0.65, f"Dose %: {stats.dose['dose']:.1f}%")
                ax2.text(col3, 0.50, f"TWA: {stats.dose['twa']:.1f} dB")
                fig.tight_layout()

            case 0: # LEVEL VS TIME
                ax = fig.add_subplot(1, 1, 1)
                t = [r['time'] for r in results]
                l = [r['lp'] for r in results]
                ax.plot(t, l)
                ax.set_title(f"Sound Pressure Level vs Time ({weighting}-weighted, {speed})")
                ax.set_xlabel("Time (s)")
                ax.set_ylabel("Level (dB)")
                ax.grid(True)

            case 2 | 3: # SPECTRAL
                ax = fig.add_subplot(1, 1, 1)
                freqs = results[0]['band_freqs']
                energy_sums = np.zeros(len(freqs))
                for r in results:
                    pressures = (10**(r['bands']/10.0)) * (20e-6**2)
                    energy_sums += pressures
                mean_db = 10 * np.log10((energy_sums / len(results)) / (20e-6**2) + 1e-30)
                
                x = np.arange(len(freqs))
                ax.bar(x, mean_db, color='#2ca02c', alpha=0.8)
                ax.set_xticks(x)
                lbls = []
                for f in freqs:
                    if f >= 1000: lbls.append(f"{f/1000:.0f}k")
                    else: lbls.append(f"{f:.0f}")
                if mode_id == 3: 
                    lbls = [l if i%3==0 else "" for i, l in enumerate(lbls)]
                ax.set_xticklabels(lbls, rotation=90)
                ax.set_title(f"Average Spectrum ({weighting}-weighted)")
                ax.set_ylabel("Level (dB)")
                ax.grid(axis='y')

        self.plot_panel.draw()

    def closeEvent(self, event):
        """
        Graceful shutdown: Prevent zombie threads if user closes window during analysis.
        """
        if self.worker is not None and self.worker.isRunning():
            reply = QMessageBox.question(
                self, 
                'Analysis Running',
                "An analysis is currently running.\nDo you want to stop it and exit?",
                QMessageBox.Yes | QMessageBox.No, 
                QMessageBox.No
            )

            if reply == QMessageBox.Yes:
                self.status_bar.showMessage("Stopping background thread...")
                self.worker.stop()
                self.worker.wait() # Block until thread cleanly exits
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