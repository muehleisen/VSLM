import traceback
from pathlib import Path
from PySide6.QtCore import QObject, Signal, Slot
from .settings_manager import SettingsManager, AppSettings
from .gui.analysis_worker import AnalysisWorker
from .result_exporter import ResultsExporter
from .constants import Weighting, ResponseSpeed, LEQ_INTERVAL_MAP

class VSLMController(QObject):
    sig_file_loaded = Signal(object, object)
    sig_analysis_started = Signal(str)
    sig_analysis_progress = Signal(int)
    sig_analysis_finished = Signal(list)
    sig_analysis_error = Signal(str)
    sig_status_message = Signal(str)
    sig_export_finished = Signal()
    
    # NEW: Signal to update progress bar maximum correctly
    sig_total_blocks = Signal(int)

    def __init__(self):
        super().__init__()
        self.settings_mgr = SettingsManager()
        self.settings = self.settings_mgr.load()
        
        self.filepath: Path | None = None
        self.start_time: float = 0.0
        self.end_time: float | None = None
        
        self.last_results: list = []
        self.worker: AnalysisWorker | None = None
        self.cal_factor = self.settings.calibration_factor
        
        self.total_blocks_estimate = 100

    def load_file(self, filepath: str):
        path = Path(filepath)
        try:
            from soundfile import info
            inf = info(str(path))
            
            self.filepath = path
            self.settings.last_directory = str(path.parent)
            self.start_time = 0.0
            self.end_time = inf.duration
            self.last_results = []
            
            self.sig_file_loaded.emit(path, inf)
            self.sig_status_message.emit("File loaded successfully.")
        except Exception as e:
            self.sig_analysis_error.emit(f"Failed to load file: {e}")

    def update_calibration(self, new_factor: float):
        self.cal_factor = new_factor
        self.settings.calibration_factor = new_factor
        self.sig_status_message.emit(f"Calibration updated: {new_factor:.4f}")

    def set_analysis_range(self, start: float, end: float):
        self.start_time = start
        self.end_time = end
        self.sig_status_message.emit(f"Analysis range set: {start:.2f}s - {end:.2f}s")

    def run_analysis(self, mode_id: int):
        if not self.filepath: return

        do_bands = (mode_id in [2, 3])
        band_res = 'third' if mode_id == 3 else 'octave'
        
        w = self.settings.weighting
        weighting_val = w.value if hasattr(w, 'value') else w
        
        s = self.settings.speed
        speed_val = s.value if hasattr(s, 'value') else s
        
        if self.worker and self.worker.isRunning():
            self.worker.stop()
        
        self.worker = AnalysisWorker(
            filepath=self.filepath,
            cal_factor=self.cal_factor,
            block_size_ms=self.settings.block_size_ms,
            weighting=weighting_val,
            do_bands=do_bands,
            band_res=band_res,
            speed=speed_val,
            band_order=self.settings.band_filter_order,
            ref_pressure=self.settings.ref_pressure
        )

        # Forward the true total blocks from worker to the UI
        self.worker.sig_total_blocks.connect(self.sig_total_blocks.emit)
        
        self.worker.sig_progress.connect(self.sig_analysis_progress.emit)
        self.worker.sig_error.connect(self.sig_analysis_error.emit)
        self.worker.sig_finished.connect(self._on_worker_finished)
        
        self.worker.start()
        self.sig_analysis_started.emit(str(speed_val))

    def stop_analysis(self):
        if self.worker:
            self.worker.stop()
            self.worker.wait()
            self.sig_status_message.emit("Analysis stopped by user.")

    def _on_worker_finished(self, results):
        if self.end_time:
            filtered = [r for r in results if self.start_time <= r['time'] <= self.end_time]
        else:
            filtered = results
            
        self.last_results = filtered
        self.sig_analysis_finished.emit(filtered)
        self.sig_status_message.emit("Analysis Complete.")

    def export_results(self, path: Path, mode_id: int, leq_interval_key):
        if not self.last_results: return

        try:
            w = self.settings.weighting
            weighting = w.value if hasattr(w, 'value') else w
            
            s = self.settings.speed
            speed = s.value if hasattr(s, 'value') else s
            
            dose_std_name = self.settings.current_dose_standard
            dose_params = self.settings.dose_standards.get(dose_std_name)
            
            if mode_id == 1: # LEQ
                 ResultsExporter.export_leq(
                    path, 
                    self.last_results, 
                    self.settings.block_size_ms, 
                    leq_interval_key, 
                    weighting,
                    dose_params,
                    dose_std_name,
                    self.settings.ref_pressure
                )
            elif mode_id == 0: # Lp
                ResultsExporter.export_lp(path, self.last_results, weighting, speed)
            elif mode_id in [2, 3]: # Spectrum
                ResultsExporter.export_spectrum(path, self.last_results, weighting, self.settings.ref_pressure)
                
            self.sig_export_finished.emit()
            self.sig_status_message.emit(f"Exported to {path.name}")
        except Exception as e:
            self.sig_analysis_error.emit(f"Export failed: {str(e)}\n{traceback.format_exc()}")

    def save_settings(self, path: Path):
        self.settings_mgr.save(self.settings, path)
        self.sig_status_message.emit("Settings saved.")

    def load_settings(self, path: Path) -> bool:
        new_settings = self.settings_mgr.load(path)
        if new_settings:
            self.settings = new_settings
            self.cal_factor = self.settings.calibration_factor
            return True
        return False

    def shutdown(self):
        if self.worker and self.worker.isRunning():
            self.worker.stop()
            self.worker.wait()
        self.settings_mgr.save(self.settings)