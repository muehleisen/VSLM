# python2/vslm/gui/workers.py
from PySide6.QtCore import QThread, Signal
from pathlib import Path
import traceback
from ..analysis_engine import StreamProcessor

class AnalysisWorker(QThread):
    # ... Signals unchanged ...
    sig_progress = Signal(int)
    sig_total_blocks = Signal(int)
    sig_finished = Signal(list)
    sig_error = Signal(str)

    def __init__(self, 
                 filepath: Path, 
                 cal_factor: float,
                 block_size_ms: float,
                 weighting: str,
                 do_bands: bool,
                 band_res: str,
                 speed: str,
                 # New Params
                 band_order: int,
                 ref_pressure: float):
        super().__init__()
        self.filepath = filepath
        self.cal_factor = cal_factor
        self.block_size_ms = block_size_ms
        self.weighting = weighting
        self.do_bands = do_bands
        self.band_res = band_res
        self.speed = speed
        self.band_order = band_order
        self.ref_pressure = ref_pressure
        self._is_running = True

    def stop(self):
        self._is_running = False

    def run(self):
        try:
            processor = StreamProcessor(self.filepath, self.cal_factor)
            total_blocks = int(processor.duration * 1000 / self.block_size_ms)
            self.sig_total_blocks.emit(total_blocks)
            
            gen = processor.run_analysis(
                block_size_ms=self.block_size_ms,
                weighting=self.weighting,
                do_band_analysis=self.do_bands,
                band_resolution=self.band_res,
                band_order=self.band_order,   # Pass
                time_weighting=self.speed,
                ref_pressure=self.ref_pressure # Pass
            )
            
            results = []
            for i, block in enumerate(gen):
                if not self._is_running: break
                results.append(block)
                if i % 10 == 0: self.sig_progress.emit(i + 1)
            
            if self._is_running:
                self.sig_progress.emit(total_blocks)
                self.sig_finished.emit(results)
            else:
                self.sig_finished.emit(results)
                
        except Exception as e:
            self.sig_error.emit(f"{str(e)}\n\n{traceback.format_exc()}")