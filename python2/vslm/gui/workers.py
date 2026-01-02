# python2/vslm/gui/workers.py
from PySide6.QtCore import QThread, Signal
from pathlib import Path
import traceback

from ..analysis_engine import StreamProcessor

class AnalysisWorker(QThread):
    """
    Background worker thread for running the VSLM analysis engine.
    Emits signals for progress updates and final results.
    """
    # Signals must be class attributes
    sig_progress = Signal(int)              # Emits current block count
    sig_total_blocks = Signal(int)          # Emits total expected blocks
    sig_finished = Signal(list)             # Emits list of result dicts
    sig_error = Signal(str)                 # Emits error message if failed
    
    def __init__(self, 
                 filepath: Path, 
                 cal_factor: float,
                 block_size_ms: float,
                 weighting: str,
                 do_bands: bool,
                 band_res: str,
                 speed: str):
        super().__init__()
        self.filepath = filepath
        self.cal_factor = cal_factor
        self.block_size_ms = block_size_ms
        self.weighting = weighting
        self.do_bands = do_bands
        self.band_res = band_res
        self.speed = speed
        
        self._is_running = True

    def stop(self):
        """Request the worker to stop processing."""
        self._is_running = False

    def run(self):
        try:
            processor = StreamProcessor(self.filepath, self.cal_factor)
            
            # Calculate total blocks for progress bar
            total_blocks = int(processor.duration * 1000 / self.block_size_ms)
            self.sig_total_blocks.emit(total_blocks)
            
            # Initialize Generator
            gen = processor.run_analysis(
                block_size_ms=self.block_size_ms,
                weighting=self.weighting,
                do_band_analysis=self.do_bands,
                band_resolution=self.band_res,
                time_weighting=self.speed
            )
            
            results = []
            
            # Consume generator
            for i, block in enumerate(gen):
                if not self._is_running:
                    break
                    
                results.append(block)
                
                # Emit progress every 10 blocks to reduce signal overhead
                if i % 10 == 0:
                    self.sig_progress.emit(i + 1)
            
            # Send final results if not crashed
            if self._is_running:
                self.sig_progress.emit(total_blocks)
                self.sig_finished.emit(results)
            else:
                # If stopped, still send what we have? Or empty?
                # Usually better to send nothing or partial. sending partial:
                self.sig_finished.emit(results)
                
        except Exception as e:
            # Format the full traceback for debugging
            error_msg = f"{str(e)}\n\n{traceback.format_exc()}"
            self.sig_error.emit(error_msg)