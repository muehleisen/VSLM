from PySide6.QtCore import QThread, Signal
from contextlib import closing
import traceback

# VSLM Imports
from ..analysis_engine import StreamProcessor

class AnalysisWorker(QThread):
    """
    Background worker thread to run the audio analysis without freezing the GUI.
    """
    # Signals to communicate with the main thread
    sig_progress = Signal(int)       # Emits current block number
    sig_total_blocks = Signal(int)   # Emits total number of blocks
    sig_finished = Signal(list)      # Emits the full list of results on completion
    sig_error = Signal(str)          # Emits error message if something fails

    def __init__(self, filepath, cal_factor, block_size_ms, weighting, 
                 do_bands, band_res, speed, band_order, ref_pressure):
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

    def run(self):
        try:
            processor = StreamProcessor(self.filepath, self.cal_factor)
            
            # Estimate total blocks for the progress bar
            # (Total Duration / Block Duration)
            total_blocks = int((processor.duration * 1000) / self.block_size_ms)
            self.sig_total_blocks.emit(total_blocks)
            
            results = []
            
            # Create the generator
            gen = processor.run_analysis(
                block_size_ms=self.block_size_ms,
                weighting=self.weighting,
                do_band_analysis=self.do_bands,
                band_resolution=self.band_res,
                time_weighting=self.speed,
                band_order=self.band_order,
                ref_pressure=self.ref_pressure
            )
            
            # --- SAFETY FIX: Use context manager to ensure generator cleanup ---
            # If the loop breaks early (Stop button), 'closing' ensures 
            # the underlying file handle in StreamProcessor is closed immediately.
            with closing(gen):
                for i, block in enumerate(gen):
                    if not self._is_running:
                        break
                    
                    results.append(block)
                    
                    # Update progress every 10 blocks to reduce overhead
                    if i % 10 == 0:
                        self.sig_progress.emit(i + 1)
            
            if self._is_running:
                # Send 100% progress and final data
                self.sig_progress.emit(total_blocks)
                self.sig_finished.emit(results)
                
        except Exception as e:
            # Capture full traceback for debugging
            error_msg = f"{str(e)}\n\nTraceback:\n{traceback.format_exc()}"
            self.sig_error.emit(error_msg)

    def stop(self):
        """Signals the thread to stop processing safely."""
        self._is_running = False