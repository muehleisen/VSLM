# python2/vslm/analysis_engine.py
import numpy as np
import soundfile as sf
import os
import warnings

# Import the new stateful filters
from .filters.weighting import WeightingFilter
from .filters.ansi import OctaveFilterBank

class StreamProcessor:
    """
    Orchestrates the streaming analysis of an audio file.
    Reads audio in blocks, applies filters, and yields time-history results.
    """
    def __init__(self, filepath, cal_factor=1.0):
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"File not found: {filepath}")
            
        self.filepath = filepath
        self.cal_factor = cal_factor
        
        # Read file info without loading data into RAM
        try:
            info = sf.info(filepath)
            self.fs = info.samplerate
            self.duration = info.duration
            self.channels = info.channels
        except Exception as e:
            raise ValueError(f"Could not read file info: {e}")
        
    def run_analysis(self, 
                     block_size_ms=100, 
                     weighting='A', 
                     do_band_analysis=False, 
                     band_resolution='octave',
                     band_order=24):
        """
        Generator that yields analysis results for each time block.
        
        Args:
            block_size_ms (float): Size of processing chunks in milliseconds.
            weighting (str): Frequency weighting ('A', 'C', 'Z').
            do_band_analysis (bool): Whether to compute octave bands.
            band_resolution (str): 'octave' or 'third'.
            band_order (int): Filter order for bands (default 24 for Class 1).
            
        Yields:
            dict: A dictionary containing 'time', 'leq', and optional 'bands'.
        """
        # 1. Initialize Stateful Filters
        # These classes now handle their own 'zi' state internally.
        weighting_filter = WeightingFilter(self.fs, weighting)
        
        band_bank = None
        if do_band_analysis:
            band_bank = OctaveFilterBank(self.fs, resolution=band_resolution, order=band_order)

        # 2. Calculate Block Size
        # Ensure block size corresponds to an integer number of samples
        block_samples = int(self.fs * (block_size_ms / 1000.0))
        if block_samples == 0:
            raise ValueError("Block size is too small for this sampling rate.")
        
        # 3. Stream Processing Loop
        current_time = 0.0
        
        # sf.blocks reads the file chunk-by-chunk efficiently
        # fill_value=0.0 ensures the last chunk is padded if it's too short
        block_gen = sf.blocks(
            self.filepath, 
            blocksize=block_samples, 
            always_2d=False, 
            fill_value=0.0
        )
        
        for chunk in block_gen:
            # Handle multi-channel: Downmix to mono (average)
            # Future versions could allow channel selection
            if chunk.ndim > 1:
                chunk = np.mean(chunk, axis=1)

            # A. Apply Calibration
            # P_pa = P_normalized * CalibrationFactor
            calibrated_chunk = chunk * self.cal_factor
            
            # B. Apply Frequency Weighting (Stateful)
            weighted_chunk = weighting_filter.process_chunk(calibrated_chunk)
            
            # C. Compute Broadband Leq (for this block)
            # Add tiny epsilon to avoid log10(0)
            ms_broadband = np.mean(weighted_chunk**2)
            leq_block = 10 * np.log10(ms_broadband / (20e-6)**2 + 1e-30)
            
            result = {
                'time': current_time,
                'leq': leq_block
            }
            
            # D. Compute Band Levels (Optional)
            if band_bank:
                # Returns matrix [n_samples, n_bands]
                filtered_bands = band_bank.process_chunk(calibrated_chunk)
                
                # Compute RMS for each band in this block
                ms_bands = np.mean(filtered_bands**2, axis=0)
                leq_bands = 10 * np.log10(ms_bands / (20e-6)**2 + 1e-30)
                
                result['bands'] = leq_bands
                result['band_freqs'] = band_bank.frequencies

            # Yield result to the GUI or Data Accumulator
            yield result
            
            current_time += (block_size_ms / 1000.0)