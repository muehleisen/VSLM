# python2/vslm/analysis_engine.py
import numpy as np
import soundfile as sf
import os
import warnings

from .filters.weighting import WeightingFilter
from .filters.ansi import OctaveFilterBank

class TimeWeightingDetector:
    """
    Applies IEC 61672-1 Time Weighting (Fast, Slow, Impulse).
    Processes every sample (no downsampling) for maximum accuracy.
    """
    def __init__(self, fs, mode='Fast'):
        self.fs = fs
        self.mode = mode
        self.state = 0.0 # Stores the previous smoothed power value
        
        # Calculate Time Constants
        # Alpha = 1 - exp(-1 / (fs * tau))
        if mode == 'Fast':
            tau = 0.125
            self.alpha_rise = 1.0 - np.exp(-1.0 / (fs * tau))
            self.alpha_fall = self.alpha_rise
        elif mode == 'Slow':
            tau = 1.0
            self.alpha_rise = 1.0 - np.exp(-1.0 / (fs * tau))
            self.alpha_fall = self.alpha_rise
        elif mode == 'Impulse':
            tau_rise = 0.035
            tau_fall = 1.5
            self.alpha_rise = 1.0 - np.exp(-1.0 / (fs * tau_rise))
            self.alpha_fall = 1.0 - np.exp(-1.0 / (fs * tau_fall))
        else:
            self.alpha_rise = 0.1
            self.alpha_fall = 0.1
        

    def process(self, chunk):
        """
        Processes a chunk of audio sample-by-sample and returns the peak Lp level (dB).
        """
        # 1. Calculate Instantaneous Power
        p2 = chunk**2
        
        # 2. Sample-by-Sample Recursive Smoothing
        current_val = self.state
        max_val_in_block = 0.0
        
        # Cache alphas to locals for loop speed optimization
        a_rise = self.alpha_rise
        a_fall = self.alpha_fall
        
        # This loop runs at full Fs (e.g. 48,000 times per second of audio)
        for s in p2:
            if s > current_val:
                # Rising Edge
                current_val = (1 - a_rise) * current_val + a_rise * s
            else:
                # Falling Edge
                current_val = (1 - a_fall) * current_val + a_fall * s
            
            if current_val > max_val_in_block:
                max_val_in_block = current_val
        
        self.state = current_val
        
        # Return dB (prevent log0)
        return 10 * np.log10(max_val_in_block / (20e-6**2) + 1e-30)

class StreamProcessor:
    """
    Orchestrates the streaming analysis of an audio file.
    """
    def __init__(self, filepath, cal_factor=1.0):
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"File not found: {filepath}")
            
        self.filepath = filepath
        self.cal_factor = cal_factor
        
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
                     band_order=24,
                     time_weighting='Fast'): 
        """
        Generator that yields analysis results for each time block.
        """
        # 1. Initialize Filters
        weighting_filter = WeightingFilter(self.fs, weighting)
        lp_detector = TimeWeightingDetector(self.fs, time_weighting)
        
        band_bank = None
        if do_band_analysis:
            band_bank = OctaveFilterBank(self.fs, resolution=band_resolution, order=band_order)

        block_samples = int(self.fs * (block_size_ms / 1000.0))
        if block_samples == 0:
            raise ValueError("Block size is too small.")
        
        with sf.SoundFile(self.filepath) as f:
            # --- Seeding Logic ---
            seed_data = f.read(block_samples, always_2d=False, fill_value=0.0)
            if seed_data.ndim > 1: seed_data = np.mean(seed_data, axis=1)
            seed_data = seed_data * self.cal_factor
            
            weighting_filter.initialize_state(seed_data)
            if band_bank: band_bank.initialize_state(seed_data)
            
            # Seed the time weighting detector (run forward only)
            lp_detector.process(seed_data) 
            
            f.seek(0)
            # ---------------------

            current_time = 0.0
            
            # Use instance method f.blocks()
            block_gen = f.blocks(blocksize=block_samples, always_2d=False, fill_value=0.0)
            
            for chunk in block_gen:
                if chunk.ndim > 1: chunk = np.mean(chunk, axis=1)

                calibrated_chunk = chunk * self.cal_factor
                weighted_chunk = weighting_filter.process_chunk(calibrated_chunk)
                
                # 1. Compute LEQ (Short-term Integration)
                ms_broadband = np.mean(weighted_chunk**2)
                leq_block = 10 * np.log10(ms_broadband / (20e-6)**2 + 1e-30)
                
                # 2. Compute Lp (Time Weighted) - Sample accurate
                lp_block = lp_detector.process(weighted_chunk)
                
                result = {
                    'time': current_time,
                    'leq': leq_block,
                    'lp': lp_block
                }
                
                if band_bank:
                    filtered_bands = band_bank.process_chunk(calibrated_chunk)
                    ms_bands = np.mean(filtered_bands**2, axis=0)
                    leq_bands = 10 * np.log10(ms_bands / (20e-6)**2 + 1e-30)
                    
                    result['bands'] = leq_bands
                    result['band_freqs'] = band_bank.frequencies

                yield result
                current_time += (block_size_ms / 1000.0)