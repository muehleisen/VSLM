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
        """
        # 1. Initialize Stateful Filters
        weighting_filter = WeightingFilter(self.fs, weighting)
        
        band_bank = None
        if do_band_analysis:
            band_bank = OctaveFilterBank(self.fs, resolution=band_resolution, order=band_order)

        # 2. Calculate Block Size
        block_samples = int(self.fs * (block_size_ms / 1000.0))
        if block_samples == 0:
            raise ValueError("Block size is too small for this sampling rate.")
        
        # 3. Open File and Seed Filters (Minimize Startup Glitch)
        # We manually manage the file handle to allow seek operations
        with sf.SoundFile(self.filepath) as f:
            
            # --- START SEEDING LOGIC ---
            # Read the first chunk to "warm up" the filters
            seed_data = f.read(block_samples, always_2d=False, fill_value=0.0)
            
            # Handle multi-channel mixing for seeding
            if seed_data.ndim > 1:
                seed_data = np.mean(seed_data, axis=1)
                
            # Apply calibration to seed data
            seed_data = seed_data * self.cal_factor
            
            # Run the Forward-Backward initialization
            weighting_filter.initialize_state(seed_data)
            if band_bank:
                band_bank.initialize_state(seed_data)
                
            # Reset file pointer to beginning for actual analysis
            f.seek(0)
            # --- END SEEDING LOGIC ---

            # 4. Stream Processing Loop
            current_time = 0.0
            
            # CORRECTED: Use the method f.blocks() on the open object
            # instead of the module function sf.blocks(f)
            block_gen = f.blocks(
                blocksize=block_samples, 
                always_2d=False, 
                fill_value=0.0
            )
            
            for chunk in block_gen:
                if chunk.ndim > 1:
                    chunk = np.mean(chunk, axis=1)

                # A. Apply Calibration
                calibrated_chunk = chunk * self.cal_factor
                
                # B. Apply Frequency Weighting (Stateful)
                weighted_chunk = weighting_filter.process_chunk(calibrated_chunk)
                
                # C. Compute Broadband Leq
                ms_broadband = np.mean(weighted_chunk**2)
                leq_block = 10 * np.log10(ms_broadband / (20e-6)**2 + 1e-30)
                
                result = {
                    'time': current_time,
                    'leq': leq_block
                }
                
                # D. Compute Band Levels (Optional)
                if band_bank:
                    filtered_bands = band_bank.process_chunk(calibrated_chunk)
                    ms_bands = np.mean(filtered_bands**2, axis=0)
                    leq_bands = 10 * np.log10(ms_bands / (20e-6)**2 + 1e-30)
                    
                    result['bands'] = leq_bands
                    result['band_freqs'] = band_bank.frequencies

                yield result
                
                current_time += (block_size_ms / 1000.0)