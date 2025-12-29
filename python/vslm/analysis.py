# vslm/analysis.py
import numpy as np
import scipy.signal
from .core import VSLMCore
from .ansi_filters import design_ansi_filter, get_exact_center_frequencies

def calculate_ansi_bands(vslm_core: VSLMCore, weighting='Z', resolution='octave', progress_callback=None):
    """
    Calculates Band Levels using time-domain ANSI S1.11 filters.
    Supports progress reporting via callback.
    """
    fs = vslm_core.fs
    
    # 1. Apply Frequency Weighting
    raw_data = vslm_core.audio_data * vslm_core.cal_factor
    data = vslm_core.apply_weighting_filter(raw_data, weighting)
    
    # 2. Get Center Frequencies
    center_freqs = get_exact_center_frequencies(resolution=resolution, base=10)

    band_levels = []
    valid_centers = []
    
    # Filter out valid bands first to know total count for progress bar
    valid_freqs = [fc for fc in center_freqs if fc < fs / 2.0]
    total_bands = len(valid_freqs)

    # 3. Filter and Compute Levels
    for i, fc in enumerate(valid_freqs):
        # --- Update Progress ---
        if progress_callback:
            # Calculate percentage (0-100)
            percent = int((i / total_bands) * 100)
            progress_callback(percent)
        # -----------------------

        try:
            sos = design_ansi_filter(fc, fs, resolution=resolution, order=24)
        except Exception as e:
            print(f"Skipping band {fc:.1f} Hz: {e}")
            continue
            
        # Apply filter
        filtered_data = scipy.signal.sosfilt(sos, data)
        
        # Compute Level
        mean_sq = np.mean(filtered_data**2)
        ref_pressure_sq = (20e-6)**2
        
        if mean_sq <= 0:
            lvl = 0
        else:
            lvl = 10 * np.log10(mean_sq / ref_pressure_sq)
            
        band_levels.append(lvl)
        valid_centers.append(fc)
        
    # Ensure progress hits 100% at end
    if progress_callback:
        progress_callback(100)
        
    return np.array(valid_centers), np.array(band_levels)

# calculate_psd remains unchanged (it is usually fast enough not to need a progress bar)
def calculate_psd(vslm_core, nfft=4096, window='hann', overlap=0.5, progress_callback=None):
    """Calculates Power Spectral Density."""
    # We add progress_callback argument just to be compatible with the Worker calling convention,
    # even if we don't use it.
    fs = vslm_core.fs
    data = vslm_core.audio_data * vslm_core.cal_factor
    
    f, pxx = scipy.signal.welch(data, fs=fs, window=window, 
                                nperseg=nfft, 
                                noverlap=int(nfft*overlap),
                                scaling='density')
    
    ref_pressure_sq = (20e-6)**2
    lpxx = 10 * np.log10(pxx / ref_pressure_sq)
    
    if progress_callback: progress_callback(100)
    return f, lpxx