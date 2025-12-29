import numpy as np
import scipy.signal
from .core import VSLMCore
# NEW: Import the compliant filter logic
from .ansi_filters import design_ansi_filter, get_exact_center_frequencies

def calculate_ansi_bands(vslm_core: VSLMCore, weighting='Z', resolution='octave'):
    """
    Calculates Band Levels using time-domain ANSI S1.11 filters.
    
    UPDATED: Now uses 'ansi_filters.py' to ensure Class 1 compliance 
    (High-order SOS + Bandwidth Correction).
    """
    fs = vslm_core.fs
    
    # 1. Apply Frequency Weighting to the raw signal first (A, C, or Z)
    #    (e.g., A-weighted Octave Bands)
    raw_data = vslm_core.audio_data * vslm_core.cal_factor
    data = vslm_core.apply_weighting_filter(raw_data, weighting)
    
    # 2. Get Exact ANSI Center Frequencies
    #    This replaces the manual "1000 * 2^n" logic with the standard-compliant list
    center_freqs = get_exact_center_frequencies(resolution=resolution, base=10)

    band_levels = []
    valid_centers = []

    # 3. Filter and Compute Levels
    for fc in center_freqs:
        # Stop if the band center is above Nyquist (Fs/2)
        if fc >= fs / 2.0:
            break
            
        # --- NEW LOGIC STARTS HERE ---
        # Generate the High-Order, Corrected Filter
        # We use order=24 to guarantee the 60dB skirts you verified in the plots.
        try:
            sos = design_ansi_filter(fc, fs, resolution=resolution, order=24)
        except Exception as e:
            # If a filter cannot be designed (e.g. too close to Nyquist), skip this band
            print(f"Skipping band {fc:.1f} Hz: {e}")
            continue
            
        # Apply the filter (Must use sosfilt for stability!)
        filtered_data = scipy.signal.sosfilt(sos, data)
        # --- NEW LOGIC ENDS HERE ---
        
        # Compute Leq for this band
        # Mean Square = Average Energy
        mean_sq = np.mean(filtered_data**2)
        ref_pressure_sq = (20e-6)**2
        
        if mean_sq <= 0:
            lvl = 0
        else:
            lvl = 10 * np.log10(mean_sq / ref_pressure_sq)
            
        band_levels.append(lvl)
        valid_centers.append(fc)
        
    return np.array(valid_centers), np.array(band_levels)

def calculate_psd(vslm_core, nfft=4096, window='hann', overlap=0.5):
    """
    Calculates Power Spectral Density.
    (Unchanged from previous version)
    """
    fs = vslm_core.fs
    data = vslm_core.audio_data * vslm_core.cal_factor
    
    f, pxx = scipy.signal.welch(data, fs=fs, window=window, 
                                nperseg=nfft, 
                                noverlap=int(nfft*overlap),
                                scaling='density')
    
    ref_pressure_sq = (20e-6)**2
    lpxx = 10 * np.log10(pxx / ref_pressure_sq)
    
    return f, lpxx