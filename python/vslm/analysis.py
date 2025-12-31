# vslm/analysis.py
import numpy as np
import scipy.signal

def calculate_lp(core, weighting='A', speed='Slow'):
    """
    Calculates Sound Pressure Level (Lp) over time.
    """
    data = core.audio_data * core.cal_factor
    fs = core.fs
    
    # 1. Apply Frequency Weighting
    if weighting != 'Z':
        data = core.apply_weighting_filter(data, weighting)
        
    # 2. Apply Time Weighting (Exponential Averaging)
    if speed == 'Slow':
        tau = 1.0
    elif speed == 'Fast':
        tau = 0.125
    else:
        tau = 0.035 # Impulse
        
    alpha = 1.0 - np.exp(-1.0 / (fs * tau))
    b = [alpha]
    a = [1.0, -(1.0 - alpha)]
    
    squared_signal = data**2
    mean_square = scipy.signal.lfilter(b, a, squared_signal)
    
    # Convert to dB
    lp_db = 10 * np.log10(np.maximum(mean_square, 1e-15)) + 94 
    t = np.arange(len(lp_db)) / fs
    
    return t, lp_db

def calculate_leq(core, weighting='A'):
    """
    Calculates Equivalent Continuous Sound Level (Leq).
    """
    data = core.audio_data * core.cal_factor
    if weighting != 'Z':
        data = core.apply_weighting_filter(data, weighting)
        
    ms_val = np.mean(data**2)
    leq = 10 * np.log10(ms_val / (20e-6)**2)
    return leq

def calculate_psd(core):
    """
    Calculates Power Spectral Density (PSD).
    """
    data = core.audio_data * core.cal_factor
    fs = core.fs
    f, Pxx = scipy.signal.welch(data, fs, nperseg=4096, scaling='density')
    Pxx_db = 10 * np.log10(np.maximum(Pxx, 1e-15))
    return f, Pxx_db

def calculate_ansi_bands(core, weighting='Z', resolution='octave'):
    """
    Calculates Octave or 1/3 Octave band levels using FFT synthesis.
    """
    data = core.audio_data * core.cal_factor
    fs = core.fs
    
    if weighting != 'Z':
        data = core.apply_weighting_filter(data, weighting)
        
    N = len(data)
    nfft = 2**int(np.ceil(np.log2(N)))
    Y = np.fft.rfft(data, n=nfft)
    f = np.fft.rfftfreq(nfft, 1/fs)
    mag_sq = np.abs(Y)**2 / (N * fs)
    
    if resolution == 'octave':
        nominal_centers = [31.5, 63, 125, 250, 500, 1000, 2000, 4000, 8000, 16000]
        factor = 2**(1/2)
    else: 
        nominal_centers = [
            25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 
            1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000, 12500, 16000, 20000
        ]
        factor = 2**(1/6)

    levels = []
    freqs_out = []
    
    for fc in nominal_centers:
        fl = fc / factor
        fu = fc * factor
        indices = np.where((f >= fl) & (f <= fu))[0]
        if len(indices) > 0:
            band_power = np.sum(mag_sq[indices])
            band_p2 = band_power * (fs/nfft) * 2 
            band_db = 10 * np.log10(np.maximum(band_p2, 1e-15) / (20e-6)**2)
            levels.append(band_db)
            freqs_out.append(fc)
            
    return freqs_out, levels

def calculate_spectrogram(core, weighting='Z', nfft=4096, slice_len=None, overlap_ratio=0.5, progress_callback=None):
    """
    Computes Spectrogram via manual STFT to support progress reporting.
    
    Args:
        core: VSLMCore instance.
        weighting: 'A', 'C', or 'Z'.
        nfft: FFT size (determines frequency resolution).
        slice_len: Target time step (determines overlap).
        progress_callback: Optional function(int) to report 0-100% progress.
    """
    # 1. Prepare Data
    sig = core.audio_data * core.cal_factor
    fs = core.fs
    
    if weighting != 'Z':
        sig = core.apply_weighting_filter(sig, weighting)
        
    # 2. Window & Step Calculation
    # We lock window size to NFFT to ensure consistent frequency resolution
    nperseg = int(nfft)
    
    if slice_len:
        # Determine overlap based on desired time slice
        samples_per_step = int(slice_len * fs)
        noverlap = nperseg - samples_per_step
        
        # Safety bounds
        if noverlap < 0: noverlap = 0 
        elif noverlap >= nperseg: noverlap = nperseg - 1
    else:
        noverlap = int(nperseg * overlap_ratio)
        
    step = nperseg - noverlap
    
    # 3. Setup STFT Loop
    win = scipy.signal.get_window('hann', nperseg)
    
    # Scale factor for Power Spectral Density (PSD) mode
    # 1 / (Fs * sum(win^2))
    scale_factor = 1.0 / (fs * np.sum(win**2))
    
    n_samples = len(sig)
    
    # Handle short files
    if n_samples < nperseg:
        sig = np.pad(sig, (0, nperseg - n_samples))
        n_samples = len(sig)
        
    num_segments = (n_samples - noverlap) // step
    
    Sxx_list = []
    
    # 4. Processing Loop
    for i in range(num_segments):
        start = i * step
        end = start + nperseg
        
        # Windowing
        segment = sig[start:end] * win
        
        # Real FFT
        spec = np.fft.rfft(segment, n=nfft)
        
        # Magnitude Squared & PSD Scaling
        psd = (np.abs(spec)**2) * scale_factor
        
        # One-sided Spectrum Correction (multiply non-DC/non-Nyquist components by 2)
        if nfft % 2 == 0:
            psd[1:-1] *= 2
        else:
            psd[1:] *= 2
            
        Sxx_list.append(psd)
        
        # Update Progress Bar (every 50 frames to reduce overhead)
        if progress_callback and i % 50 == 0:
            pct = int(100 * (i + 1) / num_segments)
            progress_callback(pct)

    # Final 100% update
    if progress_callback:
        progress_callback(100)

    # 5. Format Output
    if Sxx_list:
        Sxx = np.vstack(Sxx_list).T # Shape: (Freqs, Time)
    else:
        Sxx = np.zeros((nfft//2 + 1, 1))

    # Frequency Vector
    f = np.fft.rfftfreq(nfft, 1/fs)
    
    # Time Vector (Centered on window)
    t = (np.arange(num_segments) * step + nperseg / 2.0) / fs
    
    # Convert to dB
    Sxx_db = 10 * np.log10(np.maximum(Sxx, 1e-15))
    
    return f, t, Sxx_db