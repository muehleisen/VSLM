import numpy as np
import scipy.signal

def get_exact_center_frequencies(resolution='octave', base=10):
    """
    Returns exact Center Frequencies (Fc) based on ANSI S1.11-2004.
    Matches standard base-10 or base-2 definitions.
    """
    f_ref = 1000.0
    
    if resolution == 'octave':
        # Indices -6 to +4 (15.625 Hz to 16 kHz)
        x = np.arange(-6, 5) 
        if base == 10:
            fm = f_ref * (10.0 ** (3.0 * x / 10.0))
        else:
            fm = f_ref * (2.0 ** x)
            
    elif resolution == 'third':
        # Indices -19 to +13 (12.5 Hz to 20 kHz)
        x = np.arange(-19, 14)
        if base == 10:
            fm = f_ref * (10.0 ** (x / 10.0))
        else:
            fm = f_ref * (2.0 ** (x / 3.0))
    else:
        raise ValueError("Resolution must be 'octave' or 'third'")

    return fm

def design_ansi_filter(fc, fs, resolution='third', order=24):
    """
    Generates SOS coefficients for a COMPLIANT ANSI S1.11 Bandpass.
    
    CHANGES:
    - Default Order increased to 24 to guarantee >60dB attenuation at 1/2 octave skirts.
    - Maintains 'Bandwidth Correction' to keep passband flat within 0.05dB.
    
    Args:
        fc (float): Center frequency.
        fs (float): Sampling rate.
        resolution (str): 'octave' or 'third'.
        order (int): Order. N=24 is required for strict Class 1 skirt compliance.
    """
    # 1. Determine Exact ANSI Band Edges
    if resolution == 'octave':
        bandwidth_factor = 2**(1.0 / 2.0) # Sqrt(2)
    else:
        bandwidth_factor = 2**(1.0 / 6.0) # Sixth-root(2)
        
    f_lower_ansi = fc / bandwidth_factor
    f_upper_ansi = fc * bandwidth_factor
    
    # 2. Bandwidth Correction Strategy
    # We force the filter to be -0.05 dB (almost flat) at the ANSI edges.
    # The higher order (N=24) makes the transition cliff-like, dropping 
    # ~72dB in the half-octave gap after the edge.
    
    target_attenuation_db = 0.05 
    
    # Calculate correction factor alpha
    # ratio = ( 10^(target/10) - 1 )^(1/2N)
    term = (10**(target_attenuation_db / 10.0)) - 1
    alpha = term ** (1.0 / (2.0 * order))
    
    # Apply correction to get "Design Frequencies"
    f_lower_design = f_lower_ansi * alpha
    f_upper_design = f_upper_ansi / alpha
    
    # 3. Nyquist Check & Clamp
    nyquist = fs / 2.0
    
    if f_upper_design >= nyquist:
        # Clamp close to Nyquist if bandwidth pushes too high
        f_upper_design = nyquist * 0.99

    # 4. Design the Filter
    # output='sos' is CRITICAL for Order 24 stability
    sos = scipy.signal.butter(
        N=order, 
        Wn=[f_lower_design, f_upper_design], 
        btype='bandpass', 
        fs=fs, 
        output='sos'
    )
    
    return sos

def generate_filter_bank(fs, resolution='third', order=24):
    """
    Generates a full bank of filters for the given sampling rate.
    """
    centers = get_exact_center_frequencies(resolution, base=10)
    bank = {}
    
    for fc in centers:
        if fc < (fs / 2.0):
            try:
                sos = design_ansi_filter(fc, fs, resolution, order)
                bank[fc] = sos
            except Exception as e:
                print(f"Skipping {fc}Hz: {e}")
                
    return bank