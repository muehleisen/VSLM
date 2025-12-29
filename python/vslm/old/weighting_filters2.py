import numpy as np
import scipy.signal

# Frequency Range Constants
MIN_SUPPORTED_FS = 20000
MAX_SUPPORTED_FS = 192000

# Standard rates for UI/Testing
STANDARD_RATES = [22050, 44100, 48000, 88200, 96000, 176400, 192000]

def _get_ideal_db(freq, w_type):
    """Calculates the exact theoretical ANSI dB value for a target frequency."""
    f = freq
    f2 = f**2
    # ANSI Constants
    f1, f2_c, f3, f4 = 20.598997, 107.65265, 737.86223, 12194.217
    
    if w_type == 'A':
        num = (f4**2) * (f2**2)
        den = (f2 + f1**2) * np.sqrt((f2 + f2_c**2) * (f2 + f3**2)) * (f2 + f4**2)
    else: # C
        num = (f4**2) * f2
        den = (f2 + f1**2) * (f2 + f4**2)
        
    gain = num / den
    
    # Normalize to 1kHz
    ref = 1000.0
    r2 = ref**2
    if w_type == 'A':
        num_r = (f4**2) * (r2**2)
        den_r = (r2 + f1**2) * np.sqrt((r2 + f2_c**2) * (r2 + f3**2)) * (r2 + f4**2)
    else:
        num_r = (f4**2) * r2
        den_r = (r2 + f1**2) * (r2 + f4**2)
    
    return 20 * np.log10(gain / (num_r / den_r))

def _get_coeffs_raw(fs, weighting_type, f4_override=None):
    """
    Internal helper to generate coefficients with a specific f4 pole frequency.
    Used by the optimizer to 'nudge' the filter into compliance.
    """
    f1 = 20.598997
    f2 = 107.65265
    f3 = 737.86223
    f4 = 12194.217 if f4_override is None else f4_override
    
    # Convert to Angular Frequency
    w1 = 2 * np.pi * f1
    w2 = 2 * np.pi * f2
    w3 = 2 * np.pi * f3
    w4 = 2 * np.pi * f4
    
    if weighting_type == 'A':
        zeros_a = np.array([0, 0, 0, 0])
        poles_a = -np.array([w1, w1, w2, w3, w4, w4])
    elif weighting_type == 'C':
        zeros_a = np.array([0, 0])
        poles_a = -np.array([w1, w1, w4, w4])
    else:
        raise ValueError(f"Unknown weighting: {weighting_type}")

    # Standard Bilinear Transform
    z_d, p_d, k_d = scipy.signal.bilinear_zpk(zeros_a, poles_a, 1.0, fs)
    b, a = scipy.signal.zpk2tf(z_d, p_d, k_d)
    
    # Gain Normalization at 1000 Hz
    w_ref = 2 * np.pi * 1000.0 / fs
    _, h_ref = scipy.signal.freqz(b, a, worN=[w_ref])
    current_gain = np.abs(h_ref[0])
    
    if current_gain > 0:
        b = b / current_gain
        
    return b, a

def design_variable_weighting(fs, weighting_type='A'):
    """
    Generates ANSI S1.42 compliant coefficients.
    
    Includes an OPTIMIZATION LOOP for low sample rates (< 50kHz).
    It iteratively shifts the high-frequency pole (f4) to force the 
    16kHz response to meet the exact ANSI target.
    """
    # 1. Initial Design
    b, a = _get_coeffs_raw(fs, weighting_type)
    
    # 2. Optimization Check
    # We run this for both A and C if the sample rate is low enough to cause cramping.
    if fs < 50000:
        # Determine Check Frequency: 16kHz is usually the hardest point.
        # But for 22050Hz (Nyquist 11k), we can't check 16k.
        nyquist = fs / 2.0
        
        # If Nyquist is too low to even support the check freq, skip optimization
        if nyquist > 16000:
            f_check = 16000
        elif nyquist > 12500:
            f_check = 12000
        else:
            return b, a # Rate too low to optimize for 16k compliance anyway

        # Calculate exact ideal target at f_check
        target_db = _get_ideal_db(f_check, weighting_type)
        
        # Measure current performance
        w_check = 2 * np.pi * f_check / fs
        _, h_check = scipy.signal.freqz(b, a, worN=[w_check])
        curr_db = 20 * np.log10(np.abs(h_check[0]) + 1e-15)
        
        # If error > 0.5dB, optimize
        if abs(curr_db - target_db) > 0.5:
            # print(f"Optimizing {weighting_type}-weight @ {fs}Hz (Target {target_db:.2f}dB)...")
            
            best_b, best_a = b, a
            best_err = abs(curr_db - target_db)
            
            base_f4 = 12194.217
            
            # Expanded Search Range: 0 to +8000 Hz shift.
            # 44.1k often needs f4 to be pushed towards ~18-19k to keep 16k up.
            steps = 40
            for shift in np.linspace(0, 8000, steps):
                test_f4 = base_f4 + shift
                b_test, a_test = _get_coeffs_raw(fs, weighting_type, f4_override=test_f4)
                
                _, h_test = scipy.signal.freqz(b_test, a_test, worN=[w_check])
                db_test = 20 * np.log10(np.abs(h_test[0]) + 1e-15)
                
                err = abs(db_test - target_db)
                
                # We prefer slight 'overshoot' (lifting up) vs droop, 
                # but absolute error minimization is usually safe.
                if err < best_err:
                    best_err = err
                    best_b, best_a = b_test, a_test
                
                # Strict exit condition
                if err < 0.05:
                    break
            
            b, a = best_b, best_a

    return b, a

def apply_weighting(data, fs, weighting_type):
    """Applies A or C weighting to an audio signal."""
    if weighting_type in ['Z', 'Flat', None]:
        return data

    if not (MIN_SUPPORTED_FS <= fs <= MAX_SUPPORTED_FS):
        raise ValueError(f"Sampling rate {fs} Hz is outside supported range.")

    try:
        b, a = design_variable_weighting(fs, weighting_type)
    except Exception as e:
        raise ValueError(f"Filter design failed for {fs}Hz: {e}")
            
    return scipy.signal.lfilter(b, a, data)