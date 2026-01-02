import numpy as np
import scipy.signal

# Frequency Range Constants
MIN_SUPPORTED_FS = 20000
MAX_SUPPORTED_FS = 192000

# Standard rates for UI/Testing
STANDARD_RATES = [22050, 44100, 48000, 88200, 96000, 176400, 192000]

def _design_section_zpk(zeros_a, poles_a, fs):
    """
    Helper to create a digital SOS section from analog ZPK using standard bilinear transform.
    """
    z_d, p_d, k_d = scipy.signal.bilinear_zpk(zeros_a, poles_a, 1.0, fs)
    return scipy.signal.zpk2sos(z_d, p_d, k_d)

def _get_ideal_db(freq, w_type):
    """Calculates the exact theoretical ANSI dB value for a target frequency."""
    f = freq
    f2 = f**2
    f1, f2_c, f3, f4 = 20.598997, 107.65265, 737.86223, 12194.217
    
    if w_type == 'A':
        num = (f4**2) * (f2**2)
        den = (f2 + f1**2) * np.sqrt((f2 + f2_c**2) * (f2 + f3**2)) * (f2 + f4**2)
    else: # C
        num = (f4**2) * f2
        den = (f2 + f1**2) * (f2 + f4**2)
        
    gain = num / den
    
    # Normalize to 1kHz ideal
    r2 = 1000.0**2
    if w_type == 'A':
        num_r = (f4**2) * (r2**2)
        den_r = (r2 + f1**2) * np.sqrt((r2 + f2_c**2) * (r2 + f3**2)) * (r2 + f4**2)
    else:
        num_r = (f4**2) * r2
        den_r = (r2 + f1**2) * (r2 + f4**2)
    
    return 20 * np.log10(gain / (num_r / den_r))

def design_variable_weighting_sos(fs, weighting_type='A'):
    """
    Generates weighting filters using a SPLIT-SECTION approach.
    
    - Low-freq poles are transformed using standard bilinear (preserving bass accuracy).
    - High-freq poles are optimized independently to fix 16kHz cramping.
    - Result is an SOS cascade (Series Biquads).
    """
    # ANSI S1.42 Analog Constants
    f1 = 20.598997
    f2 = 107.65265
    f3 = 737.86223
    base_f4 = 12194.217
    
    w1 = 2 * np.pi * f1
    w2 = 2 * np.pi * f2
    w3 = 2 * np.pi * f3
    
    sos_blocks = []
    
    # --- PHASE 1: Low Frequency Sections (Standard Bilinear) ---
    # These determine the shape from 20Hz to ~1kHz. They are perfect as-is.
    
    if weighting_type == 'A':
        # Section 1: Highpass behavior (Poles at w1, Zeros at 0)
        # H1(s) = s^2 / (s+w1)^2
        sos_1 = _design_section_zpk([0, 0], [-w1, -w1], fs)
        sos_blocks.append(sos_1)
        
        # Section 2: Bandpass behavior (Poles at w2, w3, Zeros at 0)
        # H2(s) = s^2 / ((s+w2)(s+w3))
        sos_2 = _design_section_zpk([0, 0], [-w2, -w3], fs)
        sos_blocks.append(sos_2)
        
    elif weighting_type == 'C':
        # Section 1: Highpass behavior (Poles at w1, Zeros at 0)
        # H1(s) = s^2 / (s+w1)^2
        sos_1 = _design_section_zpk([0, 0], [-w1, -w1], fs)
        sos_blocks.append(sos_1)

    # --- PHASE 2: High Frequency Section (Optimized Bilinear) ---
    # This determines the rolloff at >10kHz. 
    # At low sample rates (44.1k), standard bilinear forces -inf dB at 22k.
    # We iteratively shift 'f4' UPWARDS for this section only to lift the curve at 16k.
    
    # Combine existing low-freq blocks to measure current state
    sos_low = np.vstack(sos_blocks)
    
    # Determine f4 target
    # For A/C weighting, H_high(s) = 1 / (s+w4)^2  (Pure Lowpass)
    # Note: Bilinear of 1/(s+w)^2 gives zeros at Nyquist (-1).
    
    best_f4 = base_f4
    
    # Optimization Loop (Only for potentially cramped rates)
    if fs < 50000:
        target_db = _get_ideal_db(16000, weighting_type)
        w_check = 2 * np.pi * 16000 / fs
        
        # Find optimal f4 shift
        best_err = 100.0
        
        # Search range: Nominal 12k up to 24k (pushing past Nyquist is valid math here)
        for test_f4 in np.linspace(base_f4, 24000, 50):
            w4_test = 2 * np.pi * test_f4
            
            # Generate Candidate High Section (Lowpass)
            # Zeros must be empty [] for analog lowpass, bilinear adds zeros at -1
            sos_high_cand = _design_section_zpk([], [-w4_test, -w4_test], fs)
            
            # Combine temporarily
            sos_cand = np.vstack((sos_low, sos_high_cand))
            
            # Normalize Candidate at 1kHz
            w_ref = 2 * np.pi * 1000.0 / fs
            _, h_ref = scipy.signal.sosfreqz(sos_cand, worN=[w_ref])
            gain_corr = 1.0 / (np.abs(h_ref[0]) + 1e-15)
            sos_cand[0, :3] *= gain_corr
            
            # Check 16kHz
            _, h_check = scipy.signal.sosfreqz(sos_cand, worN=[w_check])
            curr_db = 20 * np.log10(np.abs(h_check[0]) + 1e-15)
            
            err = abs(curr_db - target_db)
            if err < best_err:
                best_err = err
                best_f4 = test_f4
            
            if err < 0.05: break

    # Generate Final High Section
    w4_final = 2 * np.pi * best_f4
    sos_high = _design_section_zpk([], [-w4_final, -w4_final], fs)
    
    # Stack Final Filter
    sos_final = np.vstack((sos_low, sos_high))
    
    # --- PHASE 3: Final Normalization ---
    # Ensure 0dB at 1kHz exactly
    w_ref = 2 * np.pi * 1000.0 / fs
    _, h_ref = scipy.signal.sosfreqz(sos_final, worN=[w_ref])
    gain_correction = 1.0 / (np.abs(h_ref[0]) + 1e-15)
    
    # Apply gain to the first section only
    sos_final[0, :3] *= gain_correction
    
    return sos_final

def apply_weighting(data, fs, weighting_type):
    """
    Applies A or C weighting using SOS filtering.
    """
    if weighting_type in ['Z', 'Flat', None]:
        return data

    if not (MIN_SUPPORTED_FS <= fs <= MAX_SUPPORTED_FS):
        raise ValueError(
            f"Sampling rate {fs} Hz is outside the supported range "
            f"({MIN_SUPPORTED_FS} Hz - {MAX_SUPPORTED_FS} Hz)."
        )

    try:
        # Generate SOS
        sos = design_variable_weighting_sos(fs, weighting_type)
        # Apply filter
        return scipy.signal.sosfilt(sos, data)
        
    except Exception as e:
        raise ValueError(f"Could not design {weighting_type}-weighting filter for {fs}Hz: {e}")