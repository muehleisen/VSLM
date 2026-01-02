import numpy as np
import scipy.signal

# Standard sampling rates from the original VSLM (for which we have pre-calibrated coefficients)
SUPPORTED_FS = [22050, 44100, 48000, 96000]

# --- Golden Coefficients (extracted from vslm.m) ---
A_WEIGHTING_COEFFS = {
    22050: (
        np.array([0.411130459349308, 0.489163899078409, -1.81511690197873, 0.505982131319440, 0.408895901518515]),
        np.array([1.00000000000000, -0.820486371642917, -0.954967112698490, 0.819769137019238, -0.0308133150879904])
    ),
    44100: (
        np.array([0.218436284008832, -0.0271197371604133, -1.23075568024016, 1.66912752925371, -0.629688397944900]),
        np.array([1.00000000000000, -3.16548213635874, 3.58370028028808, -1.66952594339746, 0.251312323839929])
    ),
    48000: (
        np.array([0.180459241142510, 0.0860957271342434, -1.34233754417886, 1.70455253682352, -0.628769962323045]),
        np.array([1.00000000000000, -3.21098152161547, 3.70934390304575, -1.78459093140491, 0.286231936413823])
    ),
    96000: (
        np.array([-0.0930424680674304, 0.885061570433706, -2.09736158642518, 1.91170857911392, -0.606366095212699]),
        np.array([1.00000000000000, -3.52292032961222, 4.59276248551295, -2.61658164419795, 0.546739755260711])
    )
}

C_WEIGHTING_COEFFS = {
    22050: (
        np.array([0.403673160421101, 0.867939940009175, -0.387869031048492, -0.871337058723055, -0.0191991534215092]),
        np.array([1.00000000000000, 0.487110438920731, -1.19411745882897, -0.473099276797724, 0.208201319236022])
    ),
    44100: (
        np.array([0.205377794209128, 0.731122249663226, 0.318817527196078, -0.733113968701055, -0.526187064231238]),
        np.array([1.00000000000000, 0.720524473037467, -1.27460941890720, -0.711991070511022, 0.283142933498142])
    ),
    48000: (
        np.array([0.210715717951749, 0.718855215372106, 0.245507756808660, -0.716973939314165, -0.462053117024141]),
        np.array([1.00000000000000, 0.678016119469390, -1.28615510047168, -0.687899298538665, 0.312313266344913])
    ),
    96000: (
        np.array([0.0484465929621835, 0.151813149050837, -0.459174213032651, 0.269122623377235, -0.0102081523621007]),
        np.array([1.00000000000000, -2.93939917852114, 3.10634442611301, -1.39372064867588, 0.226775918992396])
    )
}

def design_variable_weighting(fs, weighting_type='A'):
    """
    Generates ANSI S1.42 compliant weighting coefficients (b, a) 
    for ANY sampling rate using the Bilinear Transform.
    """
    # ANSI S1.42 Analog Pole/Zero Definition (in Hz)
    f1 = 20.598997
    f2 = 107.65265
    f3 = 737.86223
    f4 = 12194.217
    
    # Convert to Angular Frequency (rad/s)
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
        raise ValueError(f"Unknown weighting type: {weighting_type}")

    # --- Bilinear Transform ---
    # FIX: bilinear_zpk returns 3 values: Zeros, Poles, Gain
    z_d, p_d, k_d = scipy.signal.bilinear_zpk(zeros_a, poles_a, 1.0, fs)
    
    # Convert ZPK to transfer function polynomials (b, a)
    b, a = scipy.signal.zpk2tf(z_d, p_d, k_d)
    
    # --- Normalization ---
    # The standard requires 0 dB (Gain=1.0) at 1000 Hz.
    w_ref = 2 * np.pi * 1000.0 / fs
    _, h_ref = scipy.signal.freqz(b, a, worN=[w_ref])
    current_gain = np.abs(h_ref[0])
    
    if current_gain > 0:
        b = b / current_gain
    
    return b, a

def apply_weighting(data, fs, weighting_type):
    """
    Applies A or C weighting to an audio signal.
    """
    if weighting_type in ['Z', 'Flat', None]:
        return data
        
    coeffs = None
    
    # 1. Try Hardcoded (Golden) Tables First
    if weighting_type == 'A' and fs in A_WEIGHTING_COEFFS:
        coeffs = A_WEIGHTING_COEFFS[fs]
    elif weighting_type == 'C' and fs in C_WEIGHTING_COEFFS:
        coeffs = C_WEIGHTING_COEFFS[fs]
        
    # 2. If not found, generate dynamically
    if coeffs is None:
        try:
            coeffs = design_variable_weighting(fs, weighting_type)
        except Exception as e:
            raise ValueError(f"Could not design {weighting_type}-weighting filter for {fs}Hz: {e}")
            
    b, a = coeffs
    return scipy.signal.lfilter(b, a, data)