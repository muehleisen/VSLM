import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal

# --- Path Hack: Allows running this script directly without installing the package ---
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
# -----------------------------------------------------------------------------------

from vslm.ansi_filters import design_ansi_filter

def plot_mask_check(fc, fs, resolution='third', order=24):
    """
    Plots the filter response against an approximate ANSI S1.11 Class 1 mask.
    Supports both 'third' and 'octave' resolutions.
    """
    sos = design_ansi_filter(fc, fs, resolution, order)
    
    # Generate Frequency Response
    w, h = scipy.signal.sosfreqz(sos, worN=16384, fs=fs) # High resolution FFT
    freqs = w
    
    # Avoid log(0) error
    mag_db = 20 * np.log10(np.abs(h) + 1e-15) 
    
    # Normalize frequency axis relative to Fc (f / fc)
    f_norm = freqs / fc
    
    plt.figure(figsize=(12, 7))
    plt.semilogx(f_norm, mag_db, 'b-', label=f'{resolution.capitalize()} Filter (N={order})', linewidth=2)
    
    # --- ANSI S1.11 Class 1 Mask (Approximate Visual Guide) ---
    if resolution == 'third':
        # 1/3 Octave Band Edges: +/- 1/6th octave (2^(1/6) approx 1.122)
        lower_edge = 2**(-1.0/6.0) # ~0.8909
        upper_edge = 2**(1.0/6.0)  # ~1.1225
        
        # Passband Tolerance
        plt.fill_between([lower_edge, upper_edge], -0.3, 0.3, color='green', alpha=0.2, label='Passband (+/- 0.3dB)')
        
        # 1/2 Octave Skirts (attenuation required > 17.5dB)
        # Ratios relative to Fc: ~0.707 and ~1.414
        plt.axvline(0.707, color='orange', linestyle=':', alpha=0.5)
        plt.axvline(1.414, color='orange', linestyle=':', alpha=0.5)
        
        # 1 Octave Skirts (attenuation required > 60dB)
        # Ratios relative to Fc: 0.5 and 2.0
        plt.plot([0.5, 0.5], [-5, -80], 'r--', linewidth=1, label='Stopband Limit')
        plt.plot([2.0, 2.0], [-5, -80], 'r--', linewidth=1)
        
    elif resolution == 'octave':
        # Full Octave Band Edges: +/- 1/2 octave (2^(1/2) approx 1.414)
        lower_edge = 2**(-1.0/2.0) # ~0.707
        upper_edge = 2**(1.0/2.0)  # ~1.414
        
        # Passband Tolerance
        plt.fill_between([lower_edge, upper_edge], -0.3, 0.3, color='green', alpha=0.2, label='Passband (+/- 0.3dB)')
        
        # For Octave bands, the stopband requirements are wider.
        # Usually checking 1 octave away from the EDGE (factor of 2 from edge).
        # Edge * 2 = 1.414 * 2 = 2.828
        # Edge / 2 = 0.707 / 2 = 0.353
        plt.plot([0.353, 0.353], [-5, -70], 'r--', linewidth=1, label='Stopband Limit')
        plt.plot([2.828, 2.828], [-5, -70], 'r--', linewidth=1)

    plt.title(f"ANSI Compliance Check: Fc={fc}Hz @ Fs={fs}Hz\nResolution: {resolution}, Order={order}", fontsize=14)
    plt.xlabel("Normalized Frequency ($f / f_c$)", fontsize=12)
    plt.ylabel("Magnitude (dB)", fontsize=12)
    plt.grid(True, which="both", linestyle='-', alpha=0.6)
    
    # Zoom in on the relevant area
    plt.xlim(0.1, 10)
    plt.ylim(-90, 5)
    
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    # Test Case 1: Standard 1kHz 1/3 Octave at 48kHz
    print("Plotting 1kHz 1/3 Octave filter...")
    plot_mask_check(fc=1000, fs=48000, resolution='third')
    
    # Test Case 2: Extreme Low Freq (20Hz) 1/3 Octave at High Sample Rate
    print("Plotting 20Hz 1/3 Octave filter at 96kHz...")
    plot_mask_check(fc=20, fs=96000, resolution='third')

    # Test Case 3: Full Octave Filter
    print("Plotting 1kHz Full Octave filter...")
    plot_mask_check(fc=1000, fs=48000, resolution='octave')