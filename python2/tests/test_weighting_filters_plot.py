# python2/tests/test_weighting_filter_plot.py
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal

# --- Path Hack to import from vslm package ---
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Import from the new location
from vslm.filters.weighting import WeightingFilter, SUPPORTED_FS

def get_ideal_weighting(freqs, type='A'):
    """Computes the theoretical ANSI S1.42 weighting curve."""
    f2 = freqs**2
    f1 = 20.598997
    f2_const = 107.65265
    f3 = 737.86223
    f4 = 12194.217
    
    if type == 'A':
        num = (f4**2) * (f2**2)
        den1 = (f2 + f1**2)
        den2 = np.sqrt((f2 + f2_const**2) * (f2 + f3**2))
        den3 = (f2 + f4**2)
        gain = num / (den1 * den2 * den3)
        # Normalize to 0dB at 1kHz
        ref_f = 1000.0
        ref_val = ((f4**2) * (ref_f**4)) / (
            ((ref_f**2) + f1**2) * np.sqrt(((ref_f**2) + f2_const**2) * ((ref_f**2) + f3**2)) * ((ref_f**2) + f4**2)
        )
        return 20 * np.log10(gain / ref_val)
        
    elif type == 'C':
        num = (f4**2) * f2
        den = (f2 + f1**2) * (f2 + f4**2)
        gain = num / den
        ref_f = 1000.0
        ref_val = ((f4**2) * (ref_f**2)) / (((ref_f**2) + f1**2) * ((ref_f**2) + f4**2))
        return 20 * np.log10(gain / ref_val)

def get_class1_tolerances(w_type='A'):
    """Returns Class 1 Tolerance Masks."""
    freqs = np.array([10, 12.5, 16, 20, 25, 31.5, 40, 50, 63, 80, 
                      100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 
                      1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 
                      8000, 10000, 12500, 16000, 20000])
    plus = []
    minus = []
    for f in freqs:
        if f < 25:
            p, m = 2.5, -np.inf 
        elif f < 100:
            p, m = 1.5, -1.5
        elif f <= 5000:
            p, m = 1.1, -1.1 
        elif f <= 16000:
            p, m = 2.5, -2.5
        else:
            p, m = 3.0, -5.0
        plus.append(p)
        minus.append(m)
    return freqs, np.array(plus), np.array(minus)

def plot_weighting(w_type='A'):
    plt.figure(figsize=(12, 8))
    
    # 1. Ideal Curve
    f_ideal = np.logspace(1, 5.3, 1000) 
    ideal_db = get_ideal_weighting(f_ideal, w_type)
    plt.semilogx(f_ideal, ideal_db, 'k--', linewidth=2, zorder=5, label=f'Ideal {w_type} (ANSI S1.42)')
    
    # 2. Tolerance Bands
    tol_f, tol_plus, tol_minus = get_class1_tolerances(w_type)
    ideal_at_tol = get_ideal_weighting(tol_f, w_type)
    upper_mask = ideal_at_tol + tol_plus
    lower_mask = ideal_at_tol + tol_minus
    lower_mask[lower_mask == -np.inf] = -200
    
    plt.fill_between(tol_f, lower_mask, upper_mask, color='green', alpha=0.2, zorder=1, label='Class 1 Tolerance')
    
    # 3. Dynamic SOS Filters for Supported Rates
    colors = plt.cm.plasma(np.linspace(0, 0.9, len(SUPPORTED_FS)))
    
    for i, fs in enumerate(SUPPORTED_FS):
        try:
            # Instantiate the filter
            wf = WeightingFilter(fs, w_type)
            
            # --- UPDATED: Access SOS coefficients and use sosfreqz ---
            sos = wf.sos
            w, h = scipy.signal.sosfreqz(sos, worN=16384, fs=fs)
            
            mag_db = 20 * np.log10(np.abs(h) + 1e-15)
            
            # Plot up to Nyquist
            mask = (w >= 10) & (w <= fs/2)
            
            plt.semilogx(w[mask], mag_db[mask], 
                         label=f'Fs={fs} Hz', 
                         color=colors[i], linewidth=1.5, alpha=0.8, zorder=10)

        except Exception as e:
            print(f"Skipping {fs}Hz: {e}")
            continue

    plt.title(f"{w_type}-Weighting Filter: MZT Hybrid SOS (High Accuracy)")
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Magnitude (dB)")
    plt.grid(True, which='both', linestyle='--', alpha=0.5)
    plt.legend(loc='lower center', ncol=4, fontsize='small')
    plt.xlim(10, 100000) 
    
    if w_type == 'A':
        plt.ylim(-80, 15) 
    else:
        plt.ylim(-20, 15)
        
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    print("Plotting A-Weighting...")
    plot_weighting('A')
    
    print("Plotting C-Weighting...")
    plot_weighting('C')