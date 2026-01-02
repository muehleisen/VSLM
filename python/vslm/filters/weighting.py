# python2/vslm/filters/weighting.py
import numpy as np
import scipy.signal
import scipy.optimize

# Standard VSLM sampling rates (for validation)
SUPPORTED_FS = [22050, 44100, 48000, 96000, 88200, 176400, 192000]

# ISO Preferred Frequencies for Optimization Targets
ISO_FREQS = np.array([
    10, 12.5, 16, 20, 25, 31.5, 40, 50, 63, 80, 
    100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 
    1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 
    8000, 10000, 12500, 16000, 20000
])

# --- Optimization Helper Functions ---

def _get_ideal_response(freqs, w_type):
    """Calculates theoretical ANSI dB gain."""
    f2 = freqs**2
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
        nr = (f4**2) * (r2**2)
        dr = (r2 + f1**2) * np.sqrt((r2 + f2_c**2) * (r2 + f3**2)) * (r2 + f4**2)
    else:
        nr = (f4**2) * r2
        dr = (r2 + f1**2) * (r2 + f4**2)
        
    return 20 * np.log10(gain / (nr/dr))

def _design_parametric_sos(f0, Q, gain_db, fs):
    """Generic Parametric EQ for fine-tuning."""
    if f0 <= 0 or f0 >= fs/2: return np.array([1., 0., 0., 1., 0., 0.])
    
    A = 10 ** (gain_db / 40.0)
    w0 = 2 * np.pi * f0 / fs
    alpha = np.sin(w0) / (2.0 * Q)
    
    b0 = 1 + alpha * A
    b1 = -2 * np.cos(w0)
    b2 = 1 - alpha * A
    a0 = 1 + alpha / A
    a1 = -2 * np.cos(w0)
    a2 = 1 - alpha / A
    
    return np.array([b0, b1, b2, a0, a1, a2]) / a0

def design_optimized_sos(fs, weighting_type='A'):
    """
    Designs an 8th-order weighting filter (4 Biquads) using Hybrid Transform.
    Uses MZT for high-frequency poles to prevent cramping near Nyquist.
    """
    
    # --- 1. Fixed Low-Frequency Sections (Bilinear) ---
    f1, f2, f3 = 20.598997, 107.65265, 737.86223
    w1, w2, w3 = [2 * np.pi * f for f in [f1, f2, f3]]
    
    fixed_sos_blocks = []
    if weighting_type == 'A':
        # H1: Highpass (w1)
        z, p, k = scipy.signal.bilinear_zpk([0, 0], [-w1, -w1], 1.0, fs)
        fixed_sos_blocks.append(scipy.signal.zpk2sos(z, p, k))
        # H2: Bandpass (w2, w3)
        z, p, k = scipy.signal.bilinear_zpk([0, 0], [-w2, -w3], 1.0, fs)
        fixed_sos_blocks.append(scipy.signal.zpk2sos(z, p, k))
    else: # C
        # H1: Highpass (w1)
        z, p, k = scipy.signal.bilinear_zpk([0, 0], [-w1, -w1], 1.0, fs)
        fixed_sos_blocks.append(scipy.signal.zpk2sos(z, p, k))

    fixed_sos = np.vstack(fixed_sos_blocks)
    
    # --- 2. Setup Targets ---
    limit_f = min(20000, 0.94 * fs / 2) # Check up to 94% of Nyquist
    valid_mask = ISO_FREQS <= limit_f
    target_freqs = ISO_FREQS[valid_mask]
    target_db = _get_ideal_response(target_freqs, weighting_type)
    
    # Add dense points > 4kHz to ensure smoothness
    if fs < 50000:
        dense_high = np.linspace(4000, limit_f, 25)
        target_freqs = np.concatenate((target_freqs, dense_high))
        target_db = np.concatenate((target_db, _get_ideal_response(dense_high, weighting_type)))
        idx = np.argsort(target_freqs)
        target_freqs = target_freqs[idx]
        target_db = target_db[idx]

    w_eval = 2 * np.pi * target_freqs / fs
    _, h_fixed = scipy.signal.sosfreqz(fixed_sos, worN=w_eval)
    fixed_resp = h_fixed

    # --- 3. Cost Function (Hybrid MZT + Minimax) ---
    def cost_func(params):
        f4_val, cp_f, cp_g, cp_q = params
        
        # A. High Freq Section (MZT Transform)
        w4 = 2 * np.pi * f4_val
        p_mzt = np.exp(-w4 / fs)
        k_mzt = (1 - p_mzt)**2
        sos_high = np.array([[k_mzt, 0, 0, 1, -2*p_mzt, p_mzt**2]])
        
        # B. Correction Biquad
        sos_corr = _design_parametric_sos(cp_f, cp_q, cp_g, fs)
        
        sos_var = np.vstack((sos_high, sos_corr))
        _, h_var = scipy.signal.sosfreqz(sos_var, worN=w_eval)
        h_total = fixed_resp * h_var
        
        # Normalize 1kHz
        idx_1k = np.argmin(np.abs(target_freqs - 1000))
        gain_1k = np.abs(h_total[idx_1k]) + 1e-15
        
        mag_db = 20 * np.log10(np.abs(h_total) / gain_1k + 1e-15)
        error = np.abs(mag_db - target_db)
        
        # Weighting
        weights = np.ones_like(target_freqs)
        weights[target_freqs > 1000] = 50
        
        return np.max(error * weights)

    # --- 4. Optimizer ---
    x0 = [12194.0, 15000.0, 0.0, 1.0]

    res = scipy.optimize.minimize(
        cost_func, x0, method='Nelder-Mead', tol=1e-5,
        options={'maxiter': 500}
    )
    
    best_f4, best_cf, best_cg, best_cq = res.x

    # --- 5. Final Construction ---
    w4_final = 2 * np.pi * best_f4
    p_mzt = np.exp(-w4_final / fs)
    k_mzt = (1 - p_mzt)**2
    sos_high = np.array([[k_mzt, 0, 0, 1, -2*p_mzt, p_mzt**2]])
    
    sos_corr = _design_parametric_sos(best_cf, best_cq, best_cg, fs)
    
    sos_final = np.vstack((fixed_sos, sos_high, sos_corr))
    
    # Final Normalization
    w_ref = 2 * np.pi * 1000.0 / fs
    _, h_ref = scipy.signal.sosfreqz(sos_final, worN=[w_ref])
    gain_corr = 1.0 / (np.abs(h_ref[0]) + 1e-15)
    sos_final[0, :3] *= gain_corr
    
    return sos_final

# --- Class Implementation ---

class WeightingFilter:
    """
    Stateful filter wrapper for streaming analysis.
    Uses dynamically generated optimized SOS filters (MZT Hybrid).
    """
    def __init__(self, fs, weighting_type='A'):
        self.fs = fs
        self.weighting_type = weighting_type.upper()
        self.sos = None
        self.zi = None
        
        # 'Z' or 'Flat' means no filtering
        if self.weighting_type in ['Z', 'FLAT', 'NONE']:
            self.passthrough = True
            return
        
        self.passthrough = False
        
        if fs not in SUPPORTED_FS:
             print(f"Warning: Non-standard sampling rate {fs} Hz. "
                   "Weighting accuracy depends on dynamic filter generation.")

        try:
            # Generate Coefficients
            self.sos = design_optimized_sos(fs, self.weighting_type)
            
            # Initialize state (zi)
            self.zi = scipy.signal.sosfilt_zi(self.sos)
            
        except Exception as e:
            raise ValueError(f"Failed to design {self.weighting_type}-weighting: {e}")

    def reset(self):
        """Clears the internal filter state."""
        if not self.passthrough:
            self.zi = scipy.signal.sosfilt_zi(self.sos)

    def initialize_state(self, chunk_data):
        """
        Seeds the filter state (zi) to minimize transient glitches.
        Method: Runs forward-backward on the chunk and uses the final backward state.
        
        Args:
            chunk_data (np.ndarray): The first block of audio to be analyzed.
        """
        if self.passthrough:
            return
            
        # 1. Forward pass (starting from zero state) -> gets state at end of chunk
        zi_init = np.zeros_like(self.zi)
        _, zi_fwd = scipy.signal.sosfilt(self.sos, chunk_data, zi=zi_init)
        
        # 2. Backward pass (starting from forward state) -> gets state at start of chunk
        _, zi_bwd = scipy.signal.sosfilt(self.sos, chunk_data[::-1], zi=zi_fwd)
        
        # 3. Set this "warmed up" state as the actual starting state
        self.zi = zi_bwd

    def process_chunk(self, chunk_data):
        """
        Filters a chunk of audio data, updating internal state.
        """
        if self.passthrough:
            return chunk_data
        
        filtered_data, self.zi = scipy.signal.sosfilt(
            self.sos, chunk_data, zi=self.zi
        )
        return filtered_data