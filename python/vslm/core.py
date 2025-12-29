import numpy as np
import scipy.signal
import soundfile as sf
from numba import jit
import warnings

# Use the new constants from weighting_filters
from .weighting_filters import apply_weighting, MIN_SUPPORTED_FS, MAX_SUPPORTED_FS

# --- Numba Optimized Routine for Impulse Meter ---
@jit(nopython=True)
def _apply_impulse_decay(signal_sq, fs):
    tau_r = 0.035
    tau_f = 1.5
    alpha_r = 1.0 - np.exp(-1 / (fs * tau_r))
    alpha_f = 1.0 - np.exp(-1 / (fs * tau_f))
    
    n_samples = len(signal_sq)
    output = np.zeros(n_samples, dtype=np.float64)
    current_val = 0.0
    
    for i in range(n_samples):
        if signal_sq[i] > current_val:
            current_val = alpha_r * signal_sq[i] + (1 - alpha_r) * current_val
        else:
            current_val = alpha_f * signal_sq[i] + (1 - alpha_f) * current_val
        output[i] = current_val
        
    return output

class VSLMCore:
    def __init__(self):
        self.audio_data = None
        self.fs = 0
        self.cal_factor = 1.0
        self.filename = "None Loaded"

    def load_file(self, filepath):
        """Loads a wav file using soundfile (supports 24-bit)."""
        data, fs = sf.read(filepath)
        
        if data.ndim > 1:
            data = data[:, 0]
            
        # Updated Range Validation
        if not (MIN_SUPPORTED_FS <= fs <= MAX_SUPPORTED_FS):
             warnings.warn(f"Sample rate {fs} Hz is outside the standard range "
                           f"({MIN_SUPPORTED_FS}-{MAX_SUPPORTED_FS} Hz). "
                           "Filter accuracy is not guaranteed.")
        
        self.audio_data = data
        self.fs = fs
        self.filename = filepath
        return True

    def set_calibration(self, db_level, ref_rms=None):
        if ref_rms is None:
            if self.audio_data is None:
                raise ValueError("No audio loaded to calibrate against.")
            ref_rms = np.sqrt(np.mean(self.audio_data**2))
            
        ref_pressure_pa = 20e-6
        target_pa = ref_pressure_pa * (10**(db_level / 20.0))
        self.cal_factor = target_pa / ref_rms
        return self.cal_factor

    def apply_weighting_filter(self, data, weighting='A'):
        return apply_weighting(data, self.fs, weighting)

    def calculate_lp(self, weighting='A', speed='Slow', dt=0.1):
        calibrated_data = self.audio_data * self.cal_factor
        weighted_data = self.apply_weighting_filter(calibrated_data, weighting)
        signal_sq = weighted_data**2
        
        if speed == 'Impulse':
            detected_sq = _apply_impulse_decay(signal_sq, self.fs)
        else:
            if speed == 'Fast':
                tau = 0.125
            else: # Slow
                tau = 1.0
            
            b1 = 1.0 - np.exp(-1.0 / (self.fs * tau))
            b_poly = [b1]
            a_poly = [1.0, -1.0 * (1.0 - b1)]
            detected_sq = scipy.signal.lfilter(b_poly, a_poly, signal_sq)

        step_samples = int(self.fs * dt)
        indices = np.arange(0, len(detected_sq), step_samples)
        
        if len(indices) == 0:
            return np.array([0]), np.array([0])
            
        downsampled_sq = detected_sq[indices]
        ref_pressure_sq = (20e-6)**2
        
        downsampled_sq[downsampled_sq < 1e-30] = 1e-30
        lp_db = 10 * np.log10(downsampled_sq / ref_pressure_sq)
        time_axis = indices / self.fs
        
        return time_axis, lp_db

    def calculate_leq(self, weighting='A', duration=None):
        calibrated_data = self.audio_data * self.cal_factor
        weighted_data = self.apply_weighting_filter(calibrated_data, weighting)
        mean_sq = np.mean(weighted_data**2)
        ref_pressure_sq = (20e-6)**2
        if mean_sq <= 0: return -np.inf
        return 10 * np.log10(mean_sq / ref_pressure_sq)