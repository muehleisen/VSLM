# vslm/core.py
import numpy as np
import scipy.signal
import soundfile as sf
from numba import jit
import warnings
from .filters import A_WEIGHTING, C_WEIGHTING, SUPPORTED_FS

# --- Numba Optimized Routine for Impulse Meter ---
@jit(nopython=True)
def _apply_impulse_decay(signal_sq, fs):
    """
    Applies the rising/falling time constants for Impulse metering.
    Ported from vslm.m logic (lines 411-430).
    Using Numba for C-like loop performance.
    """
    tau_r = 0.035  # 35ms rising
    tau_f = 1.5    # 1.5s falling
    
    # Pre-calculate decay coefficients
    alpha_r = 1.0 - np.exp(-1 / (fs * tau_r))
    alpha_f = 1.0 - np.exp(-1 / (fs * tau_f))
    
    # Initialize output array
    n_samples = len(signal_sq)
    output = np.zeros(n_samples, dtype=np.float64)
    
    current_val = 0.0
    
    for i in range(n_samples):
        # The detector logic:
        # If the new signal is higher than current detector level -> Fast Rise
        # If the new signal is lower -> Slow Fall
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
        self.cal_factor = 1.0  # Linear scale factor
        self.filename = "None Loaded"

    def load_file(self, filepath):
        """Loads a wav file using soundfile (supports 24-bit)."""
        data, fs = sf.read(filepath)
        
        # VSLM supports Mono only (Channel 1)
        if data.ndim > 1:
            data = data[:, 0]
            
        # Enforce supported sample rates for ANSI compliance
        if fs not in SUPPORTED_FS:
             # In a real app, you might resample here, but for strict compliance:
             warnings.warn(f"Sample rate {fs} not in supported list {SUPPORTED_FS}. "
                           "Filters may be undefined.")
        
        self.audio_data = data
        self.fs = fs
        self.filename = filepath
        return True

    def set_calibration(self, db_level, ref_rms=None):
        """
        Sets calibration factor based on a known dB level.
        If ref_rms is not provided, calculates it from the currently loaded file 
        (assuming the loaded file IS the cal tone).
        """
        if ref_rms is None:
            if self.audio_data is None:
                raise ValueError("No audio loaded to calibrate against.")
            ref_rms = np.sqrt(np.mean(self.audio_data**2))
            
        # Calculate required linear pressure (Pascal) for the target dB
        # Reference pressure is 20 microPascals
        ref_pressure_pa = 20e-6
        target_pa = ref_pressure_pa * (10**(db_level / 20.0))
        
        self.cal_factor = target_pa / ref_rms
        return self.cal_factor

    def apply_weighting_filter(self, data, weighting='A'):
        """Applies A or C weighting using the golden coefficients."""
        if weighting == 'Z' or weighting == 'Flat':
            return data
            
        if self.fs not in SUPPORTED_FS:
            raise ValueError(f"Fs {self.fs} not supported for {weighting}-weighting.")

        if weighting == 'A':
            coeffs = A_WEIGHTING[self.fs]
        elif weighting == 'C':
            coeffs = C_WEIGHTING[self.fs]
        else:
            raise ValueError("Unknown weighting type")
            
        b, a = coeffs
        return scipy.signal.lfilter(b, a, data)

    def calculate_lp(self, weighting='A', speed='Slow', dt=0.1):
        """
        Calculates the Sound Pressure Level profile (Lp vs Time).
        
        Args:
            weighting: 'A', 'C', or 'Z'
            speed: 'Fast', 'Slow', or 'Impulse'
            dt: Output time resolution in seconds (e.g., 0.1s for graph points)
        
        Returns:
            time_axis (np.array): Time points
            lp_db (np.array): SPL in dB
        """
        # 1. Apply Calibration
        calibrated_data = self.audio_data * self.cal_factor
        
        # 2. Apply Frequency Weighting
        weighted_data = self.apply_weighting_filter(calibrated_data, weighting)
        
        # 3. Square the signal (Pressure^2)
        signal_sq = weighted_data**2
        
        # 4. Apply Time Weighting (Detector Physics)
        if speed == 'Impulse':
            # Use the Numba-optimized routine for the non-linear impulse detector
            detected_sq = _apply_impulse_decay(signal_sq, self.fs)
        else:
            # Use standard linear IIR for Fast/Slow
            if speed == 'Fast':
                tau = 0.125
            else: # Slow
                tau = 1.0
                
            # Filter coefficient calculation (Single pole lowpass)
            # b1 = 1 - exp(-1/(fs*tau))
            b1 = 1.0 - np.exp(-1.0 / (self.fs * tau))
            a1 = b1 - 1.0
            
            # lfilter: y[n] = b*x[n] - a*y[n-1]
            # SciPy uses: a[0]*y[n] = b[0]*x[n] + b[1]*x[n-1] - a[1]*y[n-1]...
            # Equivalent to the MATLAB loop:
            b_poly = [b1]
            a_poly = [1.0, -1.0 * (1.0 - b1)] # equivalent to a2 = b1 - 1 in matlab code structure
            
            detected_sq = scipy.signal.lfilter(b_poly, a_poly, signal_sq)

        # 5. Decimate/Resample for Plotting (Simulation of "dt" spacing)
        # We average or pick samples every 'dt' seconds. 
        # The MATLAB code simply picks the point: p2((Klast+1)...) = data2(J)
        step_samples = int(self.fs * dt)
        indices = np.arange(0, len(detected_sq), step_samples)
        
        # Safety check for short files
        if len(indices) == 0:
            return np.array([0]), np.array([0])
            
        downsampled_sq = detected_sq[indices]
        
        # 6. Convert to dB
        # Ref = (20e-6)^2 = 4e-10
        # But we already calibrated to Pascals.
        # So dB = 10 * log10(p^2 / p_ref^2)
        # Note: In the MATLAB code, they use 94dB reference logic sometimes, 
        # but if we calibrate to Pascals, we use standard physics:
        ref_pressure_sq = (20e-6)**2
        
        # Avoid log(0)
        downsampled_sq[downsampled_sq < 1e-30] = 1e-30
        lp_db = 10 * np.log10(downsampled_sq / ref_pressure_sq)
        
        time_axis = indices / self.fs
        
        return time_axis, lp_db

    def calculate_leq(self, weighting='A', duration=None):
        """
        Calculates the Equivalent Continuous Sound Level (Leq).
        """
        calibrated_data = self.audio_data * self.cal_factor
        weighted_data = self.apply_weighting_filter(calibrated_data, weighting)
        
        # Leq is strictly the mean of the squared pressure
        mean_sq = np.mean(weighted_data**2)
        ref_pressure_sq = (20e-6)**2
        
        if mean_sq <= 0:
            return -np.inf
            
        leq = 10 * np.log10(mean_sq / ref_pressure_sq)
        return leq

    def calculate_octave_bands(self, weighting='Z', resolution='octave'):
        """
        Calculates Octave or 1/3 Octave bands using FFT method (VSLM "High Res FFT").
        Ref: vslm.m analyze_bandfft
        """
        calibrated_data = self.audio_data * self.cal_factor
        
        # Apply weighting if needed (usually bands are Z-weighted, but VSLM allows A-weighted bands)
        data = self.apply_weighting_filter(calibrated_data, weighting)
        
        Nfft = 65536 # 2^16, matching VSLM default
        
        # Compute Power Spectrum
        f, Pxx = scipy.signal.welch(data, fs=self.fs, window='boxcar', 
                                    nperseg=Nfft, noverlap=0, scaling='spectrum')
        
        # Define Band Centers
        # Base 2 formulation (1000 * 2^((n-30)/3)) generally
        if resolution == 'octave':
            # Standard Octave Centers: 16, 31.5, 63, 125, 250, 500, 1k, 2k, 4k, 8k, 16k
            center_freqs = np.array([16, 31.5, 63, 125, 250, 500, 1000, 2000, 4000, 8000, 16000])
            bandwidth_factor = 2**(1/2) # Full octave
        else:
            # 1/3 Octave
            # Generate from 12.5Hz to 20kHz
            # Using formula f_c = 1000 * 2^(n/3)
            n_min = -19 # approx 12.5
            n_max = 13  # approx 20k
            n = np.arange(n_min, n_max + 1)
            center_freqs = 1000 * (2**(n/3.0))
            bandwidth_factor = 2**(1/6)
            
        # Sum bins for each band
        band_levels = []
        valid_centers = []
        
        for fc in center_freqs:
            if fc > self.fs / 2:
                break
                
            fl = fc / bandwidth_factor
            fu = fc * bandwidth_factor
            
            # Find indices in FFT
            idx_start = np.searchsorted(f, fl)
            idx_end = np.searchsorted(f, fu)
            
            if idx_start == idx_end:
                # Bin too coarse for this low freq, take nearest bin
                power_sum = Pxx[idx_start] if idx_start < len(Pxx) else 0
            else:
                power_sum = np.sum(Pxx[idx_start:idx_end])
                
            ref_pressure_sq = (20e-6)**2
            if power_sum <= 0:
                lvl = 0
            else:
                lvl = 10 * np.log10(power_sum / ref_pressure_sq)
            
            band_levels.append(lvl)
            valid_centers.append(fc)
            
        return np.array(valid_centers), np.array(band_levels)