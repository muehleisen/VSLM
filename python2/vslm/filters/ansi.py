# vslm/filters/ansi.py
import numpy as np
import scipy.signal

def get_ansi_center_frequencies(resolution='octave', base=10):
    """
    Returns exact Center Frequencies (Fc) based on ANSI S1.11-2004.
    Matches standard base-10 definitions by default.
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

def design_compliant_sos(fc, fs, resolution='third', order=24):
    """
    Designs a high-order Butterworth filter compliant with ANSI S1.11 Class 1.
    
    Features:
    - Uses Order N=24 (48th order total) for steep skirts.
    - Applies 'Bandwidth Correction' to ensure the passband remains flat
      at the edges despite the high Q required.
    """
    # 1. Determine Exact ANSI Band Edges
    if resolution == 'octave':
        bandwidth_factor = 2**(1.0 / 2.0) # Sqrt(2)
    else:
        bandwidth_factor = 2**(1.0 / 6.0) # Sixth-root(2)
        
    f_lower_ansi = fc / bandwidth_factor
    f_upper_ansi = fc * bandwidth_factor
    
    # 2. Bandwidth Correction Strategy
    # Force the filter to be -0.05 dB (almost flat) at the ANSI edges.
    # Without this, a high-order Butterworth would be -3dB at the edges,
    # which violates the flatness requirement of the standard.
    target_attenuation_db = 0.05 
    
    # Calculate correction factor alpha
    # Formula: ratio = ( 10^(target/10) - 1 )^(1/2N)
    term = (10**(target_attenuation_db / 10.0)) - 1
    alpha = term ** (1.0 / (2.0 * order))
    
    # Apply correction to widen the "Design Frequencies" so the "ANSI Frequencies"
    # land at the -0.05dB point instead of the -3dB point.
    f_lower_design = f_lower_ansi * alpha
    f_upper_design = f_upper_ansi / alpha
    
    # 3. Nyquist Check & Clamp
    nyquist = fs / 2.0
    
    # If the upper design edge hits Nyquist, the filter will crash.
    # We clamp it safely below Nyquist.
    if f_upper_design >= nyquist * 0.99:
        f_upper_design = nyquist * 0.99
        # If the lower edge is also too high, this band is invalid
        if f_lower_design >= f_upper_design:
            raise ValueError(f"Band centered at {fc:.1f} Hz is too close to Nyquist ({nyquist} Hz)")

    # 4. Design the Filter (SOS output is critical for stability at Order 24)
    sos = scipy.signal.butter(
        N=order, 
        Wn=[f_lower_design, f_upper_design], 
        btype='bandpass', 
        fs=fs, 
        output='sos'
    )
    
    return sos

class OctaveFilterBank:
    """
    Manages a bank of stateful high-order bandpass filters.
    """
    def __init__(self, fs, resolution='octave', order=24):
        self.fs = fs
        self.resolution = resolution
        self.filters = [] # List of dicts: {'fc': float, 'sos': array, 'zi': array}
        
        # 1. Get Nominal Center Frequencies (Base 10 matches ANSI best)
        all_centers = get_ansi_center_frequencies(resolution, base=10)
        
        # 2. Filter out frequencies too high for this sampling rate
        # We check if the *upper edge* of the band is feasible.
        if resolution == 'octave':
            factor = 2**(1.0/2.0)
        else:
            factor = 2**(1.0/6.0)
            
        # Keep bands where the theoretical upper edge is below 95% of Nyquist
        cutoff_limit = (fs / 2.0) / factor * 0.95
        valid_centers = all_centers[all_centers < cutoff_limit]
        
        self.frequencies = valid_centers
        
        # 3. Design and Initialize Filters
        for fc in valid_centers:
            try:
                sos = design_compliant_sos(fc, fs, resolution, order)
                
                # Initialize state (zi)
                # sosfilt_zi needs shape (n_sections, 2)
                zi = scipy.signal.sosfilt_zi(sos)
                
                self.filters.append({
                    'fc': fc,
                    'sos': sos,
                    'zi': zi
                })
            except Exception as e:
                print(f"Warning: Could not design filter for {fc:.1f} Hz: {e}")

    def reset(self):
        """Resets the state of all filters in the bank."""
        for band in self.filters:
            band['zi'] = scipy.signal.sosfilt_zi(band['sos'])

    def process_chunk(self, chunk_data):
        """
        Processes a chunk of audio through the entire filter bank.
        
        Args:
            chunk_data (np.ndarray): 1D array of audio samples.
            
        Returns:
            np.ndarray: 2D array of filtered data [n_samples, n_bands].
        """
        n_samples = len(chunk_data)
        n_bands = len(self.filters)
        
        output = np.zeros((n_samples, n_bands), dtype=chunk_data.dtype)
        
        for i, band in enumerate(self.filters):
            filtered_signal, new_zi = scipy.signal.sosfilt(
                band['sos'], 
                chunk_data, 
                zi=band['zi']
            )
            band['zi'] = new_zi
            output[:, i] = filtered_signal
            
        return output