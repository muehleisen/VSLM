# python2/vslm/filters/ansi.py
import numpy as np
import scipy.signal

def get_ansi_center_frequencies(resolution='octave', base=10):
    """Returns exact Center Frequencies (Fc) based on ANSI S1.11-2004."""
    f_ref = 1000.0
    
    if resolution == 'octave':
        x = np.arange(-6, 5) 
        if base == 10:
            fm = f_ref * (10.0 ** (3.0 * x / 10.0))
        else:
            fm = f_ref * (2.0 ** x)
            
    elif resolution == 'third':
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
    """
    if resolution == 'octave':
        bandwidth_factor = 2**(1.0 / 2.0)
    else:
        bandwidth_factor = 2**(1.0 / 6.0)
        
    f_lower_ansi = fc / bandwidth_factor
    f_upper_ansi = fc * bandwidth_factor
    
    target_attenuation_db = 0.05 
    term = (10**(target_attenuation_db / 10.0)) - 1
    alpha = term ** (1.0 / (2.0 * order))
    
    f_lower_design = f_lower_ansi * alpha
    f_upper_design = f_upper_ansi / alpha
    
    nyquist = fs / 2.0
    
    if f_upper_design >= nyquist * 0.99:
        f_upper_design = nyquist * 0.99
        if f_lower_design >= f_upper_design:
            raise ValueError(f"Band centered at {fc:.1f} Hz is too close to Nyquist ({nyquist} Hz)")

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
        self.filters = [] 
        
        all_centers = get_ansi_center_frequencies(resolution, base=10)
        
        if resolution == 'octave':
            factor = 2**(1.0/2.0)
        else:
            factor = 2**(1.0/6.0)
            
        cutoff_limit = (fs / 2.0) / factor * 0.95
        valid_centers = all_centers[all_centers < cutoff_limit]
        
        self.frequencies = valid_centers
        
        for fc in valid_centers:
            try:
                sos = design_compliant_sos(fc, fs, resolution, order)
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

    def initialize_state(self, chunk_data):
        """
        Seeds all filters in the bank to minimize transient glitches.
        Runs forward-backward on the provided chunk.
        """
        for band in self.filters:
            # 1. Forward pass
            zi_init = np.zeros_like(band['zi'])
            _, zi_fwd = scipy.signal.sosfilt(band['sos'], chunk_data, zi=zi_init)
            
            # 2. Backward pass
            _, zi_bwd = scipy.signal.sosfilt(band['sos'], chunk_data[::-1], zi=zi_fwd)
            
            # 3. Update state
            band['zi'] = zi_bwd

    def process_chunk(self, chunk_data):
        """
        Processes a chunk of audio through the entire filter bank.
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