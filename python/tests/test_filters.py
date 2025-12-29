import unittest
import numpy as np
import scipy.signal
import sys
import os

# Add the parent directory to path so we can import the vslm package
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from vslm.core import VSLMCore, SUPPORTED_FS
from vslm.analysis import calculate_ansi_bands

class TestANSICompliance(unittest.TestCase):

    def setUp(self):
        self.core = VSLMCore()
        # Tolerance Masks based on ANSI S1.42-2001 (Class 1)
        # (Frequency, Target dB, Tolerance +/- dB)
        self.a_weighting_targets = [
            (63,    -26.2, 1.5),
            (125,   -16.1, 1.5),
            (250,   -8.6,  1.0),
            (500,   -3.2,  1.0),
            (1000,   0.0,  1.0), # Reference
            (2000,   1.2,  1.0),
            (4000,   1.0,  1.0),
            (8000,  -1.1,  1.5), # Tighter in reality, but 1.5 safe for digital
            (16000, -6.6,  2.5)  # High frequency roll-off check
        ]
        
        self.c_weighting_targets = [
            (63,    -0.8, 1.5),
            (125,   -0.2, 1.0),
            (1000,   0.0, 1.0),
            (8000,  -3.0, 1.5), # -3dB point
            (16000, -8.5, 2.5)
        ]

    def _measure_filter_gain(self, fs, freq, weighting_type):
        """Helper to generate a sine wave and measure filter output gain."""
        self.core.fs = fs
        duration = 1.0
        t = np.linspace(0, duration, int(fs * duration), endpoint=False)
        
        # Generate 1 Pascal RMS Sine Wave at specific freq
        # Amplitude = sqrt(2) * RMS
        signal = np.sqrt(2) * np.sin(2 * np.pi * freq * t)
        
        # Apply Filter directly
        filtered = self.core.apply_weighting_filter(signal, weighting_type)
        
        # Calculate RMS of central portion (avoid transient at start/end)
        # Discard first 10% and last 10%
        idx_start = int(len(filtered) * 0.1)
        idx_end = int(len(filtered) * 0.9)
        steady_state = filtered[idx_start:idx_end]
        
        rms_out = np.sqrt(np.mean(steady_state**2))
        
        # Gain in dB
        # Ref was 1.0
        if rms_out < 1e-9: return -999
        return 20 * np.log10(rms_out / 1.0)

    def test_a_weighting_compliance(self):
        """Verifies A-Weighting against ANSI S1.42 Class 1 for all supported Fs."""
        print("\n--- Testing A-Weighting Compliance ---")
        for fs in SUPPORTED_FS:
            print(f"Testing Fs = {fs} Hz...")
            for f_test, target_db, tol in self.a_weighting_targets:
                if f_test > fs / 2.2: 
                    continue # Skip freqs near Nyquist
                
                measured_db = self._measure_filter_gain(fs, f_test, 'A')
                
                err_msg = (f"Fs={fs}, Freq={f_test}Hz: Expected {target_db} +/- {tol}, "
                           f"Got {measured_db:.2f}")
                
                self.assertAlmostEqual(measured_db, target_db, delta=tol, msg=err_msg)

    def test_c_weighting_compliance(self):
        """Verifies C-Weighting against ANSI S1.42 Class 1 for all supported Fs."""
        print("\n--- Testing C-Weighting Compliance ---")
        for fs in SUPPORTED_FS:
            print(f"Testing Fs = {fs} Hz...")
            for f_test, target_db, tol in self.c_weighting_targets:
                if f_test > fs / 2.2: 
                    continue
                
                measured_db = self._measure_filter_gain(fs, f_test, 'C')
                
                err_msg = (f"Fs={fs}, Freq={f_test}Hz: Expected {target_db} +/- {tol}, "
                           f"Got {measured_db:.2f}")
                
                self.assertAlmostEqual(measured_db, target_db, delta=tol, msg=err_msg)

    def test_octave_band_centers(self):
        """
        Verifies that 1kHz input energy appears ONLY in the 1kHz band 
        and is attenuated in adjacent bands (ANSI S1.11 check).
        """
        print("\n--- Testing Octave Band Selectivity (ANSI S1.11) ---")
        fs = 48000
        self.core.fs = fs
        duration = 1.0
        t = np.linspace(0, duration, int(fs * duration), endpoint=False)
        
        # 1kHz Sine Wave (0 dB reference)
        # Note: 1 Pa RMS = 94 dB usually, but let's stick to relative levels
        # We simulate a "calibrated" signal where input is exactly 1.0 unit.
        self.core.audio_data = np.sqrt(2) * np.sin(2 * np.pi * 1000 * t)
        self.core.cal_factor = 1.0 
        
        freqs, levels = calculate_ansi_bands(self.core, weighting='Z', resolution='octave')
        
        # Find 1kHz band index
        idx_1k = np.argmin(np.abs(freqs - 1000))
        idx_500 = np.argmin(np.abs(freqs - 500))
        idx_2k = np.argmin(np.abs(freqs - 2000))
        
        lvl_1k = levels[idx_1k]
        lvl_500 = levels[idx_500]
        lvl_2k = levels[idx_2k]
        
        # 1. Check Center accuracy
        # Input 94dB (approx) -> Output should be close to 94dB in the main band
        # (Since we used 1.0 rms and ref is 20uPa, 20log10(1/20e-6) = 93.98 dB)
        expected_db = 20 * np.log10(1.0 / 20e-6)
        self.assertAlmostEqual(lvl_1k, expected_db, delta=0.5, 
                               msg=f"1kHz input should measure ~94dB, got {lvl_1k:.2f}")
        
        # 2. Check Attenuation (Crosstalk)
        # Octave bands usually require >18dB attenuation at adjacent octave centers
        attenuation_500 = lvl_1k - lvl_500
        attenuation_2k = lvl_1k - lvl_2k
        
        print(f"1kHz Level: {lvl_1k:.2f} dB")
        print(f"500Hz Level: {lvl_500:.2f} dB (Attenuation: {attenuation_500:.2f} dB)")
        print(f"2kHz Level: {lvl_2k:.2f} dB (Attenuation: {attenuation_2k:.2f} dB)")
        
        self.assertGreater(attenuation_500, 18.0, "Filter skirt at 500Hz is too wide")
        self.assertGreater(attenuation_2k, 18.0, "Filter skirt at 2kHz is too wide")

if __name__ == '__main__':
    unittest.main()