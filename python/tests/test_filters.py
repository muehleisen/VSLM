import unittest
import numpy as np
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from vslm.core import VSLMCore
# Use STANDARD_RATES for testing
from vslm.weighting_filters import STANDARD_RATES 
from vslm.analysis import calculate_ansi_bands

class TestANSICompliance(unittest.TestCase):

    def setUp(self):
        self.core = VSLMCore()
        # ANSI S1.42 Class 1 Tolerances
        self.a_weighting_targets = [
            (63,    -26.2, 1.5),
            (125,   -16.1, 1.5),
            (250,   -8.6,  1.0),
            (500,   -3.2,  1.0),
            (1000,   0.0,  1.0),
            (2000,   1.2,  1.0),
            (4000,   1.0,  1.0),
            (8000,  -1.1,  1.5), 
            (16000, -6.6,  2.5)  
        ]
        self.c_weighting_targets = [
            (63,    -0.8, 1.5),
            (125,   -0.2, 1.0),
            (1000,   0.0, 1.0),
            (8000,  -3.0, 1.5),
            (16000, -8.5, 2.5)
        ]

    def _measure_filter_gain(self, fs, freq, weighting_type):
        self.core.fs = fs
        duration = 1.0
        t = np.linspace(0, duration, int(fs * duration), endpoint=False)
        signal = np.sqrt(2) * np.sin(2 * np.pi * freq * t)
        filtered = self.core.apply_weighting_filter(signal, weighting_type)
        
        idx_start = int(len(filtered) * 0.25)
        idx_end = int(len(filtered) * 0.75)
        steady_state = filtered[idx_start:idx_end]
        
        rms_out = np.sqrt(np.mean(steady_state**2))
        if rms_out < 1e-12: return -999.0
        return 20 * np.log10(rms_out / 1.0)

    def test_a_weighting_compliance(self):
        print("\n--- Testing A-Weighting Compliance (Dynamic) ---")
        # Iterate over all standard rates + high res
        for fs in STANDARD_RATES:
            print(f"Testing Fs = {fs} Hz...")
            for f_test, target_db, tol in self.a_weighting_targets:
                if f_test > fs / 2.2: continue 
                measured_db = self._measure_filter_gain(fs, f_test, 'A')
                err_msg = (f"Fs={fs}, Freq={f_test}Hz: Expected {target_db} +/- {tol}, Got {measured_db:.2f}")
                self.assertAlmostEqual(measured_db, target_db, delta=tol, msg=err_msg)

    def test_c_weighting_compliance(self):
        print("\n--- Testing C-Weighting Compliance (Dynamic) ---")
        for fs in STANDARD_RATES:
            print(f"Testing Fs = {fs} Hz...")
            for f_test, target_db, tol in self.c_weighting_targets:
                if f_test > fs / 2.2: continue
                measured_db = self._measure_filter_gain(fs, f_test, 'C')
                err_msg = (f"Fs={fs}, Freq={f_test}Hz: Expected {target_db} +/- {tol}, Got {measured_db:.2f}")
                self.assertAlmostEqual(measured_db, target_db, delta=tol, msg=err_msg)

    def test_octave_band_selectivity(self):
        print("\n--- Testing Octave Band Selectivity ---")
        fs = 48000
        self.core.fs = fs
        duration = 1.0
        t = np.linspace(0, duration, int(fs * duration), endpoint=False)
        self.core.audio_data = np.sqrt(2) * np.sin(2 * np.pi * 1000 * t)
        self.core.cal_factor = 1.0 
        
        freqs, levels = calculate_ansi_bands(self.core, weighting='Z', resolution='octave')
        idx_1k = np.argmin(np.abs(freqs - 1000))
        
        lvl_1k = levels[idx_1k]
        lvl_lower = levels[idx_1k - 1]
        
        expected_db = 20 * np.log10(1.0 / 20e-6)
        self.assertAlmostEqual(lvl_1k, expected_db, delta=0.5)
        
        attenuation_lower = lvl_1k - lvl_lower
        self.assertGreater(attenuation_lower, 60.0)

if __name__ == '__main__':
    unittest.main()