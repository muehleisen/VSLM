import unittest
import numpy as np
import sys
import os

# Add the parent directory to path so we can import the vslm package
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from vslm.filters.weighting import WeightingFilter, SUPPORTED_FS
from vslm.filters.ansi import OctaveFilterBank

class TestWeightingFilter(unittest.TestCase):
    def setUp(self):
        # ANSI S1.42 Class 1 Tolerances (simplified for key frequencies)
        self.a_targets = [
            (63,    -26.2, 1.5),
            (125,   -16.1, 1.5),
            (250,   -8.6,  1.0),
            (1000,   0.0,  1.0),
            (2000,   1.2,  1.0),
            (4000,   1.0,  1.0),
            (8000,  -1.1,  1.5), 
            (16000, -6.6,  2.5)  
        ]
        self.c_targets = [
            (63,    -0.8, 1.5),
            (125,   -0.2, 1.0),
            (1000,   0.0, 1.0),
            (8000,  -3.0, 1.5)
        ]

    def _measure_gain(self, fs, freq, w_type):
        """Generates a tone, filters it, and calculates steady-state gain."""
        duration = 1.0
        t = np.linspace(0, duration, int(fs * duration), endpoint=False)
        
        # 1.0 Vrms input (amplitude = sqrt(2))
        signal = np.sqrt(2) * np.sin(2 * np.pi * freq * t)
        
        # Instantiate the new Filter Class
        wf = WeightingFilter(fs, w_type)
        
        # Process in one go for this specific test (steady state check)
        filtered = wf.process_chunk(signal)
        
        # Discard first 25% to ignore filter transient (settling time)
        idx_start = int(len(filtered) * 0.25)
        steady_state = filtered[idx_start:]
        
        rms_out = np.sqrt(np.mean(steady_state**2))
        
        # Avoid log(0)
        if rms_out < 1e-12: return -999.0
        
        # Gain in dB
        return 20 * np.log10(rms_out / 1.0)

    def test_a_weighting_compliance(self):
        print("\n--- Testing A-Weighting Compliance ---")
        for fs in SUPPORTED_FS:
            # print(f"  Testing Fs={fs} Hz")
            for f_test, target_db, tol in self.a_targets:
                if f_test > fs / 2.2: continue  # Skip if near Nyquist
                
                measured = self._measure_gain(fs, f_test, 'A')
                err_msg = (f"Fs={fs}, Freq={f_test}Hz: Expected {target_db} +/- {tol}, "
                           f"Got {measured:.2f}")
                self.assertAlmostEqual(measured, target_db, delta=tol, msg=err_msg)

    def test_c_weighting_compliance(self):
        print("\n--- Testing C-Weighting Compliance ---")
        for fs in SUPPORTED_FS:
            # print(f"  Testing Fs={fs} Hz")
            for f_test, target_db, tol in self.c_targets:
                if f_test > fs / 2.2: continue
                
                measured = self._measure_gain(fs, f_test, 'C')
                err_msg = (f"Fs={fs}, Freq={f_test}Hz: Expected {target_db} +/- {tol}, "
                           f"Got {measured:.2f}")
                self.assertAlmostEqual(measured, target_db, delta=tol, msg=err_msg)
                
    def test_streaming_continuity(self):
        print("\n--- Testing Weighting Filter State Continuity ---")
        # Ensure that processing in chunks yields identical results to processing at once
        fs = 48000
        wf_continuous = WeightingFilter(fs, 'A')
        wf_chunked = WeightingFilter(fs, 'A')
        
        # Generate random noise signal
        full_signal = np.random.randn(fs) # 1 second of noise
        
        # Case A: Continuous
        out_cont = wf_continuous.process_chunk(full_signal)
        
        # Case B: Chunked (e.g., 100ms blocks)
        chunk_size = int(fs * 0.1)
        out_chunks = []
        for i in range(0, len(full_signal), chunk_size):
            chunk = full_signal[i:i+chunk_size]
            out_chunks.append(wf_chunked.process_chunk(chunk))
        out_chunked_stitched = np.concatenate(out_chunks)
        
        # Verify strict equality (floating point tolerance)
        max_diff = np.max(np.abs(out_cont - out_chunked_stitched))
        self.assertLess(max_diff, 1e-12, "Chunked processing did not match continuous processing!")


class TestOctaveFilterBank(unittest.TestCase):
    
    def test_octave_selectivity(self):
        print("\n--- Testing Octave Band Selectivity (1kHz Tone) ---")
        fs = 48000
        
        # The default order is now 24 (High Performance)
        bank = OctaveFilterBank(fs, resolution='octave', order=24)
        
        # Generate 1kHz Tone
        duration = 1.0
        t = np.linspace(0, duration, int(fs * duration), endpoint=False)
        signal = np.sqrt(2) * np.sin(2 * np.pi * 1000 * t)
        
        # Process
        output_matrix = bank.process_chunk(signal)
        
        # Calculate RMS for each band
        start_idx = int(len(signal) * 0.5)
        band_rms = np.sqrt(np.mean(output_matrix[start_idx:]**2, axis=0))
        band_db = 20 * np.log10(band_rms + 1e-15)
        
        # Check 1kHz Band
        idx_1k = np.argmin(np.abs(bank.frequencies - 1000))
        print(f"  1kHz Band Level: {band_db[idx_1k]:.2f} dB (Expected ~0 dB)")
        self.assertAlmostEqual(band_db[idx_1k], 0.0, delta=0.5)
        
        # Check adjacent band (500Hz)
        idx_500 = np.argmin(np.abs(bank.frequencies - 500))
        attenuation = band_db[idx_1k] - band_db[idx_500]
        print(f"  500Hz Band Level: {band_db[idx_500]:.2f} dB (Attenuation: {attenuation:.2f} dB)")
        self.assertGreater(attenuation, 60.0)

    def test_third_octave_selectivity(self):
        print("\n--- Testing 1/3 Octave Band Selectivity (1kHz Tone) ---")
        fs = 48000
        
        # Initialize 1/3 Octave Bank
        bank = OctaveFilterBank(fs, resolution='third', order=24)
        
        # Generate 1kHz Tone
        duration = 1.0
        t = np.linspace(0, duration, int(fs * duration), endpoint=False)
        signal = np.sqrt(2) * np.sin(2 * np.pi * 1000 * t)
        
        # Process
        output_matrix = bank.process_chunk(signal)
        
        # Calculate RMS for each band
        start_idx = int(len(signal) * 0.5)
        band_rms = np.sqrt(np.mean(output_matrix[start_idx:]**2, axis=0))
        band_db = 20 * np.log10(band_rms + 1e-15)
        
        # 1. Check 1kHz Band (Should be 0dB)
        idx_1k = np.argmin(np.abs(bank.frequencies - 1000))
        level_1k = band_db[idx_1k]
        print(f"  1kHz Band Level: {level_1k:.2f} dB (Expected ~0 dB)")
        self.assertAlmostEqual(level_1k, 0.0, delta=0.5, 
                               msg="1kHz tone did not register correctly in 1/3 octave band")
        
        # 2. Check Adjacent Lower Band (800 Hz)
        # 1/3 Octave spacing is 2^(1/3) approx 1.26
        # 1000 / 1.2599 = 793.7 Hz (ANSI standard label is 800)
        idx_800 = np.argmin(np.abs(bank.frequencies - 793.7))
        
        # 800 Hz band edge is approx 891 Hz. 1kHz tone is well outside.
        level_800 = band_db[idx_800]
        attenuation = level_1k - level_800
        
        print(f"  800Hz Band Level: {level_800:.2f} dB (Attenuation: {attenuation:.2f} dB)")
        
        # FIX: ANSI S1.11 Class 1 requires > 17.5 dB attenuation at this point.
        # This implementation achieves ~21.5 dB, which passes.
        self.assertGreater(attenuation, 18.0, 
                           "Leakage into 800Hz 1/3 octave band is too high (Must be > 17.5dB)")

        # 3. Check Distant Band (630 Hz - 2 bands away)
        # This should be deeply attenuated (>60dB)
        idx_630 = np.argmin(np.abs(bank.frequencies - 630))
        level_630 = band_db[idx_630]
        attenuation_dist = level_1k - level_630
        
        print(f"  630Hz Band Level: {level_630:.2f} dB (Attenuation: {attenuation_dist:.2f} dB)")
        self.assertGreater(attenuation_dist, 60.0, 
                           "Leakage into 630Hz band (2 bands away) should be negligible")

    def test_streaming_continuity(self):
        print("\n--- Testing Octave Bank State Continuity ---")
        fs = 48000
        bank_continuous = OctaveFilterBank(fs, 'octave')
        bank_chunked = OctaveFilterBank(fs, 'octave')
        
        signal = np.random.randn(fs) # 1 sec noise
        
        # Continuous
        res_cont = bank_continuous.process_chunk(signal)
        
        # Chunked
        chunk_size = 4800
        chunks = []
        for i in range(0, len(signal), chunk_size):
            c = signal[i:i+chunk_size]
            chunks.append(bank_chunked.process_chunk(c))
        res_stitched = np.vstack(chunks)
        
        max_diff = np.max(np.abs(res_cont - res_stitched))
        self.assertLess(max_diff, 1e-12, "Bank chunked processing mismatch")

if __name__ == '__main__':
    unittest.main()