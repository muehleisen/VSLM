import unittest
import numpy as np
import soundfile as sf
import os
import sys

# Path Hack
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from vslm.analysis_engine import StreamProcessor

class TestAnalysisEngine(unittest.TestCase):
    
    def setUp(self):
        # Create a temporary test WAV file (1 second, 1kHz tone, 48kHz)
        self.test_file = 'temp_engine_test.wav'
        self.fs = 48000
        self.duration = 1.0
        
        t = np.linspace(0, self.duration, int(self.fs * self.duration), endpoint=False)
        
        # FIX: Generate signal with RMS = 0.5 (Peak ~= 0.707)
        # This fits safely within [-1.0, 1.0] to avoid clipping in the WAV file.
        # We will use cal_factor=2.0 later to simulate a 1.0 Pascal signal.
        self.signal_rms = 0.5
        self.signal = np.sqrt(2) * self.signal_rms * np.sin(2 * np.pi * 1000 * t)
        
        sf.write(self.test_file, self.signal, self.fs)

    def tearDown(self):
        # Consume generator fully or use try/except to ensure file release
        if os.path.exists(self.test_file):
            try:
                os.remove(self.test_file)
            except PermissionError:
                pass # File might still be locked if test crashed hard

    def test_run_analysis_broadband(self):
        print("\n--- Testing Analysis Engine: Broadband ---")
        
        # FIX: Use cal_factor=2.0 to restore the 0.5 RMS signal to 1.0 Pascal effective
        processor = StreamProcessor(self.test_file, cal_factor=2.0)
        
        # Run Analysis (100ms blocks, A-weighted)
        results = list(processor.run_analysis(block_size_ms=100, weighting='A'))
        
        self.assertEqual(len(results), 10, "Should yield exactly 10 blocks for 1s file")
        
        # Check STEADY STATE block (index 5)
        block_idx = 5
        
        # Expected: 1.0 Pascal effective -> 93.98 dB
        expected_leq = 20 * np.log10(1.0 / 20e-6)
        
        print(f"  Block {block_idx} Leq: {results[block_idx]['leq']:.2f} dB (Expected ~{expected_leq:.2f})")
        self.assertAlmostEqual(results[block_idx]['leq'], expected_leq, delta=0.1)

    def test_run_analysis_bands(self):
        print("\n--- Testing Analysis Engine: Octave Bands ---")
        
        # FIX: Use cal_factor=2.0 here too
        processor = StreamProcessor(self.test_file, cal_factor=2.0)
        
        # Run with Octave Bands (Z-weighted)
        results = list(processor.run_analysis(block_size_ms=100, weighting='Z', do_band_analysis=True))
        
        # Check Steady State (Block 5)
        result = results[5]
        
        self.assertIn('bands', result)
        self.assertIn('band_freqs', result)
        
        # Find 1kHz band
        freqs = result['band_freqs']
        levels = result['bands']
        idx_1k = np.argmin(np.abs(freqs - 1000))
        
        expected_leq = 20 * np.log10(1.0 / 20e-6)
        print(f"  1kHz Band Level: {levels[idx_1k]:.2f} dB")
        
        # Should contain nearly all energy (~94dB)
        self.assertAlmostEqual(levels[idx_1k], expected_leq, delta=0.5)
        
        # Check adjacent band (500Hz) for leakage
        idx_500 = np.argmin(np.abs(freqs - 500))
        print(f"  500Hz Band Level: {levels[idx_500]:.2f} dB")
        self.assertLess(levels[idx_500], 20.0, "Leakage into 500Hz band should be negligible")

if __name__ == '__main__':
    unittest.main()