import unittest
import numpy as np
import sys
import os

# Path Hack
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from vslm.leq import calculate_leq_analysis, LeqStats

class TestLeqAnalysis(unittest.TestCase):
    
    def create_blocks(self, duration_s, level_db, block_ms=100):
        """Helper to create a list of block results with constant level."""
        n_blocks = int(duration_s * 1000 / block_ms)
        blocks = []
        for i in range(n_blocks):
            blocks.append({
                'time': i * (block_ms/1000.0), 
                'leq': level_db
            })
        return blocks

    def test_constant_level(self):
        print("\n--- Testing LEQ: Constant 94dB Signal ---")
        # 10 seconds of 94 dB
        blocks = self.create_blocks(10.0, 94.0)
        
        stats = calculate_leq_analysis(blocks, stats_block_ms=100)
        
        # 1. Overall Leq should be exactly 94.0
        self.assertAlmostEqual(stats.overall, 94.0, delta=0.1)
        
        # 2. Max and Min should be 94.0
        self.assertAlmostEqual(stats.max, 94.0, delta=0.1)
        self.assertAlmostEqual(stats.min, 94.0, delta=0.1)
        
        # 3. All Percentiles should be 94.0
        self.assertAlmostEqual(stats.ln[10], 94.0, delta=0.1)
        self.assertAlmostEqual(stats.ln[90], 94.0, delta=0.1)

    def test_energy_averaging(self):
        print("\n--- Testing LEQ: Energy Averaging (80dB vs 100dB) ---")
        # 5 seconds of 80 dB, 5 seconds of 100 dB
        # Energy difference is factor of 100 (20 dB)
        b1 = self.create_blocks(5.0, 80.0)
        b2 = self.create_blocks(5.0, 100.0)
        blocks = b1 + b2 # Concatenate
        
        stats = calculate_leq_analysis(blocks)
        
        # Calculation:
        # P80 = ref * 10^(80/20) = ref * 10^4
        # P100 = ref * 10^(100/20) = ref * 10^5
        # P_sq_avg = (0.5 * P80^2 + 0.5 * P100^2)
        # P80^2 is negligible compared to P100^2 (1/100th size)
        # Result should be approx 100 dB - 3 dB (since active half time) = 97 dB
        
        expected_leq = 10 * np.log10(0.5 * (10**(80/10)) + 0.5 * (10**(100/10)))
        print(f"  Expected: {expected_leq:.2f} dB, Got: {stats.overall:.2f} dB")
        
        self.assertAlmostEqual(stats.overall, expected_leq, delta=0.1)
        self.assertEqual(stats.max, 100.0)
        self.assertEqual(stats.min, 80.0)

    def test_dose_niosh(self):
        print("\n--- Testing Dose: NIOSH (85dB Criterion, 3dB Exchange) ---")
        # 1 Hour of 85 dB
        # This is 1/8th of an 8-hour day.
        # Since level == criterion, accumulation rate is 100%.
        # Total dose should be 100% * (1hr / 8hr) = 12.5%
        
        blocks = self.create_blocks(3600.0, 85.0) # 1 hour
        
        stats = calculate_leq_analysis(blocks, dose_standard='NIOSH')
        
        print(f"  NIOSH Dose (1hr @ 85dB): {stats.dose['dose']:.2f}% (Expected 12.5%)")
        self.assertAlmostEqual(stats.dose['dose'], 12.5, delta=0.1)
        
        # TWA should be calculated based on the dose fraction extrapolated to 8 hours? 
        # Actually TWA formula: TWA = 10 log(D/100) + 85 (where D is percent)
        # 12.5% dose -> TWA = 10 log(0.125) + 85 = -9 + 85 = 76 dB 
        # WAIT: TWA is usually the level that would produce that dose over 8 hours.
        # If I exposed for 1 hour at 85, my 8-hour TWA is 76.
        # However, many meters display "Projected TWA".
        # Let's check the formula in leq.py:
        # TWA = 10 * log10(fraction) + criterion. 
        # log10(0.125) = -0.903. 10*-0.903 = -9.03. 85 - 9.03 = 75.97 dB.
        self.assertAlmostEqual(stats.dose['twa'], 75.97, delta=0.1)

    def test_dose_osha(self):
        print("\n--- Testing Dose: OSHA (90dB Criterion, 5dB Exchange) ---")
        # 1 Hour of 90 dB
        # This is 1/8th of a day at criterion level.
        # Dose = 12.5%
        
        blocks = self.create_blocks(3600.0, 90.0)
        stats = calculate_leq_analysis(blocks, dose_standard='OSHA')
        
        print(f"  OSHA Dose (1hr @ 90dB): {stats.dose['dose']:.2f}% (Expected 12.5%)")
        self.assertAlmostEqual(stats.dose['dose'], 12.5, delta=0.1)
        
        # Test Threshold: OSHA threshold is 80dB (Hearing Conservation) or 90dB (PEL)
        # VSLM uses 80dB threshold for both by default logic.
        # Let's test below threshold (75dB). Dose should be 0.
        blocks_quiet = self.create_blocks(60.0, 75.0)
        stats_quiet = calculate_leq_analysis(blocks_quiet, dose_standard='OSHA')
        self.assertEqual(stats_quiet.dose['dose'], 0.0)

    def test_time_aggregation(self):
        print("\n--- Testing Time History Aggregation ---")
        # 10 seconds of data in 100ms blocks (100 blocks total)
        blocks = self.create_blocks(10.0, 90.0, block_ms=100)
        
        # Request 1-second aggregation
        stats = calculate_leq_analysis(blocks, stats_block_ms=100, integration_time_s=1.0)
        
        hist_len = len(stats.history['leq'])
        print(f"  Input Blocks: {len(blocks)}, Output Intervals: {hist_len}")
        
        # Should be exactly 10 intervals
        self.assertEqual(hist_len, 10)
        
        # Each interval should be 90 dB
        self.assertAlmostEqual(stats.history['leq'][0], 90.0, delta=0.1)
        
        # Check time stamps (0, 1, 2, ... 9)
        self.assertEqual(stats.history['time'][0], 0.0)
        self.assertEqual(stats.history['time'][1], 1.0)
        self.assertEqual(stats.history['time'][-1], 9.0)

if __name__ == '__main__':
    unittest.main()