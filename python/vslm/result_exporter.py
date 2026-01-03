import csv
from pathlib import Path
import numpy as np
from . import leq_calculator
from .constants import LEQ_INTERVAL_MAP # New Import

class ResultsExporter:
    """
    Handles exporting analysis results to CSV files.
    """
    
    @staticmethod
    def export_lp(filepath: Path, results: list, weighting: str, speed: str):
        """Exports raw Time vs Lp history."""
        with open(filepath, 'w', newline='') as f:
            writer = csv.writer(f)
            header = ["Time (s)", f"Lp ({weighting}, {speed}) [dB]"]
            writer.writerow(header)
            
            for r in results:
                writer.writerow([f"{r['time']:.3f}", f"{r['lp']:.2f}"])

    @staticmethod
    def export_leq(filepath: Path, results: list, block_size_ms: float, 
                   interval_key, # Expects Enum or Key 
                   weighting: str,
                   dose_params: dict, ref_pressure: float): 
        """Exports integrated LEQ history based on the selected interval."""
        
        # --- REFACTOR: Map Lookup ---
        if interval_key in LEQ_INTERVAL_MAP:
            interval_txt, interval_sec = LEQ_INTERVAL_MAP[interval_key]
        else:
            interval_txt, interval_sec = "1 sec", 1.0

        stats = leq.calculate_leq_analysis(
            results, block_size_ms, interval_sec, dose_params, ref_pressure
        )
        
        with open(filepath, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["# VSLM Export", f"Weighting: {weighting}", f"Interval: {interval_txt}"])
            writer.writerow(["Start Time (s)", f"LEQ [dB]"])
            
            times = stats.history['time']
            levels = stats.history['leq']
            
            for t, l in zip(times, levels):
                writer.writerow([f"{t:.2f}", f"{l:.2f}"])

    @staticmethod
    def export_spectrum(filepath: Path, results: list, weighting: str, ref_pressure: float):
        """Exports the time-averaged spectrum."""
        if not results: return
        
        if 'bands' not in results[0]:
            return

        freqs = results[0]['band_freqs']
        
        energy_sums = np.zeros(len(freqs))
        count = len(results)
        for r in results:
            pressures = (10**(r['bands']/10.0)) * (ref_pressure**2)
            energy_sums += pressures
        
        if count > 0:
            mean_pressure_sq = energy_sums / count
            mean_db = 10 * np.log10(mean_pressure_sq / (ref_pressure**2) + 1e-30)
        else:
            mean_db = np.zeros(len(freqs))
            
        with open(filepath, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["# VSLM Export", f"Spectrum ({weighting})"])
            writer.writerow(["Frequency (Hz)", f"Average Level [dB]"])
            for freq, level in zip(freqs, mean_db):
                writer.writerow([f"{freq:.1f}", f"{level:.2f}"])