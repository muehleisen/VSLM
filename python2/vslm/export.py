# python2/vslm/export.py
import csv
from pathlib import Path
import numpy as np
from . import leq

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
                # Time, Lp
                writer.writerow([f"{r['time']:.3f}", f"{r['lp']:.2f}"])

    @staticmethod
    def export_leq(filepath: Path, results: list, block_size_ms: float, interval_txt: str, weighting: str):
        """Exports integrated LEQ history based on the selected interval."""
        # Map text to seconds (matching logic in GUI)
        match interval_txt:
            case "100 ms": interval = 0.1
            case "1 sec": interval = 1.0
            case "10 sec": interval = 10.0
            case "1 min": interval = 60.0
            case "15 min": interval = 900.0
            case "1 hour": interval = 3600.0
            case _: interval = 1.0

        # Calculate stats to get aggregated history
        stats = leq.calculate_leq_analysis(results, block_size_ms, interval)
        
        with open(filepath, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["Start Time (s)", f"LEQ ({weighting}, {interval_txt} interval) [dB]"])
            
            times = stats.history['time']
            levels = stats.history['leq']
            
            for t, l in zip(times, levels):
                writer.writerow([f"{t:.2f}", f"{l:.2f}"])

    @staticmethod
    def export_spectrum(filepath: Path, results: list, weighting: str):
        """Exports the time-averaged spectrum."""
        if not results: return
        
        # Check if bands exist
        if 'bands' not in results[0]:
            return

        freqs = results[0]['band_freqs']
        
        # Energy Average
        energy_sums = np.zeros(len(freqs))
        count = len(results)
        for r in results:
            # Convert dB back to Pressure^2
            pressures = (10**(r['bands']/10.0)) * (20e-6**2)
            energy_sums += pressures
        
        if count > 0:
            mean_pressure_sq = energy_sums / count
            mean_db = 10 * np.log10(mean_pressure_sq / (20e-6**2) + 1e-30)
        else:
            mean_db = np.zeros(len(freqs))
            
        with open(filepath, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["Frequency (Hz)", f"Average Level ({weighting}) [dB]"])
            for freq, level in zip(freqs, mean_db):
                writer.writerow([f"{freq:.1f}", f"{level:.2f}"])