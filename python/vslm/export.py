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
    def export_leq(filepath: Path, results: list, block_size_ms: float, 
                   interval_txt: str, weighting: str,
                   dose_params: dict, ref_pressure: float): # Updated Signature
        """Exports integrated LEQ history based on the selected interval."""
        
        # Map text to seconds
        match interval_txt:
            case "100 ms": interval = 0.1
            case "1 sec": interval = 1.0
            case "10 sec": interval = 10.0
            case "1 min": interval = 60.0
            case "15 min": interval = 900.0
            case "1 hour": interval = 3600.0
            case _: interval = 1.0

        # Calculate stats (Now passing the required args)
        stats = leq.calculate_leq_analysis(
            results, block_size_ms, interval, dose_params, ref_pressure
        )
        
        with open(filepath, 'w', newline='') as f:
            writer = csv.writer(f)
            # Add a metadata header
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
        
        # Check if bands exist
        if 'bands' not in results[0]:
            return

        freqs = results[0]['band_freqs']
        
        # Energy Average
        energy_sums = np.zeros(len(freqs))
        count = len(results)
        for r in results:
            # Convert dB back to Pressure^2 using Ref Pressure
            pressures = (10**(r['bands']/10.0)) * (ref_pressure**2)
            energy_sums += pressures
        
        if count > 0:
            mean_pressure_sq = energy_sums / count
            # Convert back to dB using Ref Pressure
            mean_db = 10 * np.log10(mean_pressure_sq / (ref_pressure**2) + 1e-30)
        else:
            mean_db = np.zeros(len(freqs))
            
        with open(filepath, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["# VSLM Export", f"Spectrum ({weighting})"])
            writer.writerow(["Frequency (Hz)", f"Average Level [dB]"])
            for freq, level in zip(freqs, mean_db):
                writer.writerow([f"{freq:.1f}", f"{level:.2f}"])