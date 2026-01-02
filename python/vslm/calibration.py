# python2/vslm/calibration.py
import numpy as np
import soundfile as sf
from pathlib import Path

# Remove hardcoded constant
# REF_PRESSURE = 20e-6 

def compute_selection_rms(filepath: Path, start_time: float, end_time: float) -> float:
    # ... (Same as before, no hardcoded constants here) ...
    if not filepath.exists(): raise FileNotFoundError(f"File not found: {filepath}")
    with sf.SoundFile(str(filepath)) as f:
        sr = f.samplerate
        start_frame = int(start_time * sr)
        end_frame = int(end_time * sr)
        duration_frames = end_frame - start_frame
        if duration_frames <= 0: raise ValueError("Invalid selection duration.")
        f.seek(start_frame)
        data = f.read(duration_frames, always_2d=True)
        if data.shape[1] > 1: data = np.mean(data, axis=1)
        else: data = data.flatten()
        return float(np.sqrt(np.mean(data**2)) + 1e-15)

def calculate_factor_from_ref(measured_rms: float, target_db: float, ref_pressure: float = 20e-6) -> float:
    """
    Calculates calibration factor using configurable reference pressure.
    """
    target_pressure = ref_pressure * (10 ** (target_db / 20.0))
    factor = target_pressure / measured_rms
    return factor