# python2/vslm/calibration.py
import numpy as np
import soundfile as sf
from pathlib import Path

REF_PRESSURE = 20e-6 # 20 microPascals

def compute_selection_rms(filepath: Path, start_time: float, end_time: float) -> float:
    """
    Reads a specific time range from a WAV file and calculates the uncalibrated RMS amplitude.
    """
    if not filepath.exists():
        raise FileNotFoundError(f"File not found: {filepath}")

    with sf.SoundFile(str(filepath)) as f:
        sr = f.samplerate
        
        # Calculate frames
        start_frame = int(start_time * sr)
        end_frame = int(end_time * sr)
        duration_frames = end_frame - start_frame
        
        if duration_frames <= 0:
            raise ValueError("Invalid selection duration.")
            
        f.seek(start_frame)
        data = f.read(duration_frames, always_2d=True)
        
        # Mix to mono if necessary
        if data.shape[1] > 1:
            data = np.mean(data, axis=1)
        else:
            data = data.flatten()
            
        # Compute RMS
        # Add epsilon to prevent divide by zero issues later
        rms = np.sqrt(np.mean(data**2)) + 1e-15
        return float(rms)

def calculate_factor_from_ref(measured_rms: float, target_db: float) -> float:
    """
    Calculates the calibration factor required to map the measured RMS 
    to the target dB SPL.
    
    Formula: K = (P_ref * 10^(L_target/20)) / RMS_measured
    """
    target_pressure = REF_PRESSURE * (10 ** (target_db / 20.0))
    factor = target_pressure / measured_rms
    return factor