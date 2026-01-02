# python2/vslm/leq.py
import numpy as np
from dataclasses import dataclass

@dataclass
class LeqStats:
    overall: float
    max: float
    min: float
    ln: dict[int, float]
    dose: dict[str, float | str]
    history: dict[str, list[float]]
    stats_block_size_ms: float

def calculate_leq_analysis(block_results: list[dict], 
                           stats_block_ms: float, 
                           integration_time_s: float, 
                           dose_params: dict,
                           ref_pressure: float = 20e-6) -> LeqStats:
    """
    Performs full LEQ analysis using configurable physics/dose parameters.
    """
    # 1. Extract raw data (dB)
    raw_db = np.array([b['leq'] for b in block_results])
    
    # Convert to pressure squared using Configurable Reference
    raw_pressure_sq = (10**(raw_db/10.0)) * (ref_pressure**2)
    
    # --- Statistics ---
    overall_msq = np.mean(raw_pressure_sq)
    overall_leq = 10 * np.log10(overall_msq / (ref_pressure**2) + 1e-30)
    
    l_max = np.max(raw_db)
    l_min = np.min(raw_db)
    
    percentiles = {
        n: np.percentile(raw_db, 100 - n) 
        for n in [10, 20, 30, 40, 50, 60, 70, 80, 90]
    }
        
    # --- Configurable Noise Dose ---
    er = dose_params.get('exchange_rate', 3.0)
    cl = dose_params.get('criterion_level', 85.0)
    tl = dose_params.get('threshold_level', 80.0)
    hours = dose_params.get('shift_hours', 8.0)
    
    q = er / np.log10(2)
    Tn = hours * 60.0 # Shift in minutes
    dt = stats_block_ms / 1000.0
    
    mask = raw_db > tl
    dose_db = raw_db[mask]
    
    if dose_db.size > 0:
        term = (dose_db - cl) / q
        accumulated = np.sum(dt * (10**term))
        dose_fraction = accumulated / Tn
        twa = (10 * np.log10(dose_fraction) + cl) if dose_fraction > 0 else 0.0
    else:
        dose_fraction = 0.0
        twa = 0.0
        
    dose_result = {
        'dose': dose_fraction * 100.0, 
        'twa': twa, 
        'standard': 'Custom' # Or passed in label
    }

    # --- Time History Aggregation ---
    blocks_per_interval = int(integration_time_s / (stats_block_ms / 1000.0))
    if blocks_per_interval < 1: 
        blocks_per_interval = 1
        
    n_total = len(raw_pressure_sq)
    n_intervals = n_total // blocks_per_interval
    
    if n_intervals > 0:
        trimmed_sq = raw_pressure_sq[:n_intervals*blocks_per_interval]
        reshaped = trimmed_sq.reshape(n_intervals, blocks_per_interval)
        means = np.mean(reshaped, axis=1)
        agg_leq = 10 * np.log10(means / (ref_pressure**2) + 1e-30)
        agg_time = np.arange(n_intervals) * integration_time_s
    else:
        agg_leq = np.array([])
        agg_time = np.array([])
        
    return LeqStats(
        overall=overall_leq, 
        max=l_max, 
        min=l_min, 
        ln=percentiles, 
        dose=dose_result, 
        history={'time': agg_time.tolist(), 'leq': agg_leq.tolist()}, 
        stats_block_size_ms=stats_block_ms
    )