# python2/vslm/leq.py
import numpy as np
from dataclasses import dataclass

@dataclass
class LeqStats:
    """Data object holding LEQ statistical analysis results."""
    overall: float
    max: float
    min: float
    ln: dict[int, float]              # e.g. {10: 85.4, 90: 45.2}
    dose: dict[str, float | str]      # e.g. {'twa': 85.0, 'standard': 'NIOSH'}
    history: dict[str, list[float]]   # e.g. {'time': [...], 'leq': [...]}
    stats_block_size_ms: float

def calculate_leq_analysis(block_results: list[dict], 
                           stats_block_ms: float = 100.0, 
                           integration_time_s: float = 1.0, 
                           dose_standard: str = 'NIOSH') -> LeqStats:
    """
    Performs full LEQ analysis on a sequence of short-term blocks.
    """
    # 1. Extract raw data
    raw_db = np.array([b['leq'] for b in block_results])
    # Convert back to pressure squared for accurate averaging
    raw_pressure_sq = (10**(raw_db/10.0)) * (20e-6**2)
    
    # --- Statistics ---
    
    # Overall LEQ
    overall_msq = np.mean(raw_pressure_sq)
    overall_leq = 10 * np.log10(overall_msq / (20e-6**2) + 1e-30)
    
    l_max = np.max(raw_db)
    l_min = np.min(raw_db)
    
    # Percentiles (Ln) using dictionary comprehension
    # L10 is the 90th percentile value
    percentiles = {
        n: np.percentile(raw_db, 100 - n) 
        for n in [10, 20, 30, 40, 50, 60, 70, 80, 90]
    }
        
    # Noise Dose Calculation
    match dose_standard:
        case 'OSHA':
            exchange_rate = 5
            criterion_level = 90
            threshold_level = 80
        case _: # NIOSH (Default)
            exchange_rate = 3
            criterion_level = 85
            threshold_level = 80
        
    q = exchange_rate / np.log10(2)
    Tn = 480 * 60.0  # 8 hours in minutes
    dt = stats_block_ms / 1000.0
    
    # Vectorized Dose Calculation
    mask = raw_db > threshold_level
    dose_db = raw_db[mask]
    
    if dose_db.size > 0:
        term = (dose_db - criterion_level) / q
        accumulated = np.sum(dt * (10**term))
        dose_fraction = accumulated / Tn
        
        twa = (10 * np.log10(dose_fraction) + criterion_level) if dose_fraction > 0 else 0.0
    else:
        dose_fraction = 0.0
        twa = 0.0
        
    dose_result = {
        'dose': dose_fraction * 100.0, 
        'twa': twa, 
        'standard': dose_standard
    }

    # --- Time History Aggregation ---
    blocks_per_interval = int(integration_time_s / (stats_block_ms / 1000.0))
    if blocks_per_interval < 1: 
        blocks_per_interval = 1
        
    # Reshape for fast averaging (truncating remainder)
    n_total = len(raw_pressure_sq)
    n_intervals = n_total // blocks_per_interval
    
    if n_intervals > 0:
        # Truncate
        trimmed_sq = raw_pressure_sq[:n_intervals*blocks_per_interval]
        # Reshape to (n_intervals, blocks_per_interval)
        reshaped = trimmed_sq.reshape(n_intervals, blocks_per_interval)
        # Average across columns (axis 1)
        means = np.mean(reshaped, axis=1)
        # Convert to dB
        agg_leq = 10 * np.log10(means / (20e-6**2) + 1e-30)
        # Time axis
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