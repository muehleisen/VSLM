import numpy as np
from dataclasses import dataclass
from typing import TYPE_CHECKING, Dict, List, Union

if TYPE_CHECKING:
    from .settings_manager import DoseStandard

@dataclass
class LeqStats:
    overall: float
    max: float
    min: float
    ln: Dict[int, float]
    dose: Dict[str, float]
    history: Dict[str, List[float]]
    stats_block_size_ms: float

def calculate_leq_analysis(block_results: list, 
                           stats_block_ms: float, 
                           integration_time_s: float, 
                           dose_params: "DoseStandard", 
                           ref_pressure: float = 20e-6) -> LeqStats:
    """
    Performs full LEQ analysis using configurable physics/dose parameters.
    """
    if not block_results:
        return LeqStats(
            overall=-100.0, max=-100.0, min=-100.0, 
            ln={n: -100.0 for n in [10, 20, 30, 40, 50, 60, 70, 80, 90]},
            dose={'dose': 0.0, 'twa': 0.0},
            history={'time': [], 'leq': []},
            stats_block_size_ms=stats_block_ms
        )

    # Extract LEQ values
    raw_db = np.array([b.get('leq', -100.0) for b in block_results])
    
    # Calculate Energy Average (Overall LEQ)
    raw_pressure_sq = (10**(raw_db/10.0)) * (ref_pressure**2)
    overall_msq = np.mean(raw_pressure_sq)
    overall_leq = 10 * np.log10(overall_msq / (ref_pressure**2) + 1e-30)
    
    # Min / Max
    l_max = np.max(raw_db)
    l_min = np.min(raw_db)
    
    # Percentiles (Ln)
    percentiles = {
        n: np.percentile(raw_db, 100 - n) 
        for n in [10, 20, 30, 40, 50, 60, 70, 80, 90]
    }
        
    # --- FIX: Use Attribute Access for Pydantic Object ---
    # dose_params is an object, not a dictionary.
    er = dose_params.exchange_rate
    cl = dose_params.criterion_level
    tl = dose_params.threshold_level
    hours = dose_params.shift_hours
    
    # Dose Calculation Logic
    q = er / np.log10(2)
    Tn = hours * 60.0 # Shift in minutes
    dt = stats_block_ms / 1000.0 # Block duration in seconds
    
    # Filter for Threshold
    mask = raw_db > tl
    dose_db = raw_db[mask]
    
    if dose_db.size > 0:
        term = (dose_db - cl) / q
        accumulated = np.sum(dt * (10**term))
        # Dose is a percentage relative to the allowed exposure
        # Note: Standard formula often uses Time in Hours for accumulation or T_criterion
        # Here we use the standard definition: D = 100 * (C/T)
        # We'll use the accumulated time-weighted values.
        
        # Simplified NIOSH/OSHA Dose Formula integration
        dose_fraction = accumulated / (Tn * 60.0) # Tn is minutes, convert to seconds
        # Actually, usually accumulated is sum(C/T), let's stick to the basic energy-like accumulation:
        # D = 100/Tc * sum(dt * 2^((L-CL)/ER))
        # 2^x = 10^(x * log10(2)) -> 10^((L-CL)/q) where q = ER/log10(2)
        # So 'term' above is correct.
        # The denominator represents the Criterion Duration (usually 8 hours)
        
        target_seconds = Tn * 60.0
        dose_percentage = 100.0 * (accumulated / target_seconds)
        
        # TWA Calculation
        # TWA = CL + q * log10(Dose/100 * (TargetDuration / MeasurementDuration)?) 
        # Or TWA = CL + q * log10(Dose/100) for an 8-hour projection?
        # Let's use the standard projection to 8 hours:
        if dose_percentage > 0:
            twa = cl + q * np.log10(dose_percentage / 100.0)
        else:
            twa = 0.0
    else:
        dose_percentage = 0.0
        twa = 0.0
        
    dose_result = {
        'dose': dose_percentage, 
        'twa': twa
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