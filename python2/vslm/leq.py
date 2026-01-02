# python2/vslm/leq.py
import numpy as np

class LeqStats:
    def __init__(self, overall_leq, l_max, l_min, percentiles, dose_result, time_history, stats_block_size_ms):
        self.overall = overall_leq
        self.max = l_max
        self.min = l_min
        self.ln = percentiles # Dict {10: val, 90: val, ...}
        self.dose = dose_result # Dict {twa: val, dose_pct: val}
        self.history = time_history # Dict {time: [], leq: []}
        self.stats_block_size_ms = stats_block_size_ms

def calculate_leq_analysis(block_results, 
                           stats_block_ms=100, 
                           integration_time_s=1.0, 
                           dose_standard='NIOSH'):
    """
    Performs full LEQ analysis on a sequence of short-term (e.g. 100ms) blocks.
    
    Args:
        block_results (list): List of dicts {'time': t, 'leq': db} from StreamProcessor.
        stats_block_ms (float): The size of blocks used for statistics (default 100ms).
        integration_time_s (float): The aggregation interval for the time history plot.
        dose_standard (str): 'NIOSH' or 'OSHA'.
        
    Returns:
        LeqStats object.
    """
    # 1. Extract raw 100ms data
    # We convert back to pressure squared for accurate averaging
    raw_db = np.array([b['leq'] for b in block_results])
    raw_pressure_sq = (10**(raw_db/10.0)) * (20e-6**2)
    
    # --- Statistics (based on 100ms blocks) ---
    
    # Overall LEQ (Energy Average of whole file)
    overall_msq = np.mean(raw_pressure_sq)
    overall_leq = 10 * np.log10(overall_msq / (20e-6**2) + 1e-30)
    
    # Min / Max
    l_max = np.max(raw_db)
    l_min = np.min(raw_db)
    
    # Percentiles (Ln)
    # L10 is the level exceeded 10% of the time (which is the 90th percentile)
    percentiles = {}
    for n in [10, 20, 30, 40, 50, 60, 70, 80, 90]:
        # numpy percentile q is 0-100 (lower to higher)
        # L10 = 90th percentile
        p_val = np.percentile(raw_db, 100 - n)
        percentiles[n] = p_val
        
    # Noise Dose Calculation
    # Defaults for NIOSH
    exchange_rate = 3
    criterion_level = 85
    threshold_level = 80
    
    if dose_standard == 'OSHA':
        exchange_rate = 5
        criterion_level = 90
        threshold_level = 80 # OSHA Hearing Conservation (80), Permissible (90)
        
    # Formula from VSLM Matlab:
    # q = Dex / log10(2)
    # D = sum(dt * 10^((L - DLc)/q)) / Tn
    # Tn = 480 * 60 (8 hours)
    
    q = exchange_rate / np.log10(2)
    Tn = 480 * 60.0
    dt = stats_block_ms / 1000.0
    
    # Filter for Threshold
    # Only levels > Threshold contribute to dose
    mask = raw_db > threshold_level
    dose_db = raw_db[mask]
    
    if len(dose_db) > 0:
        term = (dose_db - criterion_level) / q
        accumulated = np.sum(dt * (10**term))
        dose_fraction = accumulated / Tn
        
        # TWA = 10 * log10(D) + DLc
        if dose_fraction > 0:
            twa = 10 * np.log10(dose_fraction) + criterion_level
        else:
            twa = 0
    else:
        dose_fraction = 0.0
        twa = 0.0
        
    dose_result = {'dose': dose_fraction * 100.0, 'twa': twa, 'standard': dose_standard}

    # --- Time History Aggregation ---
    # Combine 100ms blocks into 'integration_time_s' chunks (e.g. 1s or 1min)
    
    blocks_per_interval = int(integration_time_s / (stats_block_ms / 1000.0))
    if blocks_per_interval < 1: 
        blocks_per_interval = 1
        
    num_intervals = len(raw_pressure_sq) // blocks_per_interval
    
    agg_time = []
    agg_leq = []
    
    for i in range(num_intervals):
        start_idx = i * blocks_per_interval
        end_idx = start_idx + blocks_per_interval
        
        chunk = raw_pressure_sq[start_idx:end_idx]
        mean_p = np.mean(chunk)
        db = 10 * np.log10(mean_p / (20e-6**2) + 1e-30)
        
        # Time stamp is the start of the interval
        t = i * integration_time_s
        
        agg_time.append(t)
        agg_leq.append(db)
        
    return LeqStats(overall_leq, l_max, l_min, percentiles, dose_result, 
                    {'time': agg_time, 'leq': agg_leq}, stats_block_ms)