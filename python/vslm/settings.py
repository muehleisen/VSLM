import yaml
from dataclasses import dataclass, field, asdict
from pathlib import Path
from datetime import datetime

DEFAULT_SETTINGS_FILE = Path.home() / ".vslm_settings.yaml"

@dataclass
class AppSettings:
    """
    Data Object representing the persistent application state.
    """
    # File I/O
    last_directory: str = str(Path.home())
    
    # Acoustics (Global)
    ref_pressure: float = 20e-6
    
    # Calibration
    calibration_factor: float = 1.0
    
    # Analysis Configuration
    weighting: str = "A"           
    speed: str = "Fast"            
    analysis_mode_index: int = 0   
    leq_interval_index: int = 1    
    block_size_ms: float = 100.0   
    
    # Advanced Analysis
    band_filter_order: int = 24    
    
    # Dose Settings
    current_dose_standard: str = 'NIOSH'
    dose_standards: dict = field(default_factory=lambda: {
        'NIOSH': {'exchange_rate': 3.0, 'criterion_level': 85.0, 'threshold_level': 80.0, 'shift_hours': 8.0},
        'OSHA':  {'exchange_rate': 5.0, 'criterion_level': 90.0, 'threshold_level': 80.0, 'shift_hours': 8.0}
    })

    # NEW: Plot Scaling Settings
    plot_autoscale: bool = True
    plot_ymin: float = 0.0
    plot_ymax: float = 120.0

class SettingsManager:
    """Handles loading and saving AppSettings to a YAML file."""
    def __init__(self, default_path: Path = DEFAULT_SETTINGS_FILE):
        self.default_path = default_path

    def load(self, filepath: Path = None) -> AppSettings:
        path_to_load = filepath if filepath else self.default_path
        
        if not path_to_load.exists():
            return AppSettings()
        
        try:
            with open(path_to_load, 'r') as f:
                data = yaml.safe_load(f)
                if not data: return AppSettings()
                valid_keys = AppSettings.__dataclass_fields__.keys()
                filtered = {k: v for k, v in data.items() if k in valid_keys}
                return AppSettings(**filtered)
        except Exception as e:
            print(f"Error loading settings: {e}")
            return AppSettings()

    def save(self, settings: AppSettings, filepath: Path = None):
        path_to_save = filepath if filepath else self.default_path
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        try:
            with open(path_to_save, 'w') as f:
                f.write(f"# Saved by vlsm.py on {timestamp}\n")
                yaml.dump(asdict(settings), f, default_flow_style=False)
        except Exception as e:
            print(f"Error saving settings: {e}")