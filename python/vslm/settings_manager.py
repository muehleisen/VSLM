import yaml
from pathlib import Path
from datetime import datetime
from pydantic import BaseModel, Field, ConfigDict
from typing import Dict

# Import Constants to use actual Enums in the model
from .constants import Weighting, ResponseSpeed

DEFAULT_SETTINGS_FILE = Path.home() / ".vslm_settings.yaml"

class DoseStandard(BaseModel):
    """Model for individual dose calculation parameters."""
    exchange_rate: float
    criterion_level: float
    threshold_level: float
    shift_hours: float

class AppSettings(BaseModel):
    """
    Pydantic Model representing the persistent application state.
    """
    # File I/O
    last_directory: str = str(Path.home())
    
    # Acoustics (Global)
    ref_pressure: float = 20e-6
    
    # Calibration
    calibration_factor: float = Field(default=1.0, ge=0)
    
    # Analysis Configuration
    weighting: Weighting = Weighting.A
    speed: ResponseSpeed = ResponseSpeed.FAST
    
    analysis_mode_index: int = 0
    leq_interval_index: int = 1
    block_size_ms: float = Field(default=100.0, gt=0)
    
    # Advanced Analysis
    band_filter_order: int = 24
    
    # Dose Settings
    current_dose_standard: str = 'NIOSH'
    dose_standards: Dict[str, DoseStandard] = Field(default_factory=lambda: {
        'NIOSH': DoseStandard(exchange_rate=3.0, criterion_level=85.0, threshold_level=80.0, shift_hours=8.0),
        'OSHA':  DoseStandard(exchange_rate=5.0, criterion_level=90.0, threshold_level=80.0, shift_hours=8.0)
    })

    # Plot Scaling Settings
    plot_autoscale: bool = True
    plot_ymin: float = 0.0
    plot_ymax: float = 120.0

    # Allow extra fields in YAML without crashing
    # validate_assignment=True forces Pydantic to convert "A" -> Weighting.A when assigned!
    model_config = ConfigDict(extra='ignore', validate_assignment=True)

class SettingsManager:
    """Handles loading and saving AppSettings to a YAML file using Pydantic."""
    def __init__(self, default_path: Path = DEFAULT_SETTINGS_FILE):
        self.default_path = default_path

    def load(self, filepath: Path = None) -> AppSettings:
        path_to_load = filepath if filepath else self.default_path
        
        if not path_to_load.exists():
            return AppSettings()
        
        try:
            with open(path_to_load, 'r') as f:
                data = yaml.safe_load(f)
                if not data: 
                    return AppSettings()
                return AppSettings.model_validate(data)
        except Exception as e:
            print(f"Error loading settings: {e}")
            return AppSettings()

    def save(self, settings: AppSettings, filepath: Path = None):
        path_to_save = filepath if filepath else self.default_path
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        try:
            settings_dict = settings.model_dump(mode='json')
            
            with open(path_to_save, 'w') as f:
                f.write(f"# Saved by VSLM on {timestamp}\n")
                yaml.dump(settings_dict, f, default_flow_style=False)
        except Exception as e:
            print(f"Error saving settings: {e}")