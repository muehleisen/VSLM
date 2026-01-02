# vslm/settings.py
from dataclasses import dataclass

@dataclass
class VSLMSettings:
    """
    Central storage for all application settings and persistent UI state.
    """

    # --- Global Analysis State (UI Buttons) ---
    frequency_weighting: str = 'A'       # Options: 'A', 'C', 'Z'
    meter_speed: str = 'Slow'            # Options: 'Slow', 'Fast', 'Impulse'
    analysis_mode: str = 'lp'            # Options: 'lp', 'leq', 'octave', 'third', 'psd', 'spec'

    # --- Lp Plot Settings ---
    plot_time_spacing: float = 1.0       # Time axis tick interval (s)
    lpplot_y_min: int = 30               # Y-Axis Min (dB)
    lpplot_y_max: int = 120              # Y-Axis Max (dB)
    lpplot_autoscale: bool = True        # Auto-scale Y-Axis

    # --- Leq & Dose Settings ---
    leq_integration_time: float = 1.0    # Integration interval (s)
    leq_percentile: int = 90             # Ln percentile (e.g., L90)
    leq_y_min: int = 30
    leq_y_max: int = 120
    leq_autoscale: bool = True
    
    # Noise Dose Criteria
    dose_exchange_rate: int = 3          # Exchange Rate (dB)
    dose_threshold: int = 80             # Threshold Level (dB)
    dose_criterion: int = 90             # Criterion Level (dB)

    # --- Band Analysis Settings ---
    band_resolution_index: int = 0       # ComboBox Index (0=1/1 Octave, 1=1/3, etc.)
    band_method_index: int = 0           # ComboBox Index (0=ANSI, 1=FFT)

    # --- PSD (Power Spectral Density) Settings ---
    psd_fft_size: int = 4096             # N Points
    psd_overlap_percent: int = 50        # Overlap %
    psd_window: str = "Hann"             # Window function name
    psd_y_min: int = 30                  # Y-Axis Min (dB)
    psd_y_max: int = 120                 # Y-Axis Max (dB)
    psd_autoscale: bool = True           # Auto-scale Y-Axis

    # --- Spectrogram Settings ---
    spec_fft_size: int = 4096            # Frequency resolution
    spec_slice_length: float = 0.1       # Time resolution (s)
    spec_overlap_percent: int = 50       # Overlap % (Alternative to Slice Length)
    spec_use_overlap: bool = False       # Flag: True=Use Overlap, False=Use Slice Length
    spec_y_min: int = 30                 # Colorbar Min (dB)
    spec_y_max: int = 120                # Colorbar Max (dB)
    spec_autoscale: bool = True          # Auto-scale Colorbar
    spec_colormap: str = "inferno"       # Matplotlib colormap name