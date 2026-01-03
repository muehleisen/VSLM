from enum import StrEnum

class AnalysisMode(StrEnum):
    """Modes used by the GUI to determine analysis and plotting behavior."""
    LP = "lp"
    LEQ = "leq"
    OCTAVE = "octave"
    THIRD_OCTAVE = "third_octave"

class Weighting(StrEnum):
    """Frequency weighting standards."""
    A = 'A'
    C = 'C'
    Z = 'Z'

class ResponseSpeed(StrEnum):
    """Time weighting detector speeds."""
    SLOW = 'Slow'
    FAST = 'Fast'
    IMPULSE = 'Impulse'

class BandResolution(StrEnum):
    """Resolution for frequency analysis."""
    OCTAVE = 'octave'
    THIRD_OCTAVE = 'third'