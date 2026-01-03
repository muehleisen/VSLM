from enum import StrEnum

class AnalysisMode(StrEnum):
    LP = "lp"
    LEQ = "leq"
    OCTAVE = "octave"
    THIRD_OCTAVE = "third_octave"

class Weighting(StrEnum):
    A = 'A'
    C = 'C'
    Z = 'Z'

class ResponseSpeed(StrEnum):
    SLOW = 'Slow'
    FAST = 'Fast'
    IMPULSE = 'Impulse'

class BandResolution(StrEnum):
    OCTAVE = 'octave'
    THIRD_OCTAVE = 'third'

class LeqInterval(StrEnum):
    MS_100 = "ms_100"
    SEC_1 = "sec_1"
    SEC_10 = "sec_10"
    MIN_1 = "min_1"
    MIN_15 = "min_15"
    HOUR_1 = "hour_1"

# Centralized Definition: Key -> (Display Text, Seconds)
LEQ_INTERVAL_MAP = {
    LeqInterval.MS_100: ("100 ms", 0.1),
    LeqInterval.SEC_1:  ("1 sec",  1.0),
    LeqInterval.SEC_10: ("10 sec", 10.0),
    LeqInterval.MIN_1:  ("1 min",  60.0),
    LeqInterval.MIN_15: ("15 min", 900.0),
    LeqInterval.HOUR_1: ("1 hour", 3600.0),
}

class DoseKeys(StrEnum):
    EXCHANGE_RATE = "exchange_rate"
    CRITERION_LEVEL = "criterion_level"
    THRESHOLD_LEVEL = "threshold_level"
    SHIFT_HOURS = "shift_hours"
    NAME = "name"