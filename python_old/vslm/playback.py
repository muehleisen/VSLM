# vslm/playback.py
import sounddevice as sd
import numpy as np

class AudioPlayer:
    def __init__(self):
        self.stream = None
        self.is_playing = False

    def play(self, data, fs):
        if self.is_playing:
            self.stop()
            
        self.is_playing = True
        
        # Normalize for playback (don't blast calibrated levels!)
        # Calibrated data might be > 1.0 (Pascals). Audio cards want -1.0 to 1.0.
        # Find max and normalize just for hearing purposes.
        mx = np.max(np.abs(data))
        if mx > 1.0:
            play_data = data / mx
        else:
            play_data = data
            
        sd.play(play_data, fs)

    def stop(self):
        sd.stop()
        self.is_playing = False