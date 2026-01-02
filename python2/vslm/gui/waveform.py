# python2/vslm/gui/waveform.py
import numpy as np
import soundfile as sf
import pyqtgraph as pg
from PySide6.QtWidgets import (QDialog, QVBoxLayout, QPushButton, 
                               QHBoxLayout, QLabel, QDialogButtonBox, QWidget)
from PySide6.QtCore import Qt, Signal

class WaveformViewer(QWidget):
    """
    Widget to display the amplitude envelope of a large audio file 
    and allow time-range selection.
    """
    selection_changed = Signal(float, float)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.layout = QVBoxLayout(self)
        self.layout.setContentsMargins(0, 0, 0, 0)
        
        # PyQtGraph Setup
        self.plot_widget = pg.PlotWidget()
        self.plot_widget.setBackground('w')
        self.plot_widget.showGrid(x=True, y=True, alpha=0.3)
        self.plot_widget.setLabel('left', 'Amplitude')
        self.plot_widget.setLabel('bottom', 'Time', units='s')
        self.plot_widget.setMouseEnabled(x=True, y=False) # Only pan/zoom time
        self.layout.addWidget(self.plot_widget)
        
        # Data Containers
        self.time_axis = None
        self.min_envelope = None
        self.max_envelope = None
        self.duration = 0.0
        
        # Selection Region (The "Audacity" selector)
        self.region = pg.LinearRegionItem()
        self.region.setZValue(10)
        self.region.sigRegionChanged.connect(self._on_region_changed)
        self.plot_widget.addItem(self.region)
        
        # Curves
        # We fill between min and max to create the "solid" waveform look
        self.curve_min = self.plot_widget.plot(pen=pg.mkPen('#1f77b4', width=1))
        self.curve_max = self.plot_widget.plot(pen=pg.mkPen('#1f77b4', width=1))
        self.fill = pg.FillBetweenItem(self.curve_min, self.curve_max, brush=pg.mkBrush('#1f77b450'))
        self.plot_widget.addItem(self.fill)

    def load_file(self, filepath):
        """
        Reads the file and computes a decimated min/max envelope for efficient plotting.
        """
        info = sf.info(filepath)
        self.duration = info.duration
        fs = info.samplerate
        
        # Target ~10,000 points for the visualization (enough for 4k screens)
        target_points = 10000
        total_frames = info.frames
        
        # Calculate step size to achieve target points
        step = max(1, total_frames // target_points)
        
        # Arrays to hold the envelope
        mins = []
        maxs = []
        
        # Read and decimate
        # Using blocks keeps RAM usage low during load
        with sf.SoundFile(filepath) as f:
            while f.tell() < total_frames:
                # Read a chunk (e.g., 10 steps worth)
                chunk_frames = step * 10
                data = f.read(chunk_frames, always_2d=True)
                
                # If stereo, mix to mono for visualization
                if data.shape[1] > 1:
                    data = np.mean(data, axis=1)
                else:
                    data = data.flatten()
                
                # Compute min/max for each 'step' in this chunk
                # Reshape to (-1, step) to vectorize the min/max calc
                # Handle edge case where data len isn't multiple of step
                trim = len(data) % step
                if trim > 0:
                    data = data[:-trim]
                
                if len(data) > 0:
                    reshaped = data.reshape(-1, step)
                    mins.append(np.min(reshaped, axis=1))
                    maxs.append(np.max(reshaped, axis=1))
        
        self.min_envelope = np.concatenate(mins)
        self.max_envelope = np.concatenate(maxs)
        
        # Time axis
        num_points = len(self.min_envelope)
        self.time_axis = np.linspace(0, self.duration, num_points)
        
        # Update Plot
        self.curve_min.setData(self.time_axis, self.min_envelope)
        self.curve_max.setData(self.time_axis, self.max_envelope)
        
        # Reset Region to full file
        self.region.setRegion([0, self.duration])
        self.plot_widget.autoRange()

    def _on_region_changed(self):
        min_t, max_t = self.region.getRegion()
        self.selection_changed.emit(min_t, max_t)

    def get_selection(self):
        return self.region.getRegion()


class WaveformDialog(QDialog):
    """
    Modal dialog containing the WaveformViewer.
    """
    def __init__(self, filepath, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Select File Section")
        self.resize(900, 400)
        
        layout = QVBoxLayout(self)
        
        # Viewer
        self.viewer = WaveformViewer()
        self.viewer.load_file(filepath)
        layout.addWidget(self.viewer)
        
        # Info Label
        self.lbl_info = QLabel("Selected: 0.00s - 0.00s (Duration: 0.00s)")
        self.lbl_info.setAlignment(Qt.AlignCenter)
        layout.addWidget(self.lbl_info)
        self.viewer.selection_changed.connect(self.update_label)
        
        # Buttons
        btns = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btns.accepted.connect(self.accept)
        btns.rejected.connect(self.reject)
        layout.addWidget(btns)
        
        # Initialize label
        start, end = self.viewer.get_selection()
        self.update_label(start, end)

    def update_label(self, start, end):
        # Clamp to bounds
        start = max(0, start)
        end = min(self.viewer.duration, end)
        self.lbl_info.setText(f"Selected: {start:.2f}s - {end:.2f}s (Duration: {end-start:.2f}s)")

    def get_selection(self):
        start, end = self.viewer.get_selection()
        # Clamp
        start = max(0, start)
        end = min(self.viewer.duration, end)
        return start, end