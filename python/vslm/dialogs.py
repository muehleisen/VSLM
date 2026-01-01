# vslm/dialogs.py
from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QFormLayout, QDialogButtonBox, 
    QLabel, QSpinBox, QDoubleSpinBox, QComboBox, 
    QMessageBox, QCheckBox
)

class BaseSettingsDialog(QDialog):
    """Base class for settings dialogs with OK/Cancel buttons."""
    def __init__(self, title, parent=None):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.layout = QVBoxLayout(self)
        self.form_layout = QFormLayout()
        self.layout.addLayout(self.form_layout)
        
        self.buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        self.layout.addWidget(self.buttons)

class AboutDialog:
    """Wrapper for the About message box."""
    @staticmethod
    def show(parent):
        QMessageBox.about(
            parent, 
            "About VSLM", 
            "Virtual Sound Level Meter (VSLM)\n"
            "Python Port Version 1.0\n\n"
            "Original MATLAB Source by Ralph Muehleisen.\n"
            "Ported to Python/PySide6."
        )

class CalibrationDialog(BaseSettingsDialog):
    def __init__(self, current_val=94.0, parent=None):
        super().__init__("Set Calibration", parent)
        self.spin_level = QDoubleSpinBox()
        self.spin_level.setRange(0.0, 194.0)
        self.spin_level.setDecimals(2)
        self.spin_level.setValue(current_val)
        self.spin_level.setSuffix(" dB")
        self.form_layout.addRow("Calibrator Level:", self.spin_level)

    def get_value(self):
        return self.spin_level.value()

# --- Plotting & Scales ---

class PlotTimeSpacingDialog(BaseSettingsDialog):
    def __init__(self, current_val=1.0, parent=None):
        super().__init__("Plot Time Spacing", parent)
        self.spin_time = QDoubleSpinBox()
        self.spin_time.setRange(0.01, 3600.0)
        self.spin_time.setDecimals(2)
        self.spin_time.setValue(current_val)
        self.spin_time.setSuffix(" s")
        self.form_layout.addRow("Major Tick Spacing:", self.spin_time)

    def get_value(self):
        return self.spin_time.value()

class PlotScalesDialog(BaseSettingsDialog):
    def __init__(self, current_min=30, current_max=120, is_autoscale=False, parent=None):
        super().__init__("Plot Scales", parent)
        
        # Autoscale Checkbox
        self.chk_autoscale = QCheckBox("Autoscale")
        self.chk_autoscale.setChecked(is_autoscale)
        self.chk_autoscale.toggled.connect(self.toggle_inputs)
        
        # Spinboxes
        self.spin_min = QSpinBox()
        self.spin_min.setRange(-200, 200)
        self.spin_min.setValue(current_min)
        self.spin_min.setSuffix(" dB")
        
        self.spin_max = QSpinBox()
        self.spin_max.setRange(-200, 200)
        self.spin_max.setValue(current_max)
        self.spin_max.setSuffix(" dB")
        
        self.form_layout.addRow("", self.chk_autoscale)
        self.form_layout.addRow("Y-Axis Min:", self.spin_min)
        self.form_layout.addRow("Y-Axis Max:", self.spin_max)
        
        # Initialize State
        self.toggle_inputs(is_autoscale)

    def toggle_inputs(self, checked):
        self.spin_min.setEnabled(not checked)
        self.spin_max.setEnabled(not checked)

    def get_values(self):
        return self.spin_min.value(), self.spin_max.value(), self.chk_autoscale.isChecked()

# --- Leq Settings ---

class LeqIntegrationDialog(BaseSettingsDialog):
    def __init__(self, current_val=1.0, parent=None):
        super().__init__("LEQ Integration Time", parent)
        self.spin_time = QDoubleSpinBox()
        self.spin_time.setRange(0.1, 86400.0)
        self.spin_time.setValue(current_val)
        self.spin_time.setSuffix(" s")
        self.form_layout.addRow("Integration Interval:", self.spin_time)

    def get_value(self):
        return self.spin_time.value()

class LeqPercentileDialog(BaseSettingsDialog):
    def __init__(self, current_val=90, parent=None):
        super().__init__("LEQ Percentile", parent)
        self.spin_perc = QSpinBox()
        self.spin_perc.setRange(1, 99)
        self.spin_perc.setValue(current_val)
        self.spin_perc.setPrefix("L")
        self.form_layout.addRow("Percentile Exceeded:", self.spin_perc)

    def get_value(self):
        return self.spin_perc.value()

class NoiseDoseDialog(BaseSettingsDialog):
    def __init__(self, exchange=3, threshold=80, criterion=90, parent=None):
        super().__init__("Noise Dose Criterion", parent)
        self.combo_exchange = QComboBox()
        self.combo_exchange.addItems(["3 dB", "4 dB", "5 dB"])
        idx = max(0, exchange - 3)
        self.combo_exchange.setCurrentIndex(idx if idx < 3 else 0)
        self.spin_thresh = QSpinBox()
        self.spin_thresh.setRange(0, 140)
        self.spin_thresh.setValue(threshold)
        self.spin_thresh.setSuffix(" dB")
        self.spin_crit = QSpinBox()
        self.spin_crit.setRange(0, 140)
        self.spin_crit.setValue(criterion)
        self.spin_crit.setSuffix(" dB")
        self.form_layout.addRow("Exchange Rate:", self.combo_exchange)
        self.form_layout.addRow("Threshold Level:", self.spin_thresh)
        self.form_layout.addRow("Criterion Level:", self.spin_crit)

    def get_values(self):
        rate = self.combo_exchange.currentIndex() + 3
        return rate, self.spin_thresh.value(), self.spin_crit.value()

# --- Band Settings ---

class BandResolutionDialog(BaseSettingsDialog):
    def __init__(self, current_idx=0, parent=None):
        super().__init__("Band Resolution", parent)
        self.combo_res = QComboBox()
        self.combo_res.addItems(["1/1 Octave", "1/3 Octave", "1/6 Octave", "1/12 Octave"])
        self.combo_res.setCurrentIndex(current_idx)
        self.form_layout.addRow("Resolution:", self.combo_res)

    def get_value(self):
        return self.combo_res.currentText()

class BandMethodDialog(BaseSettingsDialog):
    def __init__(self, current_idx=0, parent=None):
        super().__init__("Band Calculation Method", parent)
        self.combo_method = QComboBox()
        self.combo_method.addItems(["ANSI S1.11 (Time Domain)", "FFT Synthesis (Freq Domain)"])
        self.combo_method.setCurrentIndex(current_idx)
        self.form_layout.addRow("Method:", self.combo_method)

    def get_value(self):
        return self.combo_method.currentText()

# --- PSD & Spectrogram Settings ---

class FFTSizeDialog(BaseSettingsDialog):
    def __init__(self, current_size=4096, parent=None):
        super().__init__("FFT Size", parent)
        self.combo_size = QComboBox()
        sizes = [str(2**i) for i in range(8, 17)]
        self.combo_size.addItems(sizes)
        idx = self.combo_size.findText(str(current_size))
        if idx >= 0: self.combo_size.setCurrentIndex(idx)
        else: self.combo_size.setCurrentIndex(4)
        self.form_layout.addRow("Points (N):", self.combo_size)

    def get_value(self):
        return int(self.combo_size.currentText())

class OverlapDialog(BaseSettingsDialog):
    def __init__(self, current_val=50, parent=None):
        super().__init__("FFT Overlap", parent)
        self.spin_ov = QSpinBox()
        self.spin_ov.setRange(0, 99)
        self.spin_ov.setValue(current_val)
        self.spin_ov.setSuffix(" %")
        self.form_layout.addRow("Overlap Percentage:", self.spin_ov)

    def get_value(self):
        return self.spin_ov.value()

class WindowDialog(BaseSettingsDialog):
    def __init__(self, current_win="Hann", parent=None):
        super().__init__("Window Function", parent)
        self.combo_win = QComboBox()
        self.combo_win.addItems(["Hann", "Hamming", "Blackman", "Bartlett", "Rectangular"])
        idx = self.combo_win.findText(current_win)
        if idx >= 0: self.combo_win.setCurrentIndex(idx)
        self.form_layout.addRow("Window Type:", self.combo_win)

    def get_value(self):
        return self.combo_win.currentText()

class SliceLengthDialog(BaseSettingsDialog):
    def __init__(self, current_val=0.1, parent=None):
        super().__init__("Slice Length", parent)
        self.spin_slice = QDoubleSpinBox()
        self.spin_slice.setRange(0.001, 10.0)
        self.spin_slice.setDecimals(3)
        self.spin_slice.setValue(current_val)
        self.spin_slice.setSuffix(" s")
        self.form_layout.addRow("Time Slice:", self.spin_slice)

    def get_value(self):
        return self.spin_slice.value()

class SpectrogramViewDialog(BaseSettingsDialog):
    """Dialog for Spectrogram 'View' settings (Colormap only)."""
    def __init__(self, cmap="inferno", parent=None):
        super().__init__("Spectrogram View", parent)
        
        self.combo_cmap = QComboBox()
        # Common scientific colormaps
        cmaps = [
            "inferno", "magma", "plasma", "viridis", 
            "cividis", "twilight", 
            "jet", "turbo", "rainbow",               
            "ocean", "gist_earth", "terrain",
            "gray", "bone", "hot", "cool"
        ]
        self.combo_cmap.addItems(cmaps)
        
        idx = self.combo_cmap.findText(cmap)
        if idx >= 0: 
            self.combo_cmap.setCurrentIndex(idx)
        else:
            self.combo_cmap.setCurrentIndex(0)
        
        self.form_layout.addRow("Color Map:", self.combo_cmap)

    def get_value(self):
        return self.combo_cmap.currentText()