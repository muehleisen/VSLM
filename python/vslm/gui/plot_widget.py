import matplotlib
matplotlib.use('QtAgg') 
from PySide6.QtWidgets import (QWidget, QVBoxLayout, QDialog, QCheckBox, 
                               QDoubleSpinBox, QHBoxLayout, QPushButton,
                               QFormLayout, QToolButton)
from PySide6.QtGui import QAction, QIcon, QPainter, QPen, QPixmap, QColor, QFont
from PySide6.QtCore import Signal, QSize, Qt
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

class PlotSettingsDialog(QDialog):
    """Dialog to set plot scaling limits."""
    def __init__(self, autoscale, ymin, ymax, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Plot Scaling")
        self.resize(300, 150)
        
        layout = QVBoxLayout(self)
        
        self.chk_auto = QCheckBox("Autoscale Y-Axis")
        self.chk_auto.setChecked(autoscale)
        self.chk_auto.toggled.connect(self._toggle_inputs)
        layout.addWidget(self.chk_auto)
        
        form = QFormLayout()
        self.spin_min = QDoubleSpinBox()
        self.spin_min.setRange(-200, 200)
        self.spin_min.setValue(ymin)
        form.addRow("Min dB:", self.spin_min)
        
        self.spin_max = QDoubleSpinBox()
        self.spin_max.setRange(-200, 200)
        self.spin_max.setValue(ymax)
        form.addRow("Max dB:", self.spin_max)
        layout.addLayout(form)
        
        # Dialog Buttons
        btns = QHBoxLayout()
        ok = QPushButton("OK")
        ok.clicked.connect(self.accept)
        cancel = QPushButton("Cancel")
        cancel.clicked.connect(self.reject)
        btns.addWidget(ok)
        btns.addWidget(cancel)
        layout.addLayout(btns)
        
        self._toggle_inputs(autoscale)

    def _toggle_inputs(self, checked):
        self.spin_min.setEnabled(not checked)
        self.spin_max.setEnabled(not checked)
        
    def get_values(self):
        return self.chk_auto.isChecked(), self.spin_min.value(), self.spin_max.value()

class CustomToolbar(NavigationToolbar):
    """Custom Toolbar that adds a Scaling button using a Unicode character."""
    sig_open_scaling = Signal()

    def __init__(self, canvas, parent=None):
        super().__init__(canvas, parent)
        
        # 1. Create the Custom Button
        self.btn_scale = QToolButton(self)
        self.btn_scale.setToolTip("Configure Plot Scaling (Min/Max)")
        
        # 2. Set the Unicode character (⇕) and make it bold
        self.btn_scale.setText("⇕") 
        font = self.btn_scale.font()
        font.setBold(True)
        
        # Use QFont.Weight enum for PySide6 compatibility
        font.setWeight(QFont.Weight.Black) 
        
        # Updated font size to 20
        font.setPointSize(24) 
        self.btn_scale.setFont(font)
        
        self.btn_scale.clicked.connect(self.sig_open_scaling.emit)

        # 3. Position logic
        target_action = None
        for action in self.actions():
            if self.widgetForAction(action) == self.locLabel:
                target_action = action
                break
        
        if target_action:
            self.insertSeparator(target_action)
            self.insertWidget(target_action, self.btn_scale)
        else:
            self.addSeparator()
            self.addWidget(self.btn_scale)

class MatplotlibWidget(QWidget):
    # Emits (autoscale, ymin, ymax) when user changes settings in the dialog
    sig_scaling_changed = Signal(bool, float, float)

    def __init__(self, parent=None):
        super().__init__(parent)
        
        # Use local variable 'layout' to avoid shadowing QWidget.layout()
        layout = QVBoxLayout(self) 
        
        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        
        # Use Custom Toolbar
        self.toolbar = CustomToolbar(self.canvas, self)
        self.toolbar.sig_open_scaling.connect(self.open_scaling_dialog)
        
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        
        # Internal state for the dialog
        self._autoscale = True
        self._ymin = 0.0
        self._ymax = 120.0

    def set_plot_settings(self, autoscale, ymin, ymax):
        """Called by MainWindow to initialize state."""
        self._autoscale = autoscale
        self._ymin = ymin
        self._ymax = ymax

    def open_scaling_dialog(self):
        dlg = PlotSettingsDialog(self._autoscale, self._ymin, self._ymax, self)
        if dlg.exec():
            # Update local state and notify MainWindow
            auto, ymin, ymax = dlg.get_values()
            self._autoscale = auto
            self._ymin = ymin
            self._ymax = ymax
            self.sig_scaling_changed.emit(auto, ymin, ymax)

    def draw(self):
        self.canvas.draw()