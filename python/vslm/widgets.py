# vslm/widgets.py
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QLabel, QPushButton, QSizePolicy
)
from PySide6.QtCore import Qt, QSize
from PySide6.QtGui import QCursor

# Matplotlib Imports
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
# Essential for '3d' projection to work
from mpl_toolkits.mplot3d import Axes3D 

class VSLMPlotWidget(QWidget):
    """
    A self-contained widget handling the Matplotlib Figure, Canvas, and Toolbar.
    """
    def __init__(self, parent=None):
        super().__init__(parent)
        
        # Create Layout
        self.layout = QVBoxLayout(self)
        self.layout.setContentsMargins(0, 0, 0, 0)
        
        # Initialize Matplotlib Figure
        self.figure = Figure(figsize=(5, 4), dpi=100)
        self.canvas = FigureCanvas(self.figure)
        
        # Add Toolbar
        self.toolbar = NavigationToolbar(self.canvas, self)
        
        # Add to Layout
        self.layout.addWidget(self.toolbar)
        self.layout.addWidget(self.canvas)
        
        # Initialize default view
        self.show_placeholder()

    def show_placeholder(self):
        """Clears the figure and shows the placeholder text (Axis OFF)."""
        self.figure.clear()
        self.ax = self.figure.add_subplot(111)
        self.ax.axis('off')
        self.ax.text(0.5, 0.5, "Load a file to begin", ha='center', va='center')
        self.canvas.draw()

    def prepare_plot(self, projection=None):
        """
        Clears the figure and prepares for data plotting.
        Args:
            projection (str): '3d' for 3D plots, None for standard 2D.
        """
        self.figure.clear()
        self.ax = self.figure.add_subplot(111, projection=projection)
        if not projection: 
            self.ax.axis('on') 
        return self.ax

    def draw(self):
        """Forces a canvas redraw."""
        self.canvas.draw()


class LabelledToggleButton(QWidget):
    """
    A custom widget with a Label stacked vertically above a Rectangular Toggle Button.
    """
    BTN_STYLE = """
        QPushButton {
            border: 2px solid #aaa;
            border-radius: 8px;
            background-color: #f5f5f5;
        }
        QPushButton:checked {
            background-color: #3b82f6;
            border-color: #1d4ed8;
        }
        QPushButton:hover {
            border-color: #3b82f6;
        }
    """
    
    def __init__(self, label_text, btn_id, button_group, parent=None):
        super().__init__(parent)
        
        self.layout = QVBoxLayout(self)
        self.layout.setContentsMargins(0, 5, 0, 5)
        self.layout.setSpacing(4)
        
        self.label = QLabel(label_text)
        self.label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.label.setStyleSheet("font-weight: bold; color: #444; font-size: 11px;")
        
        self.button = QPushButton("")
        self.button.setCheckable(True)
        self.button.setFixedSize(QSize(50, 30))
        self.button.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self.button.setStyleSheet(self.BTN_STYLE)
        
        if button_group:
            button_group.addButton(self.button, btn_id)
            
        self.layout.addWidget(self.label)
        self.layout.addWidget(self.button, alignment=Qt.AlignmentFlag.AlignCenter)

    def setChecked(self, checked):
        self.button.setChecked(checked)