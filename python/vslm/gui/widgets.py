# python2/vslm/gui/widgets.py
from PySide6.QtWidgets import QWidget, QVBoxLayout
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

class MatplotlibWidget(QWidget):
    """
    A simple wrapper to display Matplotlib plots in PySide6.
    """
    def __init__(self, parent=None):
        super().__init__(parent)
        self.layout = QVBoxLayout(self)
        self.layout.setContentsMargins(0,0,0,0)
        
        self.figure = Figure(figsize=(5, 4), dpi=100)
        self.canvas = FigureCanvas(self.figure)
        self.layout.addWidget(self.canvas)
        
        self.ax = self.figure.add_subplot(111)
        self.ax.text(0.5, 0.5, "Ready", ha='center', va='center')
        self.ax.axis('off')

    def reset(self):
        self.figure.clear()
        self.ax = self.figure.add_subplot(111)

    def draw(self):
        self.canvas.draw()