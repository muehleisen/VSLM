import matplotlib
matplotlib.use('QtAgg') # Ensure QtAgg backend
from PySide6.QtWidgets import QWidget, QVBoxLayout
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

class MatplotlibWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.layout = QVBoxLayout(self)
        
        # Create the Figure and Canvas
        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        
        # Create the Toolbar
        # The toolbar requires the canvas and the parent widget
        self.toolbar = NavigationToolbar(self.canvas, self)
        
        # Add widgets to layout: Toolbar first (top), then Canvas
        self.layout.addWidget(self.toolbar)
        self.layout.addWidget(self.canvas)

    def draw(self):
        self.canvas.draw()