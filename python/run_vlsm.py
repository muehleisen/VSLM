# python2/run_vslm.py
import sys
import os
from PySide6.QtWidgets import QApplication

# Ensure the current directory is strictly in the Python path
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

from vslm.gui.main_window import MainWindow

if __name__ == "__main__":
    app = QApplication(sys.argv)
    
    # Optional: Set a clean fusion style for Windows/Linux consistency
    app.setStyle("Fusion")
    
    window = MainWindow()
    window.show()
    
    sys.exit(app.exec())