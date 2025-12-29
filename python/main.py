# main.py
import sys
from PyQt6.QtWidgets import QApplication
from vslm.gui import MainWindow

def main():
    # Create the Application
    app = QApplication(sys.argv)
    app.setApplicationName("VSLM Python")

    # Create and Show the Window
    window = MainWindow()
    window.show()

    # Run the Event Loop
    sys.exit(app.exec())

if __name__ == "__main__":
    main()