import sys
from PySide6.QtWidgets import QApplication
from vslm.gui import MainWindow

def main():
    # Initialize the Application
    app = QApplication(sys.argv)
    app.setStyle('Fusion') # Optional: Fusion style looks good on all platforms
    
    # Create and Show the Main Window
    window = MainWindow()
    window.show()
    
    # Start the Event Loop
    sys.exit(app.exec())

if __name__ == "__main__":
    main()