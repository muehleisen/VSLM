# python2/vslm/gui/about_dialog.py
from PySide6.QtWidgets import QDialog, QVBoxLayout, QLabel, QPushButton, QFrame
from PySide6.QtCore import Qt

class AboutDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("About VSLM")
        self.resize(300, 250)
        self._init_ui()

    def _init_ui(self):
        layout = QVBoxLayout(self)
        
        # App Name & Version
        lbl_title = QLabel("VSLM 2.0 (Python)")
        lbl_title.setStyleSheet("font-size: 18px; font-weight: bold; color: #1e3a8a;")
        lbl_title.setAlignment(Qt.AlignCenter)
        layout.addWidget(lbl_title)
        
        lbl_ver = QLabel("Version 2.0.0-alpha")
        lbl_ver.setStyleSheet("font-size: 12px; color: #555;")
        lbl_ver.setAlignment(Qt.AlignCenter)
        layout.addWidget(lbl_ver)
        
        # Horizontal Line
        line = QFrame()
        line.setFrameShape(QFrame.HLine)
        line.setFrameShadow(QFrame.Sunken)
        layout.addWidget(line)
        
        # Credits / License
        credits_text = (
            "<b>License:</b> MIT License<br><br>"
            "<b>Credits:</b><br>"
            "Original MATLAB Code: [Name]<br>"
            "Python Port: [Name]<br>"
            "GUI Framework: PySide6 (Qt)<br>"
            "Plotting: Matplotlib"
        )
        lbl_credits = QLabel(credits_text)
        lbl_credits.setAlignment(Qt.AlignCenter)
        lbl_credits.setWordWrap(True)
        layout.addWidget(lbl_credits)
        
        layout.addStretch()
        
        # Close Button
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.accept)
        layout.addWidget(btn_close)