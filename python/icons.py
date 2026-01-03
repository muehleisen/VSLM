import sys
from PySide6.QtWidgets import QApplication, QWidget, QGridLayout, QLabel, QStyle
from PySide6.QtCore import Qt

class IconPreview(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Qt Standard Icons")
        layout = QGridLayout(self)
        
        # Get all StandardPixmap enums
        style = self.style()
        meta = QStyle.StandardPixmap
        
        # Loop through common enums (0 to ~70)
        row, col = 0, 0
        for i in range(70):
            try:
                # Convert integer to Enum member if possible for name, 
                # or just use the integer if specific enum mapping is tricky in loop
                enum_val = QStyle.StandardPixmap(i)
                name = enum_val.name.replace("SP_", "")
                
                icon = style.standardIcon(enum_val)
                if not icon.isNull():
                    lbl_icon = QLabel()
                    lbl_icon.setPixmap(icon.pixmap(32, 32))
                    lbl_icon.setAlignment(Qt.AlignCenter)
                    
                    layout.addWidget(lbl_icon, row, col)
                    layout.addWidget(QLabel(name), row + 1, col)
                    
                    col += 1
                    if col > 6:
                        col = 0
                        row += 2
            except:
                pass

if __name__ == "__main__":
    app = QApplication(sys.argv)
    w = IconPreview()
    w.show()
    sys.exit(app.exec())