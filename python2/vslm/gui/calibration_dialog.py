# python2/vslm/gui/calibration_dialog.py
from PySide6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QLabel, 
                               QLineEdit, QPushButton, QGroupBox, QDoubleSpinBox,
                               QMessageBox, QFormLayout)
from PySide6.QtCore import Qt
from pathlib import Path

from ..calibration import compute_selection_rms, calculate_factor_from_ref

class CalibrationDialog(QDialog):
    def __init__(self, current_factor: float, filepath: Path = None, start: float = 0.0, end: float = 0.0, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Calibration Manager")
        self.resize(400, 350)
        
        self.result_factor = current_factor
        self.filepath = filepath
        self.selection_start = start
        self.selection_end = end
        
        self._init_ui()
        
    def _init_ui(self):
        layout = QVBoxLayout(self)
        
        # --- Current Status ---
        self.lbl_current = QLabel(f"Current Calibration Factor: {self.result_factor:.4f}")
        self.lbl_current.setStyleSheet("font-size: 14px; font-weight: bold; color: #1e40af; margin-bottom: 10px;")
        self.lbl_current.setAlignment(Qt.AlignCenter)
        layout.addWidget(self.lbl_current)
        
        # --- Method 1: Manual Entry ---
        grp_manual = QGroupBox("Method 1: Manual Entry")
        lay_manual = QHBoxLayout()
        
        self.spin_manual = QDoubleSpinBox()
        self.spin_manual.setRange(0.0001, 10000.0)
        self.spin_manual.setDecimals(4)
        self.spin_manual.setValue(self.result_factor)
        self.spin_manual.setSuffix(" (Linear)")
        
        btn_apply_manual = QPushButton("Apply Manual")
        btn_apply_manual.clicked.connect(self.on_apply_manual)
        
        lay_manual.addWidget(QLabel("Factor:"))
        lay_manual.addWidget(self.spin_manual)
        lay_manual.addWidget(btn_apply_manual)
        grp_manual.setLayout(lay_manual)
        layout.addWidget(grp_manual)
        
        # --- Method 2: From Selection ---
        grp_ref = QGroupBox("Method 2: From Selected Audio")
        lay_ref = QVBoxLayout()
        
        info_txt = "No file loaded or selection too short."
        has_valid_selection = False
        if self.filepath and (self.selection_end - self.selection_start > 0.1):
            dur = self.selection_end - self.selection_start
            info_txt = f"Selection: {self.selection_start:.2f}s - {self.selection_end:.2f}s ({dur:.2f}s)"
            has_valid_selection = True
            
        lay_ref.addWidget(QLabel(info_txt))
        
        form_ref = QFormLayout()
        self.spin_ref_level = QDoubleSpinBox()
        self.spin_ref_level.setRange(0.0, 194.0)
        self.spin_ref_level.setValue(94.0) # Standard Pistonphone
        self.spin_ref_level.setSuffix(" dB")
        
        form_ref.addRow("Reference Level:", self.spin_ref_level)
        lay_ref.addLayout(form_ref)
        
        self.btn_calc_ref = QPushButton("Calculate & Apply")
        self.btn_calc_ref.setStyleSheet("background-color: #dcfce7; font-weight: bold;")
        self.btn_calc_ref.clicked.connect(self.on_calculate_ref)
        self.btn_calc_ref.setEnabled(has_valid_selection)
        
        lay_ref.addWidget(self.btn_calc_ref)
        grp_ref.setLayout(lay_ref)
        layout.addWidget(grp_ref)
        
        # --- Close ---
        layout.addStretch()
        btn_close = QPushButton("Done")
        btn_close.clicked.connect(self.accept)
        layout.addWidget(btn_close)

    def on_apply_manual(self):
        self.result_factor = self.spin_manual.value()
        self.update_display()
        QMessageBox.information(self, "Updated", f"Calibration factor set to {self.result_factor:.4f}")

    def on_calculate_ref(self):
        try:
            # 1. Compute RMS of uncalibrated file selection
            rms = compute_selection_rms(self.filepath, self.selection_start, self.selection_end)
            
            # 2. Calculate needed factor
            target_db = self.spin_ref_level.value()
            new_factor = calculate_factor_from_ref(rms, target_db)
            
            # 3. Update UI
            self.result_factor = new_factor
            self.spin_manual.setValue(new_factor)
            self.update_display()
            
            QMessageBox.information(self, "Success", 
                                    f"Measured RMS: {rms:.6f}\n"
                                    f"Target: {target_db} dB\n"
                                    f"New Factor: {new_factor:.4f}")
            
        except Exception as e:
            QMessageBox.critical(self, "Calibration Failed", str(e))

    def update_display(self):
        self.lbl_current.setText(f"Current Calibration Factor: {self.result_factor:.4f}")
    
    def get_factor(self):
        return self.result_factor