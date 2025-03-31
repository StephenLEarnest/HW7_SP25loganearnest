import sys
from pyXSteam.XSteam import XSteam
from PyQt5.QtWidgets import (QApplication, QWidget, QVBoxLayout, QHBoxLayout,
                             QGroupBox, QLabel, QComboBox, QLineEdit, QPushButton,
                             QRadioButton)
from scipy.optimize import fsolve


# Minimal UnitConversion module (replace with your UC if available)
class UC:
    psi_to_bar = 0.0689476
    bar_to_psi = 14.5038

    @staticmethod
    def F_to_C(f): return (f - 32) * 5 / 9

    @staticmethod
    def C_to_F(c): return c * 9 / 5 + 32

    btuperlb_to_kJperkg = 2.326
    kJperkg_to_btuperlb = 0.429923
    btuperlbF_to_kJperkgC = 4.1868
    kJperkgC_to_btuperlbF = 0.2388459
    ft3perlb_to_m3perkg = 0.062428
    m3perkg_to_ft3perlb = 16.0185


# Reusing thermoSatProps and thermoState from your original code (simplified here)
class thermoSatProps:
    def __init__(self, p=None, t=None, SI=True):
        self.steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS if SI else XSteam.UNIT_SYSTEM_FLS)
        if p is not None:
            self.getSatProps(p, SI)
        elif t is not None:
            self.getSatProps(self.steamTable.psat_t(t), SI)

    def getSatProps(self, p, SI=True):
        self.pSat = p
        self.tSat = self.steamTable.tsat_p(p)
        self.vf = self.steamTable.vL_p(p)
        self.vg = self.steamTable.vV_p(p)
        self.hf = self.steamTable.hL_p(p)
        self.hg = self.steamTable.hV_p(p)
        self.uf = self.steamTable.uL_p(p)
        self.ug = self.steamTable.uV_p(p)
        self.sf = self.steamTable.sL_p(p)
        self.sg = self.steamTable.sV_p(p)
        self.vgf = self.vg - self.vf
        self.hgf = self.hg - self.hf
        self.sgf = self.sg - self.sf
        self.ugf = self.ug - self.uf


class thermoState:
    def __init__(self, p=None, t=None, v=None, u=None, h=None, s=None, x=None):
        self.steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS)
        self.region = "saturated"
        self.p = p
        self.t = t
        self.v = v
        self.u = u
        self.h = h
        self.s = s
        self.x = x

    def computeProperties(self):
        if self.region == "two-phase":
            self.u = self.steamTable.uL_p(self.p) + self.x * (
                        self.steamTable.uV_p(self.p) - self.steamTable.uL_p(self.p))
            self.h = self.steamTable.hL_p(self.p) + self.x * (
                        self.steamTable.hV_p(self.p) - self.steamTable.hL_p(self.p))
            self.s = self.steamTable.sL_p(self.p) + self.x * (
                        self.steamTable.sV_p(self.p) - self.steamTable.sL_p(self.p))
            self.v = self.steamTable.vL_p(self.p) + self.x * (
                        self.steamTable.vV_p(self.p) - self.steamTable.vL_p(self.p))
        else:
            self.u = self.steamTable.u_pt(self.p, self.t)
            self.h = self.steamTable.h_pt(self.p, self.t)
            self.s = self.steamTable.s_pt(self.p, self.t)
            self.v = self.steamTable.v_pt(self.p, self.t)
            self.x = 1.0 if self.region == "super-heated vapor" else 0.0

    def setState(self, stProp1, stProp2, stPropVal1, stPropVal2, SI=True):
        self.steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS if SI else XSteam.UNIT_SYSTEM_FLS)
        SP = [stProp1.lower(), stProp2.lower()]
        f1, f2 = float(stPropVal1), float(stPropVal2)

        if SP[0] == 'p' or SP[1] == 'p':
            oFlipped = SP[0] != 'p'
            SP1 = SP[0] if oFlipped else SP[1]
            self.p = f1 if not oFlipped else f2
            tSat = self.steamTable.tsat_p(self.p)
            if SP1 == 't':
                self.t = f2 if not oFlipped else f1
                if self.t < tSat:
                    self.region = "sub-cooled liquid"
                elif self.t > tSat:
                    self.region = "super-heated vapor"
                else:
                    self.region = "two-phase"
                    self.x = 0.5
            elif SP1 == 'v':
                self.v = f2 if not oFlipped else f1
                vf, vg = self.steamTable.vL_p(self.p), self.steamTable.vV_p(self.p)
                if self.v < vf:
                    self.region = "sub-cooled liquid"
                    self.t = fsolve(lambda T: self.v - self.steamTable.v_pt(self.p, T), tSat - 1)[0]
                elif self.v > vg:
                    self.region = "super-heated vapor"
                    self.t = fsolve(lambda T: self.v - self.steamTable.v_pt(self.p, T), tSat + 1)[0]
                else:
                    self.region = "two-phase"
                    self.x = (self.v - vf) / (vg - vf)
                    self.t = tSat
            # Add remaining cases (ph, pu, ps, px) as needed from your original code
        # Placeholder: Implement other 20 cases similarly (Tv, Th, etc.)
        self.computeProperties()

    def __sub__(self, other):
        delta = thermoState()
        delta.p = self.p - other.p
        delta.t = self.t - other.t
        delta.h = self.h - other.h
        delta.u = self.u - other.u
        delta.s = self.s - other.s
        delta.v = self.v - other.v
        return delta


class MainWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS)
        self.currentUnits = 'SI'
        self.initUI()
        self.setUnits()
        self.show()

    def initUI(self):
        self.setWindowTitle("Thermodynamic State Calculator")
        main_layout = QVBoxLayout()

        # Units Selection
        units_box = QGroupBox("Units")
        units_layout = QHBoxLayout()
        self.rdo_SI = QRadioButton("SI")
        self.rdo_English = QRadioButton("English")
        self.rdo_SI.setChecked(True)
        units_layout.addWidget(self.rdo_SI)
        units_layout.addWidget(self.rdo_English)
        units_box.setLayout(units_layout)
        main_layout.addWidget(units_box)

        # Specified Properties
        spec_box = QGroupBox("Specified Properties")
        spec_layout = QVBoxLayout()

        # State 1
        self.state1_box = QGroupBox("State 1")
        state1_layout = QVBoxLayout()
        self.cmb_s1_p1 = QComboBox()
        self.le_s1_p1 = QLineEdit()
        self.lbl_s1_p1_units = QLabel()
        self.cmb_s1_p2 = QComboBox()
        self.le_s1_p2 = QLineEdit()
        self.lbl_s1_p2_units = QLabel()
        for w in [self.cmb_s1_p1, self.le_s1_p1, self.lbl_s1_p1_units,
                  self.cmb_s1_p2, self.le_s1_p2, self.lbl_s1_p2_units]:
            state1_layout.addWidget(w)
        self.state1_box.setLayout(state1_layout)
        spec_layout.addWidget(self.state1_box)

        # State 2
        self.state2_box = QGroupBox("State 2")
        state2_layout = QVBoxLayout()
        self.cmb_s2_p1 = QComboBox()
        self.le_s2_p1 = QLineEdit()
        self.lbl_s2_p1_units = QLabel()
        self.cmb_s2_p2 = QComboBox()
        self.le_s2_p2 = QLineEdit()
        self.lbl_s2_p2_units = QLabel()
        for w in [self.cmb_s2_p1, self.le_s2_p1, self.lbl_s2_p1_units,
                  self.cmb_s2_p2, self.le_s2_p2, self.lbl_s2_p2_units]:
            state2_layout.addWidget(w)
        self.state2_box.setLayout(state2_layout)
        spec_layout.addWidget(self.state2_box)

        spec_box.setLayout(spec_layout)
        main_layout.addWidget(spec_box)

        # Populate combo boxes
        properties = ["Pressure (p)", "Temperature (T)", "Internal Energy (u)",
                      "Enthalpy (h)", "Entropy (s)", "Specific Volume (v)", "Quality (x)"]
        for cmb in [self.cmb_s1_p1, self.cmb_s1_p2, self.cmb_s2_p1, self.cmb_s2_p2]:
            cmb.addItems(properties)

        # State Properties
        props_box = QGroupBox("State Properties")
        props_layout = QHBoxLayout()
        self.s1_props_box = QGroupBox("State 1")
        self.s2_props_box = QGroupBox("State 2")
        self.delta_props_box = QGroupBox("State Change")
        self.lbl_s1_props = QLabel()
        self.lbl_s2_props = QLabel()
        self.lbl_delta_props = QLabel()
        self.s1_props_box.setLayout(QVBoxLayout())
        self.s2_props_box.setLayout(QVBoxLayout())
        self.delta_props_box.setLayout(QVBoxLayout())
        self.s1_props_box.layout().addWidget(self.lbl_s1_props)
        self.s2_props_box.layout().addWidget(self.lbl_s2_props)
        self.delta_props_box.layout().addWidget(self.lbl_delta_props)
        props_layout.addWidget(self.s1_props_box)
        props_layout.addWidget(self.s2_props_box)
        props_layout.addWidget(self.delta_props_box)
        props_box.setLayout(props_layout)
        main_layout.addWidget(props_box)

        # Calculate Button and Warning
        self.btn_calculate = QPushButton("Calculate")
        self.lbl_warning = QLabel("")
        main_layout.addWidget(self.btn_calculate)
        main_layout.addWidget(self.lbl_warning)

        self.setLayout(main_layout)

        # Signals
        self.rdo_SI.clicked.connect(self.setUnits)
        self.rdo_English.clicked.connect(self.setUnits)
        for cmb in [self.cmb_s1_p1, self.cmb_s1_p2, self.cmb_s2_p1, self.cmb_s2_p2]:
            cmb.currentIndexChanged.connect(self.setUnits)
        self.btn_calculate.clicked.connect(self.calculateProperties)

    def setUnits(self):
        SI = self.rdo_SI.isChecked()
        newUnits = 'SI' if SI else 'EN'
        UnitChange = self.currentUnits != newUnits
        self.currentUnits = newUnits
        self.steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS if SI else XSteam.UNIT_SYSTEM_FLS)

        if SI:
            self.p_Units, self.t_Units = "bar", "C"
            self.u_Units, self.h_Units = "kJ/kg", "kJ/kg"
            self.s_Units, self.v_Units = "kJ/kg*C", "m^3/kg"
        else:
            self.p_Units, self.t_Units = "psi", "F"
            self.u_Units, self.h_Units = "btu/lb", "btu/lb"
            self.s_Units, self.v_Units = "btu/lb*F", "ft^3/lb"

        # Update State 1
        self.updateUnitsAndValues(self.cmb_s1_p1, self.le_s1_p1, self.lbl_s1_p1_units,
                                  self.cmb_s1_p2, self.le_s1_p2, self.lbl_s1_p2_units, SI, UnitChange)
        # Update State 2
        self.updateUnitsAndValues(self.cmb_s2_p1, self.le_s2_p1, self.lbl_s2_p1_units,
                                  self.cmb_s2_p2, self.le_s2_p2, self.lbl_s2_p2_units, SI, UnitChange)

    def updateUnitsAndValues(self, cmb1, le1, lbl1, cmb2, le2, lbl2, SI, UnitChange):
        props = [cmb1.currentText(), cmb2.currentText()]
        try:
            vals = [float(le1.text() or 0), float(le2.text() or 0)]
        except ValueError:
            vals = [0.0, 0.0]

        for i, (prop, val, lbl) in enumerate([(props[0], vals[0], lbl1), (props[1], vals[1], lbl2)]):
            if 'Pressure' in prop:
                lbl.setText(self.p_Units)
                if UnitChange:
                    vals[i] = val * UC.psi_to_bar if SI else val * UC.bar_to_psi
            elif 'Temperature' in prop:
                lbl.setText(self.t_Units)
                if UnitChange:
                    vals[i] = UC.F_to_C(val) if SI else UC.C_to_F(val)
            elif 'Energy' in prop:
                lbl.setText(self.u_Units)
                if UnitChange:
                    vals[i] = val * UC.btuperlb_to_kJperkg if SI else val * UC.kJperkg_to_btuperlb
            elif 'Enthalpy' in prop:
                lbl.setText(self.h_Units)
                if UnitChange:
                    vals[i] = val * UC.btuperlb_to_kJperkg if SI else val * UC.kJperkg_to_btuperlb
            elif 'Entropy' in prop:
                lbl.setText(self.s_Units)
                if UnitChange:
                    vals[i] = val * UC.btuperlbF_to_kJperkgC if SI else val * UC.kJperkgC_to_btuperlbF
            elif 'Volume' in prop:
                lbl.setText(self.v_Units)
                if UnitChange:
                    vals[i] = val * UC.ft3perlb_to_m3perkg if SI else val * UC.m3perkg_to_ft3perlb
            elif 'Quality' in prop:
                lbl.setText("")

        le1.setText("{:.3f}".format(vals[0]))
        le2.setText("{:.3f}".format(vals[1]))

    def makeLabel(self, state):
        return (f"Region = {state.region}\n"
                f"Pressure = {state.p:.3f} ({self.p_Units})\n"
                f"Temperature = {state.t:.3f} ({self.t_Units})\n"
                f"Internal Energy = {state.u:.3f} ({self.u_Units})\n"
                f"Enthalpy = {state.h:.3f} ({self.h_Units})\n"
                f"Entropy = {state.s:.3f} ({self.s_Units})\n"
                f"Specific Volume = {state.v:.3f} ({self.v_Units})\n"
                f"Quality = {state.x:.3f}")

    def makeDeltaLabel(self, state1, state2):
        delta = state2.__sub__(state1)
        return (f"Property Change:\n"
                f"ΔP = {delta.p:.3f} ({self.p_Units})\n"
                f"ΔT = {delta.t:.3f} ({self.t_Units})\n"
                f"Δu = {delta.u:.3f} ({self.u_Units})\n"
                f"Δh = {delta.h:.3f} ({self.h_Units})\n"
                f"Δs = {delta.s:.3f} ({self.s_Units})\n"
                f"Δv = {delta.v:.3f} ({self.v_Units})")

    def calculateProperties(self):
        self.state1 = thermoState()
        self.state2 = thermoState()
        SI = self.rdo_SI.isChecked()
        self.lbl_warning.setText("")

        # State 1
        SP1 = [self.cmb_s1_p1.currentText()[-2:-1].lower(), self.cmb_s1_p2.currentText()[-2:-1].lower()]
        if SP1[0] == SP1[1]:
            self.lbl_warning.setText("Warning: Same property specified twice for State 1.")
            return
        try:
            f1 = [float(self.le_s1_p1.text()), float(self.le_s1_p2.text())]
        except ValueError:
            self.lbl_warning.setText("Error: Invalid numeric input for State 1.")
            return
        self.state1.setState(SP1[0], SP1[1], f1[0], f1[1], SI)

        # State 2
        SP2 = [self.cmb_s2_p1.currentText()[-2:-1].lower(), self.cmb_s2_p2.currentText()[-2:-1].lower()]
        if SP2[0] == SP2[1]:
            self.lbl_warning.setText("Warning: Same property specified twice for State 2.")
            return
        try:
            f2 = [float(self.le_s2_p1.text()), float(self.le_s2_p2.text())]
        except ValueError:
            self.lbl_warning.setText("Error: Invalid numeric input for State 2.")
            return
        self.state2.setState(SP2[0], SP2[1], f2[0], f2[1], SI)

        # Display results
        self.lbl_s1_props.setText(self.makeLabel(self.state1))
        self.lbl_s2_props.setText(self.makeLabel(self.state2))
        self.lbl_delta_props.setText(self.makeDeltaLabel(self.state1, self.state2))


def main():
    app = QApplication(sys.argv)
    window = MainWindow()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()