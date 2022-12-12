import sys

import PySide6.QtCore
from PySide6.QtWidgets import QApplication, QMenu
from PySide6.QtWidgets import QDialog, QPushButton,QGridLayout,QToolTip,QCheckBox,QComboBox
from PySide6.QtWidgets import QTreeView,QMainWindow,QVBoxLayout,QHBoxLayout,QSizePolicy
from PySide6.QtWidgets import QTreeWidget,QTreeWidgetItem,QDockWidget,QWidget,QGroupBox
from PySide6.QtWidgets import QLabel, QSpinBox, QLineEdit, QDoubleSpinBox, QTableWidgetItem
from PySide6.QtGui import QCursor
from PySide6 import QtCore, QtGui, QtWidgets, QtUiTools

import os
from PySide6.QtCore import Slot,Qt,QPoint
#d3v imports
from commands import Command
from iohandlers import IOHandler
from signals import Signals
from typing import Dict,List, Tuple
from selinfo import SelectionInfo
from core import Geometry
from grillage.grillage_model import *
from grillage.grillage_visualizer import GrillageGeometry
from grillage.grillage_mesher import *
from timeit import default_timer as timer

from core import geometry_manager as manager
import logging
from tests.test_mesh import *

def tmp_fun_gen_hc_var1():
    # New hatch cover (L_overall [m], B_overall [m], N_longitudinal, N_transverse)
    hc_var_1 = Grillage(18.54, 18.18, 5, 5)

    # Materials list
    ST24 = MaterialProperty(1, 210000, 0.3, 7850, 235, "ST24")
    AH32 = MaterialProperty(2, 210000, 0.3, 7850, 315, "AH32")
    AH36 = MaterialProperty(3, 210000, 0.3, 7850, 355, "AH36")

    hc_var_1.add_material(ST24)
    hc_var_1.add_material(AH32)
    hc_var_1.add_material(AH36)

    tc = CorrosionAddition(1, 2)
    hc_var_1.add_corrosion_addition(tc)

    # Beam Properties
    initial_longitudinal_beam = TBeamProperty(1, 1089, 10, 230, 16, ST24)
    initial_transverse_beam = TBeamProperty(2, 1089, 10, 545, 40, ST24)
    initial_edge_beam = LBeamProperty(3, 1089, 10, 150, 16, ST24)
    initial_stiffener = HatBeamProperty(4, 150, 10, 200, 50, AH36)
    center_girder = TBeamProperty(5, 1089, 10, 560, 40, AH32)
    FB_beam = FBBeamProperty(6, 1089, 10, ST24)
    # initial_stiffener = BulbBeamProperty(4, 240, 10, AH36)
    # initial_stiffener = TBeamProperty(4, 250, 8, 90, 12, ST24)
    # initial_stiffener = FBBeamProperty(4, 250, 8, ST24)

    hc_var_1.add_beam_prop(initial_longitudinal_beam)
    hc_var_1.add_beam_prop(initial_transverse_beam)
    hc_var_1.add_beam_prop(initial_edge_beam)
    hc_var_1.add_beam_prop(initial_stiffener)
    hc_var_1.add_beam_prop(center_girder)
    hc_var_1.add_beam_prop(FB_beam)

    # Plate property
    plateprop1 = PlateProperty(1, 10, ST24)   # Initial plate property
    hc_var_1.add_plate_prop(plateprop1)

    # Stiffener layouts
    stifflayout1 = StiffenerLayout(1, initial_stiffener, DefinitionType.SPACING, 0.873)
    stifflayout2 = StiffenerLayout(2, initial_stiffener, DefinitionType.SPACING, 0.935)
    hc_var_1.add_stiffener_layout(stifflayout1)
    hc_var_1.add_stiffener_layout(stifflayout2)

    stiff_dir = BeamDirection.TRANSVERSE  # Initial beam orientation

    # Grillage generation
    hc_var_1.generate_prim_supp_members()
    hc_var_1.generate_segments(initial_longitudinal_beam, initial_transverse_beam, initial_edge_beam)
    hc_var_1.generate_plating(plateprop1, stifflayout1, stiff_dir)
    hc_var_1.generate_elementary_plate_panels()

    # Group setting of intercostal stiffeners
    # hc_var_1.plating()[6].set_intercostal_stiffeners(4, FB_beam)
    # hc_var_1.plating()[7].set_intercostal_stiffeners(4, FB_beam)
    # hc_var_1.plating()[10].set_intercostal_stiffeners(4, FB_beam)
    # hc_var_1.plating()[11].set_intercostal_stiffeners(4, FB_beam)

    # Individual placement of intercostal stiffeners
    # hc_var_1.plating()[1].elementary_plate_panels[1].intercostal_stiffener_num = 1
    # hc_var_1.plating()[1].elementary_plate_panels[1].beam_prop = FB_beam

    # Intercostal stiffener deletion from plating zone 6
    # hc_var_1.plating()[6].regenerate_elementary_plate_panel()

    hc_var_1.assign_symmetric_members()
    hc_var_1.assign_symmetric_plating()
    hc_var_1.assign_symmetric_segments()

    # Set all Primary Supporting Member spacing, [m]
    hc_var_1.set_all_longitudinal_PSM(4.5, 4.59, 4.59)
    hc_var_1.set_all_transverse_PSM(4.325, 4.935, 4.935)

    # Set different stiffener layout for selected plate and all other plating
    # zones between Primary Supporting Members defined by that selected plate
    hc_var_1.set_plating_prop_transversals(1, "stiff_layout", stifflayout2)
    hc_var_1.set_plating_prop_transversals(4, "stiff_layout", stifflayout2)

    # Simultaneous change of all symmetric plating properties, possible selection
    # hc_var_1.set_plating_prop_symmetric(1, "stiff_dir", BeamDirection.LONGITUDINAL)

    # Set same Beam Property for the entire transverse Primary Supporting Member
    # hc_var_1.set_tran_member_beam_property(1, FB_beam)
    hc_var_1.set_tran_member_beam_property(3, center_girder)

    return hc_var_1


# noinspection PyUnresolvedReferences
class GenerateNewHC(QDialog):
    def __init__(self, parent):
        super().__init__(parent)
        loader = QtUiTools.QUiLoader()
        ui_file = QtCore.QFile("../commands/hc_GUI.ui")  # Qt designer ui file
        ui_file.open(QtCore.QFile.ReadOnly)
        self.hc_gui = loader.load(ui_file, parent)
        ui_file.close()
        self.hc_gui.setWindowTitle("Generate New Grillage Structure")
        self.hc_gui.show()

        self.table_materials = self.hc_gui.tableMaterials
        self.table_plate = self.hc_gui.tablePlate
        self.table_beams = self.hc_gui.tableBeamProps
        self.table_layouts = self.hc_gui.tableStiffenerLayouts

        self.cbox_plate_mat = self.hc_gui.cBox_PlateMaterial
        self.cbox_beam_mat = self.hc_gui.cBox_BeamMaterial
        self.cbox_beam_type = self.hc_gui.cBox_BeamType
        self.cbox_layout_beam = self.hc_gui.cBox_LayoutBeamProp
        self.cbox_long_beam = self.hc_gui.cBox_LongitudinalBeamProp
        self.cbox_tran_beam = self.hc_gui.cBox_TransverseBeamProp
        self.cbox_stiff_beam = self.hc_gui.cBox_StiffenerBeamProp
        self.cbox_layout_type = self.hc_gui.cBox_LayoutDefinitionType

        self.cbox_beam_type.currentIndexChanged.connect(self.update_beam_combobox)

        self.hc_gui.btnAddNewMaterial.clicked.connect(self.add_new_material_prop)
        self.hc_gui.btnAddNewPlate.clicked.connect(self.add_new_plate_prop)
        self.hc_gui.btnAddNewBeam.clicked.connect(self.add_new_beam)
        self.hc_gui.btnAddNewLayout.clicked.connect(self.add_new_layout)
        self.hc_gui.btnGenerateGrillage.clicked.connect(self.gnerate_new_hc)

        self.load_materials_list(self.default_materials())
        self.update_material_combobox()
        self.update_beam_combobox()

    def gnerate_new_hc(self):
        # row_count = self.table_materials.rowCount()
        # for row in range(0, row_count):
        #     mat_id = int(self.table_materials.item(row, 0).text())
        #     name = self.table_materials.item(row, 1).text()
        #     E = float(self.table_materials.item(row, 2).text())
        #     v = float(self.table_materials.item(row, 3).text())
        #     Reh = float(self.table_materials.item(row, 4).text())
        #     ro = float(self.table_materials.item(row, 5).text())
        #
        # return hc_variant
        pass

    @staticmethod
    def default_materials():
        def_mats = [{"ID": 1, "E": 210000, "v": 0.3, "ro": 7.85 * 10**(-9), "Reh": 235, "name": "ST24"},
                    {"ID": 1, "E": 210000, "v": 0.3, "ro": 7.85 * 10**(-9), "Reh": 315, "name": "AH32"},
                    {"ID": 3, "E": 210000, "v": 0.3, "ro": 7.85 * 10**(-9), "Reh": 355, "name": "AH36"}]
        return def_mats

    def load_materials_list(self, default_material_list):
        """
        Load default materials list into table_materials.
        """
        self.table_materials.setRowCount(3)
        row_index = 0
        for material in default_material_list:
            ID = QTableWidgetItem(str(material["ID"]))
            name = QTableWidgetItem(material["name"])
            E = QTableWidgetItem(str(material["E"]))
            v = QTableWidgetItem(str(material["v"]))
            Reh = QTableWidgetItem(str(material["Reh"]))
            ro = QTableWidgetItem(str(material["ro"]))

            ID.setTextAlignment(Qt.AlignCenter)
            name.setTextAlignment(Qt.AlignCenter)
            E.setTextAlignment(Qt.AlignCenter)
            v.setTextAlignment(Qt.AlignCenter)
            Reh.setTextAlignment(Qt.AlignCenter)
            ro.setTextAlignment(Qt.AlignCenter)

            self.table_materials.setItem(row_index, 0, ID)
            self.table_materials.setItem(row_index, 1, name)
            self.table_materials.setItem(row_index, 2, E)
            self.table_materials.setItem(row_index, 3, v)
            self.table_materials.setItem(row_index, 4, Reh)
            self.table_materials.setItem(row_index, 5, ro)
            row_index += 1

    def add_new_material_prop(self):
        """
        Add new material property to table_materials.
        """
        mat_ID_input = self.hc_gui.lineEdit_materialID.text()
        name_input = self.hc_gui.lineEdit_materialName.text()
        E_input = self.hc_gui.lineEdit_MaterialE.text()
        v_input = self.hc_gui.lineEdit_Materialv.text()
        Reh_input = self.hc_gui.lineEdit_MaterialReh.text()
        ro_input = self.hc_gui.lineEdit_MaterialRo.text()

        ID = QTableWidgetItem(mat_ID_input)
        name = QTableWidgetItem(name_input)
        E = QTableWidgetItem(E_input)
        v = QTableWidgetItem(v_input)
        Reh = QTableWidgetItem(Reh_input)
        ro = QTableWidgetItem(ro_input)

        ID.setTextAlignment(Qt.AlignCenter)
        name.setTextAlignment(Qt.AlignCenter)
        E.setTextAlignment(Qt.AlignCenter)
        v.setTextAlignment(Qt.AlignCenter)
        Reh.setTextAlignment(Qt.AlignCenter)
        ro.setTextAlignment(Qt.AlignCenter)

        row_count = self.table_materials.rowCount()
        self.table_materials.insertRow(row_count)
        self.table_materials.setItem(row_count, 0, ID)
        self.table_materials.setItem(row_count, 1, name)
        self.table_materials.setItem(row_count, 2, E)
        self.table_materials.setItem(row_count, 3, v)
        self.table_materials.setItem(row_count, 4, Reh)
        self.table_materials.setItem(row_count, 5, ro)
        self.update_material_combobox()

    def get_materials_list(self):
        """
        :return: List of material names entered into table_materials.
        """
        row_count = self.table_materials.rowCount()
        mat_list = [None] * row_count
        for row in range(0, row_count):
            name = self.table_materials.item(row, 1).text()
            mat_list[row] = name
        return mat_list

    def update_material_combobox(self):
        """
        Updates material comboboxes for new plate and beam generation.
        """
        self.cbox_plate_mat.clear()
        self.cbox_beam_mat.clear()
        self.cbox_plate_mat.addItems(self.get_materials_list())
        self.cbox_beam_mat.addItems(self.get_materials_list())

    def update_stiffener_beam_combobox(self):
        """
        Updates stiffener beam combobox for new layout and grillage generation.
        """
        self.cbox_layout_beam.clear()
        self.cbox_stiff_beam.clear()
        self.cbox_layout_beam.addItems(self.get_beam_list())
        self.cbox_stiff_beam.addItems(self.get_beam_list())

    def add_new_plate_prop(self):
        """
        Add new plate property to table_plate.
        """
        plate_id_input = self.hc_gui.lineEdit_plateID.text()
        plate_name_input = self.hc_gui.lineEdit_plateName.text()
        tp_input = self.hc_gui.lineEdit_plateThickness.text()
        plate_mat_input = self.cbox_plate_mat.currentText()

        plate_id = QTableWidgetItem(plate_id_input)
        plate_name = QTableWidgetItem(plate_name_input)
        tp = QTableWidgetItem(tp_input)
        plate_mat = QTableWidgetItem(plate_mat_input)

        plate_id.setTextAlignment(Qt.AlignCenter)
        plate_name.setTextAlignment(Qt.AlignCenter)
        tp.setTextAlignment(Qt.AlignCenter)
        plate_mat.setTextAlignment(Qt.AlignCenter)

        row_count = self.table_plate.rowCount()
        self.table_plate.insertRow(row_count)
        self.table_plate.setItem(row_count, 0, plate_id)
        self.table_plate.setItem(row_count, 1, plate_name)
        self.table_plate.setItem(row_count, 2, tp)
        self.table_plate.setItem(row_count, 3, plate_mat)

    def get_T_beam_dims(self):
        """
        :return: T beam dimensions string.
        """
        hw = self.hc_gui.lineEdit_hw_T.text()
        tw = self.hc_gui.lineEdit_tw_T.text()
        bf = self.hc_gui.lineEdit_bf_T.text()
        tf = self.hc_gui.lineEdit_tf_T.text()
        dim_str = hw + "x" + tw + "/"
        dim_str += bf + "x" + tf
        return dim_str

    def get_L_beam_dims(self):
        """
        :return: L beam dimensions string.
        """
        hw = self.hc_gui.lineEdit_hw_L.text()
        tw = self.hc_gui.lineEdit_tw_L.text()
        bf = self.hc_gui.lineEdit_bf_L.text()
        tf = self.hc_gui.lineEdit_tf_L.text()
        dim_str = hw + "x" + tw + "/"
        dim_str += bf + "x" + tf
        return dim_str

    def get_FB_beam_dims(self):
        """
        :return: FB beam dimensions string.
        """
        hw = self.hc_gui.lineEdit_hw_FB.text()
        tw = self.hc_gui.lineEdit_tw_FB.text()
        dim_str = hw + "x" + tw
        return dim_str

    def get_Bulb_beam_dims(self):
        """
        :return: Bulb beam dimensions string.
        """
        hw = self.hc_gui.lineEdit_hw_Bulb.text()
        tw = self.hc_gui.lineEdit_tw_Bulb.text()
        dim_str = hw + "x" + tw
        return dim_str

    def get_Hat_beam_dims(self):
        """
        :return: L beam dimensions string.
        """
        h = self.hc_gui.lineEdit_h_Hat.text()
        t = self.hc_gui.lineEdit_t_Hat.text()
        bf = self.hc_gui.lineEdit_bf_Hat.text()
        fi = self.hc_gui.lineEdit_fi_Hat.text()
        dim_str = h + "x" + t + "x"
        dim_str += bf + "x" + fi
        return dim_str

    def get_beam_dimensions(self):
        beam_type = self.cbox_beam_type.currentText()
        dim_str = None
        if beam_type == "T":
            dim_str = self.get_T_beam_dims()
        elif beam_type == "L":
            dim_str = self.get_L_beam_dims()
        elif beam_type == "FB":
            dim_str = self.get_FB_beam_dims()
        elif beam_type == "Bulb":
            dim_str = self.get_Bulb_beam_dims()
        elif beam_type == "Hat":
            dim_str = self.get_Hat_beam_dims()
        return dim_str

    def add_new_beam(self):
        """
        Add new beam property to table_beams.
        """
        beam_id_input = self.hc_gui.lineEdit_beamID.text()
        beam_name_input = self.hc_gui.lineEdit_beamName.text()
        beam_type_input = self.cbox_beam_type.currentText()
        beam_dims_input = self.get_beam_dimensions()
        beam_mat_input = self.cbox_beam_mat.currentText()

        beam_id = QTableWidgetItem(beam_id_input)
        beam_name = QTableWidgetItem(beam_name_input)
        beam_type = QTableWidgetItem(beam_type_input)
        beam_dims = QTableWidgetItem(beam_dims_input)
        beam_mat = QTableWidgetItem(beam_mat_input)

        beam_id.setTextAlignment(Qt.AlignCenter)
        beam_name.setTextAlignment(Qt.AlignCenter)
        beam_type.setTextAlignment(Qt.AlignCenter)
        beam_dims.setTextAlignment(Qt.AlignCenter)
        beam_mat.setTextAlignment(Qt.AlignCenter)

        row_count = self.table_beams.rowCount()
        self.table_beams.insertRow(row_count)
        self.table_beams.setItem(row_count, 0, beam_id)
        self.table_beams.setItem(row_count, 1, beam_name)
        self.table_beams.setItem(row_count, 2, beam_type)
        self.table_beams.setItem(row_count, 3, beam_dims)
        self.table_beams.setItem(row_count, 4, beam_mat)
        self.update_stiffener_beam_combobox()

    def get_beam_list(self):
        row_count = self.table_beams.rowCount()
        beam_list = [None] * row_count
        for row in range(0, row_count):
            name = self.table_beams.item(row, 1).text()
            beam_list[row] = name
        return beam_list

    def add_new_layout(self):
        """
        Add new stiffener layout to table_layouts.
        """
        layout_id_input = self.hc_gui.lineEdit_LayoutID.text()
        layout_name_input = self.hc_gui.lineEdit_LayoutName.text()
        layout_beam_input = self.cbox_layout_beam.currentText()
        layout_type_input = self.cbox_layout_type.currentText()
        layout_value_input = self.hc_gui.lineEdit_LayoutValue.text()

        layout_id = QTableWidgetItem(layout_id_input)
        layout_name = QTableWidgetItem(layout_name_input)
        layout_beam = QTableWidgetItem(layout_beam_input)
        layout_type = QTableWidgetItem(layout_type_input)
        layout_value = QTableWidgetItem(layout_value_input)

        layout_id.setTextAlignment(Qt.AlignCenter)
        layout_name.setTextAlignment(Qt.AlignCenter)
        layout_beam.setTextAlignment(Qt.AlignCenter)
        layout_type.setTextAlignment(Qt.AlignCenter)
        layout_value.setTextAlignment(Qt.AlignCenter)

        row_count = self.table_layouts.rowCount()
        self.table_layouts.insertRow(row_count)
        self.table_layouts.setItem(row_count, 0, layout_id)
        self.table_layouts.setItem(row_count, 1, layout_name)
        self.table_layouts.setItem(row_count, 2, layout_beam)
        self.table_layouts.setItem(row_count, 3, layout_type)
        self.table_layouts.setItem(row_count, 4, layout_value)

    def update_beam_combobox(self):
        """
        Hides all overlapping group boxes (gb) and shows the selected one.
        """
        beam_type = self.cbox_beam_type.currentText()
        gb_T = self.hc_gui.groupBox_T_Scantlings
        gb_L = self.hc_gui.groupBox_L_Scantlings
        gb_FB = self.hc_gui.groupBox_FB_Scantlings
        gb_Bulb = self.hc_gui.groupBox_Bulb_Scantlings
        gb_Hat = self.hc_gui.groupBox_Hat_Scantlings

        gb_T.hide()
        gb_L.hide()
        gb_FB.hide()
        gb_Bulb.hide()
        gb_Hat.hide()

        if beam_type == "T":
            gb_T.show()
        elif beam_type == "L":
            gb_L.show()
        elif beam_type == "FB":
            gb_FB.show()
        elif beam_type == "Bulb":
            gb_Bulb.show()
        elif beam_type == "Hat":
            gb_Hat.show()


class SGDCommand(Command):
    def __init__(self):
        super().__init__()
        self._app = QApplication.instance()
        importer=SGDImporter()
        self.app.registerIOHandler(importer)
        #self._tree: QTreeView = self.mainwin.window.findChild(QTreeView, "geometryTree")
        #self._tree.hide()
        self._grillgeo:GrillageGeometry=None
        self.selected_geometries=[]
        self.menuMain = QMenu("Grillage")
        mb = self.mainwin.menuBar()

        self.menuModel = QMenu("&Model")
        self.menuMain.addMenu(self.menuModel)
        self.menuAnalysis = QMenu("&Analysis")
        self.menuMain.addMenu(self.menuAnalysis)
        self.menuTests = QMenu("&Tests")
        self.menuMain.addMenu(self.menuTests)
        mb.addMenu(self.menuMain)

        actionNewHatch = self.menuModel.addAction("&New Hatch Cover")
        actionNewHatch.triggered.connect(self.onNewHatchCoverGUI)

        actionGenerateFEM = self.menuAnalysis.addAction("&Generate FEM...")
        actionGenerateFEM.triggered.connect(self.onGenerateFEM)

        actionRunTest = self.menuTests.addAction("&Generate Test Mesh V1")
        actionRunTest.triggered.connect(self.onActionGenerateTestMeshV1)
        
        actionRunTest = self.menuTests.addAction("&Generate Test Mesh V2")
        actionRunTest.triggered.connect(self.onActionGenerateTestMeshV2)

        actionRunTest = self.menuTests.addAction("&Test Edge Node Spacing")
        actionRunTest.triggered.connect(Test_edge_node_spacing)

        try:
            manager.selected_geometry_changed.connect(self.onSelectedGeometryChanged)
            manager.geometry_created.connect(self.onGeometryCreated)
        except BaseException as error:
            print('An exception occurred: {}'.format(error))
        except:
            print('Unknown exception occurred during signals connection')
        self.mainwin.update()
        self.mainwin.resize(1000, 700)

    def onActionGenerateTestMeshV1(self):
        grill_fem = generate_test_mesh_v1()
        grill_fem.regenerate()
        if grill_fem is not None:
            manager.add_geometry([grill_fem])
            manager.show_geometry([grill_fem])
        pass

    def onActionGenerateTestMeshV2(self):
        grill_fem = generate_test_mesh_v2()
        grill_fem.regenerate()
        if grill_fem is not None:
            manager.add_geometry([grill_fem])
            manager.show_geometry([grill_fem])
        pass

    def onGenerateFEM(self):
        analysis_gui = GrillageAnalysisGUI(self.mainwin, self._grillgeo)
        analysis_gui.mesh_parameters_gui()
        analysis_gui.exec()
        pass

    def onNewHatchCover(self):
        QApplication.changeOverrideCursor(QCursor(Qt.WaitCursor))
        old_grillgeo = self._grillgeo
        grill= tmp_fun_gen_hc_var1()
        self._grillgeo = GrillageGeometry(grill,'New hatch cover name')
        if old_grillgeo is not None:
            manager.remove_geometry([old_grillgeo])
        if self._grillgeo is not None:
            manager.add_geometry([self._grillgeo])
            manager.show_geometry([self._grillgeo])
        QApplication.restoreOverrideCursor()

    def onNewHatchCoverGUI(self):
        new_hc = GenerateNewHC(self.mainwin)
        return new_hc

    @Slot()
    def onGeometryCreated(self, geometries: List[Geometry]):
        for g in geometries:
            if isinstance(g, GrillageGeometry):
                self._grillgeo = g
                break

    @Slot()
    def onSelectedGeometryChanged(self, visible: List[Geometry], loaded: List[Geometry], selected: List[Geometry]):
        self.selected_geometries=selected

    @property
    def app(self):
        return self._app

    @property
    def mainwin(self):
        return self.app.mainFrame

    @property
    def glwin(self):
        return self.mainwin.glWin


class SGDImporter(IOHandler):
    def __init__(self):
        super().__init__()

    def do_import_geometry(self, fileName):
        if len(fileName) < 1:
            return
        filename_noext, file_extension = os.path.splitext(fileName)
        if file_extension not in self.getImportFormats():
            return
        QApplication.changeOverrideCursor(QCursor(Qt.WaitCursor))
        grillage = GrillageModelData(fileName).read_file()
        if grillage != None:
            os.chdir(os.path.dirname(fileName))
            g=GrillageGeometry(grillage,filename_noext)
            logging.debug("do_import_geometry: {}".format(g.guid))
            QApplication.restoreOverrideCursor()
            return g
        else:
            QApplication.restoreOverrideCursor()

    def getImportFormats(self):
        return (".gin")


class GrillageAnalysisGUI(QDialog):
    def __init__(self, parent, grillgeo: GrillageGeometry):
        super().__init__(parent)
        self.mainwin = parent
        self._grillgeo = grillgeo
        self.setWindowTitle("Grillage Finite Element Mesh Generation")
        self.setFixedSize(450, 310)
        self.main_layout = QGridLayout()
        self.setLayout(self.main_layout)

        self.mesh_variant = QComboBox(self)
        self.aos_override = QComboBox(self)
        self.name_input = QLineEdit(self)
        self.ebs_input = QSpinBox(self)
        self.eweb_input = QSpinBox(self)
        self.eaf_input = QSpinBox(self)
        self.far_input = QDoubleSpinBox(self)
        self.par_input = QDoubleSpinBox(self)
        self.dpar_input = QDoubleSpinBox(self)

    def grillage_mesh_variant(self):
        grill_var = self._grillgeo.grillage
        aos_override = self.get_aos()
        mesh_var = self.mesh_variant.currentText()
        if mesh_var == "V1":
            return MeshVariantV1(grill_var, aos_override)
        elif mesh_var == "V2":
            return MeshVariantV2(grill_var, aos_override)

    def get_aos(self):
        aos = self.aos_override.currentText()
        if aos == "Automatic Axis of Symmetry discovery":
            return None
        elif aos == "Longitudinal Axis of Symmetry override":
            return AOS.LONGITUDINAL
        elif aos == "Transverse Axis of Symmetry override":
            return AOS.TRANSVERSE
        elif aos == "No Axis of Symmetry override":
            return AOS.NONE

    def generate_grill_mesh(self):
        QApplication.changeOverrideCursor(QCursor(Qt.WaitCursor))
        name = self.name_input.text()
        ebs = self.ebs_input.value()
        eweb = self.eweb_input.value()
        eaf = self.eaf_input.value()
        far = self.far_input.value()
        par = self.par_input.value()
        dpar = self.dpar_input.value()
        self.close()

        mesh = self.grillage_mesh_variant()
        symmetry = mesh.mesh_extent.axis_of_symm
        grill_fem = mesh.generate_grillage_mesh(name, ebs, eweb, eaf, far, par, dpar)

        pressure = 0.0343    # N/mm2
        mesh.generate_loadcase(grill_fem, symmetry, pressure)

        gravity = -9810.0    # mm/s2
        weight_val = [0, 0, gravity, 0, 0, 0]
        mesh.generate_self_weight(grill_fem, weight_val)

        grill_fem.regenerate()
        if grill_fem is not None:
            manager.add_geometry([grill_fem])
            manager.show_geometry([grill_fem])
        QApplication.restoreOverrideCursor()
        return grill_fem

    def mesh_var_label(self):
        label = QLabel(self)
        label.setText("Mesh variant")
        self.main_layout.addWidget(label, 0, 0)

    def mesh_var_widget(self):
        self.mesh_variant.addItem("V1")
        self.mesh_variant.addItem("V2")
        self.mesh_variant.setMaximumWidth(70)
        self.main_layout.addWidget(self.mesh_variant, 0, 1)

    def aos_override_widget(self):
        self.aos_override.addItem("Automatic Axis of Symmetry discovery")
        self.aos_override.addItem("Longitudinal Axis of Symmetry override")
        self.aos_override.addItem("Transverse Axis of Symmetry override")
        self.aos_override.addItem("No Axis of Symmetry override")
        self.main_layout.addWidget(self.aos_override, 1, 0)

    def name_widget(self):
        self.name_input.setText("Grillage mesh")
        self.main_layout.addWidget(self.name_input, 2, 0)

    def user_input_labels(self):
        input_labels = {3: "Number of elements between stiffeners",
                        4: "Number of elements along the height of PSM web",
                        5: "Number of elements across PSM flange",
                        6: "Max PSM flange aspect ratio",
                        7: "Max plating and PSM web aspect ratio",
                        8: "Desired plating aspect ratio"}

        for key, val in input_labels.items():
            label = QLabel(self)
            label.setText(val)
            self.main_layout.addWidget(label, key, 0)

    def ar_input_widgets(self):
        self.far_input.setMinimum(1)
        self.far_input.setValue(7)
        self.far_input.setDecimals(1)
        self.far_input.setSingleStep(0.1)
        self.main_layout.addWidget(self.far_input, 6, 1)

        self.par_input.setMinimum(1)
        self.par_input.setValue(3)
        self.par_input.setDecimals(1)
        self.par_input.setSingleStep(0.1)
        self.main_layout.addWidget(self.par_input, 7, 1)

        self.dpar_input.setMinimum(1)
        self.dpar_input.setValue(2)
        self.dpar_input.setDecimals(1)
        self.dpar_input.setSingleStep(0.1)
        self.main_layout.addWidget(self.dpar_input, 8, 1)

    def element_num_widgets(self):
        self.ebs_input.setValue(1)
        self.main_layout.addWidget(self.ebs_input, 3, 1)

        self.eweb_input.setMinimum(1)
        self.eweb_input.setValue(3)
        self.main_layout.addWidget(self.eweb_input, 4, 1)

        self.eaf_input.setMinimum(1)
        self.eaf_input.setValue(1)
        self.main_layout.addWidget(self.eaf_input, 5, 1)

    def button_generate(self):
        generate_button = QPushButton("Generate mesh", self)
        self.main_layout.addWidget(generate_button, 9, 3)
        generate_button.clicked.connect(self.generate_grill_mesh)

    def mesh_parameters_gui(self):
        self.mesh_var_label()
        self.mesh_var_widget()
        self.aos_override_widget()
        self.name_widget()
        self.user_input_labels()
        self.ar_input_widgets()
        self.element_num_widgets()
        self.button_generate()


def createCommand():
    return SGDCommand()