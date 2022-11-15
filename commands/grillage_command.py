from PySide6.QtWidgets import QApplication, QMenu
from PySide6.QtWidgets import QDialog, QPushButton,QGridLayout,QToolTip,QCheckBox,QComboBox
from PySide6.QtWidgets import QTreeView,QMainWindow,QVBoxLayout,QHBoxLayout,QSizePolicy
from PySide6.QtWidgets import QTreeWidget,QTreeWidgetItem,QDockWidget,QWidget,QGroupBox
from PySide6.QtGui import QCursor


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
        actionNewHatch.triggered.connect(self.onNewHatchCover)

        actionGenerateFEM = self.menuAnalysis.addAction("&Generate FEM")
        actionGenerateFEM.triggered.connect(self.onGenerateFEM)

        actionRunTest = self.menuMain.addAction("&Generate Test Mesh")
        actionRunTest.triggered.connect(self.onActionGenerateTestMesh)

        actionRunTest = self.menuTests.addAction("&Calculate Mesh Dimensions")
        actionRunTest.triggered.connect(self.onActionCalculate_mesh_dimensions)

        actionRunTest = self.menuTests.addAction("&Show Beam Property")
        actionRunTest.triggered.connect(self.onActionTestProperty)

        actionRunTest = self.menuTests.addAction("&Full Model Overlap Check")
        actionRunTest.triggered.connect(self.onActionFullOverlapCheck)

        actionRunTest = self.menuTests.addAction("&Current TEST")
        actionRunTest.triggered.connect(self.onActionCurrentTest)

        actionRunTest = self.menuTests.addAction("&Test Mesh V2 base size")
        actionRunTest.triggered.connect(self.onActionTestV2basesize)

        actionRunTest = self.menuTests.addAction("&Test Mesh V2 transition mesh")
        actionRunTest.triggered.connect(self.onActionTestV2transition)

        try:
            manager.selected_geometry_changed.connect(self.onSelectedGeometryChanged)
            manager.geometry_created.connect(self.onGeometryCreated)
        except BaseException as error:
            print('An exception occurred: {}'.format(error))
        except:
            print('Unknown exception occurred during signals connection')
        self.mainwin.update()

    def onActionGenerateTestMesh(self):
        grill_fem = generate_test_mesh()
        grill_fem.regenerate()
        if grill_fem is not None:
            manager.add_geometry([grill_fem])
            manager.show_geometry([grill_fem])
        pass

    def onActionCalculate_mesh_dimensions(self):
        Test_calculate_mesh_dimensions()

    def onActionFullOverlapCheck(self):
        TestFullModelOverlap()

    def onActionTestProperty(self):
        grill_fem = generate_test_mesh()
        Test_GeoFEM_T_L_beam_property(grill_fem)

    def onActionCurrentTest(self):
        Test_current()

    def onActionTestV2basesize(self):
        Test_MeshVariant_V2_basesize()

    def onActionTestV2transition(self):
        Test_MeshVariant_V2_transition()

    def onGenerateFEM(self):
        QApplication.changeOverrideCursor(QCursor(Qt.WaitCursor))
        tart = timer()

        # mesh_extent = MeshExtent(self._grillgeo.grillage, AOS.NONE)    # Calculate mesh extents with Axis of Symmetry override
        mesh_extent = MeshExtent(self._grillgeo.grillage)                # Calculate mesh extents with automatic Axis of Symmetry discovery

        mesher = ElementSizeV1(mesh_extent)     # Calculate mesh dimensions for mesh variant V1
        # mesher = ElementSizeV2(mesh_extent)    # Calculate mesh dimensions for mesh variant V2

        # Mesh Control
        mesher.min_num_ebs = 1              # Minimum number of elements between stiffeners
        mesher.min_num_eweb = 3             # Minimum number of elements along psm web height
        mesher.num_eaf = 1                  # Number of elements across the psm flange
        mesher.flange_aspect_ratio = 8      # Max flange aspect ratio
        mesher.plate_aspect_ratio = 4       # Max plate aspect ratio
        mesher.des_plate_aspect_ratio = 3   # Desired plate aspect ratio

        grillage_mesh = GrillageMesh(mesher)
        grill_fem = grillage_mesh.generate_grillage_mesh_v1('Mesh variant V1 test')
        grill_fem.merge_coincident_nodes()

        grill_fem.regenerate()
        if grill_fem is not None:
            manager.add_geometry([grill_fem])
            manager.show_geometry([grill_fem])
        QApplication.restoreOverrideCursor()

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


def createCommand():
    return SGDCommand()