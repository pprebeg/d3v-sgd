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

def tmp_fun_gen_hc_var1():
    # Eksplicitno zadana topologija hc_var_1 za provjeru:  18.54m x 18.18m,  mreza nosaca 5 x 5
    hc_var_1 = Grillage(18.54, 18.18, 5, 5)

    # Lista materijala
    ST24 = MaterialProperty(1, 210000, 0.3, 7850, 235,
                            "ST24")  # PITANJE: Treba li ID biti dodjeljen automatski ili ga je ok ovako rucno zadavati?
    AH32 = MaterialProperty(2, 210000, 0.3, 7850, 315, "AH32")
    AH36 = MaterialProperty(3, 210000, 0.3, 7850, 355, "AH36")

    hc_var_1.add_material(ST24)
    hc_var_1.add_material(AH32)
    hc_var_1.add_material(AH36)

    # Korozijski dodatak
    tc = CorrosionAddition(1, 2)
    hc_var_1.add_corrosion_addition(tc)

    # hc_var_1.add_corrosion_addition(1, 2)     # ALTERNATIVNO ZADAVANJE I SPREMANJE

    # Beam property
    initial_longitudinal_beam = TBeamProperty(1, 1089, 10, 230, 16, ST24)  # inicijalni longitudinal T beam prop
    initial_transverse_beam = TBeamProperty(2, 1089, 10, 545, 40, ST24)  # inicijalni transverse T beam prop
    initial_edge_beam = LBeamProperty(3, 1089, 10, 150, 16, ST24)  # inicijalni rubni L beam prop
    # initial_edge_beam = FBBeamProperty(3, 1089, 10, ST24)               # inicijalni rubni FB beam prop
    initial_stiffener = HatBeamProperty(4, 220, 6, 220, 80, AH36)  # inicijalna ukrepa
    center_girder = TBeamProperty(5, 1089, 10, 560, 40, AH32)

    hc_var_1.add_beam_prop(
        initial_longitudinal_beam)  # PITANJE: Može li se ovo dodavanje preko add_property bolje izvesti?
    hc_var_1.add_beam_prop(initial_transverse_beam)
    hc_var_1.add_beam_prop(initial_edge_beam)
    hc_var_1.add_beam_prop(initial_stiffener)
    hc_var_1.add_beam_prop(center_girder)

    # Plate property
    plateprop1 = PlateProperty(1, 10, ST24)  # inicijalni plate property za cijeli poklopac
    plateprop2 = PlateProperty(2, 10, AH32)
    plateprop3 = PlateProperty(3, 9, ST24)
    plateprop4 = PlateProperty(4, 9, ST24)
    hc_var_1.add_plate_prop(plateprop1)
    hc_var_1.add_plate_prop(plateprop2)  # Dodatni plate property za testiranje identifikacije jednakih svojstava
    hc_var_1.add_plate_prop(plateprop3)
    hc_var_1.add_plate_prop(plateprop4)

    # Stiffener layouts
    stifflayout1 = StiffenerLayout(1, initial_stiffener, "spacing", 0.873)
    stifflayout2 = StiffenerLayout(2, initial_stiffener, "spacing", 0.935)
    hc_var_1.add_stiffener_layout(stifflayout1)  # dodavanje stiffener layouta u dictionary
    hc_var_1.add_stiffener_layout(stifflayout2)

    stiff_dir = BeamDirection.TRANSVERSE  # inicijalna orijentacija ukrepa na svim zonama oplate

    # Generacija topologije
    hc_var_1.generate_prim_supp_members()  # Generacija svih jakih nosaca
    hc_var_1.generate_segments(initial_longitudinal_beam, initial_transverse_beam,
                               initial_edge_beam)  # Generacija svih segmenata
    hc_var_1.generate_plating(plateprop1, stifflayout1, stiff_dir)  # Generacija oplate

    # Pridruzivanje simetricnih elemenata
    hc_var_1.assign_symmetric_members()
    hc_var_1.assign_symmetric_plating()
    hc_var_1.assign_symmetric_segments()

    # Izmjena polozaja jakih nosaca
    Grillage.set_all_longitudinal_PSM(hc_var_1, 4.5, 4.59, 4.59)
    Grillage.set_all_transverse_PSM(hc_var_1, 4.325, 4.935, 4.935)

    # Izmjene plating property - postavljanje drugačijih svojstava svim poljima oplate između poprečnih nosača koji definiraju polja oplate 1 i 4
    Grillage.set_plating_prop_transversals(hc_var_1, 1, "stiff_layout", stifflayout2)
    Grillage.set_plating_prop_transversals(hc_var_1, 4, "stiff_layout", stifflayout2)

    # Grillage.set_plating_prop_symmetric(hc_var_1, 6, "stiff_dir", BeamOrientation.LONGITUDINAL)
    # hc_var_1.plating()[5].stiff_layout = stifflayout1   # Unos drugacijeg layouta da hc_check javi gresku

    # Izmjena svojstava nosaca
    Grillage.set_tran_member_beam_property(hc_var_1, 3, center_girder)
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
        mb.addMenu(self.menuMain)

        actionNewHatch = self.menuModel.addAction("&New Hatch Cover")
        actionNewHatch.triggered.connect(self.onNewHatchCover)

        actionGenerateFEM = self.menuAnalysis.addAction("&Generate FEM")
        actionGenerateFEM.triggered.connect(self.onGenerateFEM)

        try:
            manager.selected_geometry_changed.connect(self.onSelectedGeometryChanged)
            manager.geometry_created.connect(self.onGeometryCreated)
        except BaseException as error:
            print('An exception occurred: {}'.format(error))
        except:
            print('Unknown exception occurred during signals connection')
        self.mainwin.update()

    def onGenerateFEM(self):
        QApplication.changeOverrideCursor(QCursor(Qt.WaitCursor))
        tart = timer()
        extents = MeshExtent(self._grillgeo.grillage, AOS.NONE)                # Opseg izrade mreže uz ručni odabir simetrije
        # extents = MeshExtent(self._grillgeo.grillage)  # Opseg izrade mreže uz automatsko prepoznavanje simetrije
        mesher = MeshV1(extents)  # Izračun dimenzija mreže za V1
        # mesher = MeshV2(extents)       # Izračun dimenzija mreže za V2

        # Kontrola mreže
        mesher.min_num_ebs = 1  # Postavljanje minimalnog broja elemenata između ukrepa
        mesher.min_num_eweb = 3  # Postavljanje minimalnog broja elemenata duž visine struka
        mesher.num_eaf = 1  # Postavljanje broja elemenata u smjeru širine prirubnice
        mesher.flange_aspect_ratio = 7  # Postavljanje aspektnog odnosa elemenata prirubnica jakih nosača i oplate uz struk jakih nosača
        mesher.plate_aspect_ratio = 4  # Postavljanje aspektnog odnosa elemenata oplate i strukova jakih nosača
        mesher.des_plate_aspect_ratio = 3  # Postavljanje poželjnog aspektnog odnosa elemenata oplate

        # Potrebno računati dimenzije mreže za sve testove osim generate_mesh_V1()
        mesher.calculate_mesh_dimensions()  # Izračun svih dimenzija za odabranu mrežu

        gm = GrillageMesh(mesher)
        grill_fem = gm.generate_mesh('naziv modela')
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