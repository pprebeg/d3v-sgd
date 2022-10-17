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
from grillage.grillage_model import Plate,PrimarySuppMem,Grillage,Segment, BeamDirection as bd
from grillage.grillage_model import TBeamProperty,BulbBeamProperty,HatBeamProperty
from grillage.grillage_model import GrillageModelData
from grillage.grillage_model import Grillage
from grillage.grillagemesher import GrillageGeometry

from core import geometry_manager as manager
import logging

class SGDCommand(Command):
    def __init__(self):
        super().__init__()
        self._app = QApplication.instance()
        importer=SGDImporter()
        self.app.registerIOHandler(importer)
        #self._tree: QTreeView = self.mainwin.window.findChild(QTreeView, "geometryTree")
        #self._tree.hide()
        self._grillgeo=None
        self.selected_geometries=[]
        self.menuMain = QMenu("Grillage")
        mb = self.mainwin.menuBar()
        mb.addMenu(self.menuMain)
        self.menuModel = QMenu("&Model")
        self.menuFEM = QMenu("&FEM")
        self.menuAnalysis = QMenu("&Analysis")
        self.menuMain.addMenu(self.menuModel)
        self.menuMain.addMenu(self.menuFEM)
        self.menuMain.addMenu(self.menuAnalysis)


        self.menuMain.show()

        actionNewHatch = self.menuModel.addAction("&New Hatch Cover")
        actionNewHatch.triggered.connect(self.onNewHatchCover)

        try:
            Signals.get().importGeometry.connect(self.register_grillage_geometry)
            manager.selected_geometry_changed.connect(self.onSelectedGeometryChanged)
            print('OK')
        except BaseException as error:
            print('An exception occurred: {}'.format(error))
        except:
            print('Unknown exception occurred during signals connection')


    def onNewHatchCover(self):
        grill= Grillage()
        self._grillgeo = GrillageGeometry(grill)

    @Slot()
    def register_grillage_geometry(self, grillgeo):
        if isinstance(grillgeo, GrillageGeometry):
            self._grillgeo = grillgeo

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