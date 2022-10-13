from PySide6.QtWidgets import QApplication, QMenu
from PySide6.QtWidgets import QDialog, QPushButton,QGridLayout,QToolTip,QCheckBox,QComboBox
from PySide6.QtWidgets import QTreeView,QMainWindow,QVBoxLayout,QHBoxLayout,QSizePolicy
from PySide6.QtWidgets import QTreeWidget,QTreeWidgetItem,QDockWidget,QWidget,QGroupBox


import os
from PySide6.QtCore import Slot,Qt,QPoint
#d3v imports
from commands import Command
from iohandlers import IOHandler
from signals import Signals
from typing import Dict,List, Tuple
from selinfo import SelectionInfo
from core import Geometry
from grillage.grillage_model import GrillageModelData,Grillage

class SGDCommand(Command):
    def __init__(self):
        super().__init__()
        self._app = QApplication.instance()
        importer=SGDImporter()
        self.app.registerIOHandler(importer)
        self._tree: QTreeView = self.mainwin.window.findChild(QTreeView, "geometryTree")
        #self._tree.hide()
        self._grillage=None
        self.si=0
        self.menuMain = QMenu("Grillage")
        self.menuModel = QMenu("&Model")
        self.menuFEM = QMenu("&FEM")
        self.menuAnalysis = QMenu("&Analysis")
        self.menuMain.addMenu(self.menuModel)
        self.menuMain.addMenu(self.menuFEM)
        self.menuMain.addMenu(self.menuAnalysis)

        actionNewHatch = self.menuModel.addAction("&New Hatch Cover")
        actionNewHatch.triggered.connect(self.onNewHatchCover)



        Signals.get().geometryImported.connect(self.registerGrillage)
        Signals.get().selectionChanged.connect(self.registerSelection)

    def onNewHatchCover(self):
        self._grillage = Grillage()
        pass

    @Slot()
    def registerGrillage(self, grillage):
        if isinstance(grillage, Geometry):
            self._grillage = grillage

    @Slot()
    def registerSelection(self, si):
        self.si = si
        if si.isEmpty():
            pass
        else:
            pass
            # currDBB = self.si.getGeometry()
            # print(self.dbbproblem)
            # if isinstance(currDBB, DBBBaseAll):
            #     pos: QPoint = self.mainwin.pos()
            #     pos.setX(pos.x() + self.mainwin.glWin.dragInfo.wStartPos.x() + 20)
            #     pos.setY(pos.y() + self.mainwin.glWin.size().height() - self.mainwin.glWin.dragInfo.wStartPos.y())
            #     msg = currDBB.get_info()
            #     QApplication.instance().clipboard().setText(str(msg))
            #     QToolTip.showText(pos, msg, msecShowTime=10)

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

    def importGeometry(self, fileName):
        if len(fileName) < 1:
            return
        filename_noext, file_extension = os.path.splitext(fileName)
        if file_extension not in self.getImportFormats():
            return
        grillage = GrillageModelData(fileName).read_file()
        if grillage != None:
            os.chdir(os.path.dirname(fileName))
            Signals.get().geometryImported.emit(grillage)

    def getImportFormats(self):
        return (".gin")


def createCommand():
    return SGDCommand()