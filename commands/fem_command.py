try:
    from PySide6.QtWidgets import QApplication, QMenu
    from PySide6.QtWidgets import QDialog, QPushButton,QGridLayout,QToolTip,QCheckBox,QComboBox
    from PySide6.QtWidgets import QTreeView,QMainWindow,QVBoxLayout,QHBoxLayout,QSizePolicy
    from PySide6.QtWidgets import QTreeWidget,QTreeWidgetItem,QDockWidget,QWidget,QGroupBox
    from PySide6.QtWidgets import QInputDialog,QLineEdit,QFormLayout,QTextEdit
    from PySide6.QtGui import QCursor
    import importlib

    import os
    from PySide6.QtCore import Slot,Qt,QPoint,QObject
    #d3v imports
    from commands import Command
    from iohandlers import IOHandler
    from signals import Signals
    from typing import Dict,List, Tuple
    from selinfo import SelectionInfo
    from core import Geometry

    from femdir.oofemwrap import *

    from femdir.geooofem import GeoOOFEM
    from femdir.geofem import GeoFEM
    from core import geometry_manager as manager
    import logging
except BaseException as error:
    print('An exception occurred: {}'.format(error))
except:
    print('Unknown exception occurred')

class FEMCommand(Command):
    def __init__(self):
        null = None
        super().__init__()
        self._app = QApplication.instance()
        importer=FEMImporter()
        self.app.registerIOHandler(importer)

        self.femmdl: GeoFEM = null
        self.si = 0
        self._curlc = null
        self._combo_lc: QComboBox = null
        self._menuMain = null
        self._show_oofem_analysis_menu = False
        # path to oofem
        self._oofem_input_filepath = ''
        self._oofem_idset_outtypes = null
        oofempymodfind = importlib.util.find_spec("oofempy") is not None
        if oofempymodfind:
            self._show_oofem_analysis_menu = True

        self.add_toolbars_and_menus()


        try:
            #manager.selected_geometry_changed.connect(self.onSelectedGeometryChanged)
            self.glwin.selector.selection_info_changled.connect(self.registerSelection)
            manager.geometry_created.connect(self.onGeometryCreated)
        except BaseException as error:
            print('An exception occurred: {}'.format(error))
        except:
            print('Unknown exception occurred during signals connection')
        self.mainwin.update()

    def add_toolbars_and_menus(self):
        mb = self.mainwin.menuBar()

        self._menuMain = QMenu("FEM")
        self._menuViewType = QMenu("&View Type")
        self._menuMain.addMenu(self._menuViewType)
        self._menuResults = QMenu("&Results")
        self._menuMain.addMenu(self._menuResults)
        self._menuOOFEM = QMenu("&OOFEM")
        self._menuMain.addMenu(self._menuOOFEM)

        menuExportOOFEM = self._menuOOFEM.addAction("&Export input file")
        menuExportOOFEM.triggered.connect(self.onExportOOFEMin)

        menuMeshControl = self._menuMain.addAction("View Control")
        menuMeshControl.triggered.connect(self.onMeshControl)

        mb.addMenu(self._menuMain)
        self._combo_lc = self.init_combo_lc()
        if self._show_oofem_analysis_menu:
            self.add_OOFEM_analysis_toolbars_and_menus()

    def add_OOFEM_analysis_toolbars_and_menus(self):
        menu_anylyse_oofem = self._menuOOFEM.addAction("&Analyse")
        menu_anylyse_oofem.triggered.connect(self.on_analyse_oofem)

    def init_combo_lc(self)->QComboBox:
        combo = QComboBox()
        if self.femmdl == None or self.femmdl.loadcases ==None:
            return combo
        for lc in self.femmdl.loadcases.values():
            combo.addItem(str(lc.id)+': '+ lc.name,lc)
        combo.currentIndexChanged.connect(self.loadcasechanged)
        return combo
    def loadcasechanged(self,index:int):

        text = self._combo_lc.itemText(index)
        self._curlc = self._combo_lc.itemData(index, Qt.UserRole)
        pass


    @Slot()
    def onGeometryCreated(self, geometries:List[Geometry]):
        for g in geometries:
            if isinstance(g, GeoFEM):
                self.femmdl = g
                self.addViewTypeAndResulsMenus(self.femmdl,self._menuViewType, self._menuResults)
                break
            if isinstance(g, GeoOOFEM):
                self._oofem_input_filepath = g._model.filename

    def addViewTypeAndResulsMenus(self,femmdl:GeoFEM, menuViewType:QMenu, menuResults:QMenu):
        menuViewType.clear()
        menuResults.clear()
        for key, atr in femmdl.attrib_val_functions.items():
            menuResul = menuViewType.addAction(key)
            menuResul.triggered.connect(self.onColorControlMenu)
        for key, res in femmdl.element_results.items():
            menuResul = menuResults.addAction(key)
            menuResul.triggered.connect(self.onColorControlMenu)


    def onColorControlMenu(self):
        action = QObject.sender(self)
        key = action.iconText()
        self.femmdl.prepareModelForVisualization(key)
        self.meshctrl.setCurrentGeoFEM(self.femmdl)
        self.femmdl.emit_geometries_rebuild()
        pass

    @Slot()
    def registerSelection(self, si: SelectionInfo):
        self.si = si
        if si.isEmpty():
            self.femmdl.unselect()
        else:
            pos: QPoint = self._mainwin.pos()
            pos.setX(pos.x() + self._mainwin.glWin.dragInfo.wStartPos.x() + 20)
            pos.setY(pos.y() + self._mainwin.glWin.size().height() - self._mainwin.glWin.dragInfo.wStartPos.y())
            msg = self.femmdl.selected_entitiy.get_info()
            QApplication.instance().clipboard().setText(str(msg))
            QToolTip.showText(pos, msg, msecShowTime=10)

    def onOptimize(self):
        if isinstance(self.femmdl, GeoFEM):
            pass

    def onMeshControl(self):
        if isinstance(self.femmdl, GeoFEM):
            self.meshctrl.setCurrentGeoFEM(self.femmdl)
            self.meshctrl.exec()

    def onExportOOFEMin(self):
        if isinstance(self.femmdl, GeoFEM):
            eltypes = get_initial_element_types_dict()
            dir_path = os.path.dirname(self.femmdl.get_input_file_path())
            text, ok = QInputDialog.getText(self._mainwin, 'GeoOOFEM Export input file', 'Use next file name?',
                                            QLineEdit.Normal, 'oofem_exp_input.in')
            if ok:
                self._oofem_input_filepath = os.path.abspath(dir_path) + '\\' + text
                self._oofem_idset_outtypes = generate_OOFEM_input_file(self._oofem_input_filepath, self.femmdl, eltypes,
                                                                       True)
            pass

    def on_analyse_oofem(self):
        if os.path.exists(self._oofem_input_filepath):
            if False:
                analyse_with_OOFEM(self._oofem_input_filepath, self._oofem_idset_outtypes, self.femmdl)
            else:
                dict_idset_outtypes: Dict[OutputElementType, int] = {}  # output
                dict_idset_outtypes[4] = OutputElementType.Shell
                analyse_with_OOFEM(self._oofem_input_filepath, dict_idset_outtypes, self.femmdl)

    def get_selected_element(self):
        return self.femmdl.selected_entitiy

    def get_selected_node(self):
        return None

    def get_selected_load(self):
        return None
    @property
    def app(self):
        return self._app

    @property
    def mainwin(self):
        return self.app.mainFrame

    @property
    def glwin(self):
        return self.mainwin.glWin


class FEMImporter(IOHandler):
    def __init__(self):
        super().__init__()

    def do_import_geometry(self, fileName):
        if len(fileName) < 1:
            return
        filename_noext, file_extension = os.path.splitext(fileName)
        if file_extension not in self.getImportFormats():
            return
        QApplication.changeOverrideCursor(QCursor(Qt.WaitCursor))
        femmdl = None
        os.chdir(os.path.dirname(fileName))
        if file_extension == '.in':
            if GeoOOFEM.isOOFEMFile(fileName):
                femmdl = GeoOOFEM(fileName,filename_noext)
        logging.debug("do_import_geometry: {}".format(femmdl.guid))
        QApplication.restoreOverrideCursor()
        if femmdl is not None:
            return femmdl

    def getImportFormats(self):
        return (".in",".xml",".sin")


class DialogMeshControl(QDialog):
    def __init__(self, parent):
        super().__init__(parent)
        self.mainwin = parent


        self.mainLayout = QVBoxLayout()

        self.setLayout(self.mainLayout)
        self.currfemmdl:GeoFEM = None
        self.mc=0
        self.txtMCupptresh=0
        self.txtMClowtresh = 0


    def createButton(self, text, member):
        button = QPushButton(text)
        button.clicked.connect(member)
        return button

    def applyMeshControl(self):

        self.mc.lowertreshold=float(self.txtMClowtresh.text())
        self.mc.uppertreshold = float(self.txtMCupptresh.text())
        self.currfemmdl.regenerateusingcolor()
        self.currfemmdl.emit_geometries_rebuild()

    def setCurrentGeoFEM(self, geofem):
        self.currfemmdl = geofem
        self.mc = geofem.mc
        if len(self.windowTitle()) < 10:
            self.setWindowTitle("Mesh Control")
            self.btnModify = self.createButton("&Modify", self.applyMeshControl)
            self.btnModify.setFixedWidth(50)
            #self.btnModify.setFixedHeight(20)

            propLayout=QFormLayout()
            self.mainLayout.addLayout(propLayout)
            self.mainLayout.addWidget(self.btnModify)
            self.txtMCupptresh = QLineEdit()

            propLayout.addRow("&Upper treshold:", self.txtMCupptresh)

            self.txtMClowtresh = QLineEdit()


            propLayout.addRow("&Lower treshold:", self.txtMClowtresh)

        self.txtMClowtresh.setText(str(self.mc.lowertreshold))
        self.txtMCupptresh.setText(str(self.mc.uppertreshold))

class DialogSelectionInfo(QDialog):
    def __init__(self, parent):
        super().__init__(parent)
        self.mainwin = parent
        self.mainLayout = QVBoxLayout()
        self.setLayout(self.mainLayout)
        self.currfemmdl:GeoFEM = None
        self.si=0
        self.txtInfo = QTextEdit()
        self.mainLayout.addWidget(self.txtInfo)
        self.setWindowTitle("Selection info")


    def setCurrentSelection(self, geofem:GeoFEM, si:SelectionInfo):
        self.currfemmdl = geofem
        self.si = si
        if geofem.selected_entitiy is not None:
            self.txtInfo.setText(str(geofem.selected_entitiy.get_info()))
        elif geofem.selected_entitiy is not None:
            self.txtInfo.setText(str(geofem.selected_entitiy))

def createCommand():
    return FEMCommand()