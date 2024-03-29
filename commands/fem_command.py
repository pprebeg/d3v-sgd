try:
    from PySide6.QtWidgets import QApplication, QMenu
    from PySide6.QtWidgets import QDialog, QPushButton,QGridLayout,QToolTip,QCheckBox,QComboBox
    from PySide6.QtWidgets import QTreeView,QMainWindow,QVBoxLayout,QHBoxLayout,QSizePolicy
    from PySide6.QtWidgets import QTreeWidget,QTreeWidgetItem,QDockWidget,QWidget,QGroupBox
    from PySide6.QtWidgets import QInputDialog,QLineEdit,QFormLayout,QTextEdit
    from PySide6.QtWidgets import QLabel, QSpinBox,QToolBar
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
    from mainwin import MainFrame
    from femdir.oofemwrap import *

    from femdir.geooofem import GeoOOFEM
    from femdir.geofem import GeoFEM
    from core import geometry_manager as manager
    import logging
    from femdir.oofemanalysis import is_oofempy_loaded
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
        self.view_type_key = "Element Type"
        self.femmdl: GeoFEM = null
        self.si = 0
        self._curlc = null
        self._combo_lc: QComboBox = null
        self._menuMain = null
        self._show_oofem_analysis_menu = False
        # path to oofem
        self._oofem_input_filepath = ''
        self._oofem_idset_outtypes = null
        #oofempymodfind = importlib.util.find_spec("oofempy") is not None
        if is_oofempy_loaded():
            self._show_oofem_analysis_menu = True

        self.add_toolbars_and_menus()
        self.meshctrl = DialogMeshControl(self.mainwin)
        self.analysis_output = None

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
        self._menuOOFEM = QMenu("&OOFEM")
        self._menuMain.addMenu(self._menuOOFEM)
        self._menuOOFEMresults = QMenu("&Analysis results")
        self._menuMain.addMenu(self._menuOOFEMresults)

        menuExportOOFEM = self._menuOOFEM.addAction("&Export input file")
        menuExportOOFEM.triggered.connect(self.onExportOOFEMin)

        menuMeshControl = self._menuMain.addAction("View Control")
        menuMeshControl.triggered.connect(self.onMeshControl)

        menuAnalysisResults = self._menuOOFEMresults.addAction("&Show displacement...")
        menuAnalysisResults.triggered.connect(self.onShowDisplacement)

        menuAnalysisResults = self._menuOOFEMresults.addAction("&Print element numbers")
        menuAnalysisResults.triggered.connect(self.onPrintElementNumbers)

        menuAnalysisResults = self._menuOOFEMresults.addAction("&Print displacements")
        menuAnalysisResults.triggered.connect(self.onPrintDisplacements)

        menuAnalysisResults = self._menuOOFEMresults.addAction("&Maximum displacement")
        menuAnalysisResults.triggered.connect(self.onPrintMaxDisplacement)

        mb.addMenu(self._menuMain)
        self.inittoolbar()
        if self._show_oofem_analysis_menu:
            self.add_OOFEM_analysis_toolbars_and_menus()

    def get_node_displacement(self):
        """
        :return: Nodal displacement dictionary
        """
        nodal_displacement = {}
        results = self.analysis_output
        node_out, shell_out, beam_out, react_out, node_id_out = results

        loadcase_id = 0
        node_ids = node_id_out[loadcase_id]
        displacement = node_out[loadcase_id]
        n_nodes = len(node_ids)
        for node in range(0, n_nodes):
            nodal_displacement[node_ids[node]] = displacement[node]
        return nodal_displacement

    def onShowDisplacement(self):
        displacement_gui = DialogShowDeformation(self.mainwin, self)
        displacement_gui.show_deformation_gui()
        displacement_gui.exec()

    def onPrintElementNumbers(self):
        num_quad = 0
        num_beam = 0
        num_tria = 0
        num_total = self.femmdl.num_elements
        for element in self.femmdl.elements.values():
            if element.get_type() is FEMElementType.Quad:
                num_quad += 1
            elif element.get_type() is FEMElementType.Beam:
                num_beam += 1
            elif element.get_type() is FEMElementType.Tria:
                num_tria += 1
        print("Total number of elements:", num_total,
              ", number of quads:", num_quad,
              ", number of beams:", num_beam,
              ", number of triangles:", num_tria)

    def onPrintDisplacements(self):
        node_displacements = self.get_node_displacement()
        for key, val in node_displacements.items():
            print("Node ID:", key, ",nodal displacement:", val)

    def onPrintMaxDisplacement(self):
        node_displacements = self.get_node_displacement()
        z_displacement = {}
        for key, val in node_displacements.items():
            z_displacement[key] = val[2]
        max_displacement = min(z_displacement.values())
        print("Maximum displacement:", "{:.2f}".format(max_displacement), "mm")

    def add_OOFEM_analysis_toolbars_and_menus(self):
        menu_anylyse_oofem = self._menuOOFEM.addAction("&Analyse")
        menu_anylyse_oofem.triggered.connect(self.on_analyse_oofem)

    def init_combo_lc(self)->QComboBox:
        combo = QComboBox()
        self.reset_combo_lc_items(combo)
        combo.currentIndexChanged.connect(self.loadcasechanged)
        return combo
    def reset_combo_lc_items(self,combo:QComboBox):
        combo.clear()
        if self.femmdl == None or self.femmdl.loadcases ==None:
            pass
        else:
            for lcID in sorted(self.femmdl.loadcases):
                lc=self.femmdl.getLoadCase(lcID)
                combo.addItem(str(lc.id)+': '+ lc.name,lc)
            combo.setCurrentIndex(0)

    def inittoolbar(self):
        self._combo_lc = self.init_combo_lc()
        #self.mainwin.addToolBarBreak(Qt.RightToolBarArea)
        #self.mainwin.addToolBarBreak(Qt.TopToolBarArea)
        path = "modules\\commands\\femdir\\images"
        #Top toolbar
        area = Qt.TopToolBarArea
        tb = QToolBar()
        self.mainwin.addToolBar(area, tb)
        # add loadcase combo
        tb.addWidget(QLabel("LC: "))
        tb.addWidget(self._combo_lc)
        return

        imgpath =path+'\\test.svg'
        act = QAction(QIcon(imgpath), 'TestT', mw)
        act.setShortcut('Ctrl+T')
        act.triggered.connect(self.testaction)
        tb.addAction(act)

        toolButton = QToolButton()
        toolButton.setIcon(QIcon(imgpath))
        toolButton.setCheckable(True)
        toolButton.toggled.connect(self.testaction)
        tb.addWidget(toolButton)


        #tb = QToolBar()

        #mw.addToolBar(area, tb)

        #Right toolbar
        area = Qt.RightToolBarArea
        imgpath = path + '\\test.svg'
        act = QAction(QIcon(imgpath), 'TestR', mw)
        act.setShortcut('Ctrl+R')
        act.triggered.connect(self.testaction)
        tb = QToolBar()
        tb.addAction(act)
        self.mainwin.addToolBar(area, tb)
        pass
    def loadcasechanged(self,index:int):

        text = self._combo_lc.itemText(index)
        self._curlc = self._combo_lc.itemData(index, Qt.UserRole)
        if self._curlc is not None and self.view_type_key != "Element Type":
            self.prep_for_visualization(self.view_type_key,self._curlc.id)
        pass


    @Slot()
    def onGeometryCreated(self, geometries:List[Geometry]):
        for g in geometries:
            if isinstance(g, GeoFEM):
                is_new_fem = self.femmdl is None
                if not is_new_fem:
                    if self.femmdl.guid is not g.guid:
                        is_new_fem=True
                if is_new_fem:
                    self.femmdl = g
                    self.reset_combo_lc_items(self._combo_lc)
                    self.add_viewtype_and_results_menus()
                    self.view_type_key = "Element Type"
                    # self.addViewTypeAndResulsMenus(self.femmdl,self._menuViewType, self._menuResults)
                    break
            if isinstance(g, GeoOOFEM):
                self._oofem_input_filepath = g.filename

    def add_viewtype_and_results_menus(self):
        femmdl = self.femmdl
        menuViewType = self._menuViewType
        menuResults = QMenu("&Results")
        menuViewType.clear()
        menuResul = menuViewType.addAction("Element Type")
        menuResul.triggered.connect(self.onElementType)
        for key, atr in femmdl.attrib_val_functions.items():
            menuResul = menuViewType.addAction(key)
            menuResul.triggered.connect(self.onColorControlMenu)
        if len(femmdl.element_results) > 0:
            self._menuViewType.addMenu(menuResults)
        for key, res in femmdl.element_results.items():
            menuResul = menuResults.addAction(key)
            menuResul.triggered.connect(self.onColorControlMenu)

    def onElementType(self):
        self.femmdl.regenerate()
        self.femmdl.emit_geometries_rebuild()
        self.view_type_key="Element Type"

    def onColorControlMenu(self):
        action = QObject.sender(self)
        key = action.iconText()
        self.view_type_key = key
        lcID=-1
        if self._curlc is not None:
            lcID = self._curlc.id
        self.prep_for_visualization(key,lcID)
        pass

    def prep_for_visualization(self,view_key,lcID):
        self.femmdl.prepareModelForVisualization(view_key, lcID)
        self.meshctrl.setCurrentGeoFEM(self.femmdl)
        self.femmdl.emit_geometries_rebuild()

    @Slot()
    def registerSelection(self, si: SelectionInfo):
        self.si = si
        if self.femmdl is not None:
            if si.isEmpty():
                self.femmdl.unselect()
            else:
                self.femmdl.onSelected(si)
                pos: QPoint = self.mainwin.pos()
                pos.setX(pos.x() + self.mainwin.glWin.dragInfo.wStartPos.x() + 20)
                pos.setY(pos.y() + self.mainwin.glWin.size().height() - self.mainwin.glWin.dragInfo.wStartPos.y())
                msg = self.femmdl.selected_entitiy.get_info()
                QApplication.instance().clipboard().setText(str(msg))
                #QToolTip.showText(pos, msg, msecShowTime=10000)
                QToolTip.showText(pos, msg)

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
            text, ok = QInputDialog.getText(self.mainwin, 'GeoOOFEM Export input file', 'Use next file name?',
                                            QLineEdit.Normal, 'oofem_exp_input.in')
            if ok:
                self._oofem_input_filepath = os.path.abspath(dir_path) + '\\' + text
                self._oofem_idset_outtypes = generate_OOFEM_input_file(self._oofem_input_filepath, self.femmdl, eltypes,
                                                                       True)
            pass

    def on_analyse_oofem(self):
        if isinstance(self.femmdl, GeoFEM):
            eltypes = get_initial_element_types_dict()
            dir_path = os.path.dirname(self.femmdl.get_input_file_path())
            self._oofem_input_filepath = os.path.abspath(dir_path) + '\\' + 'oofem_run_file.in'
            self._oofem_idset_outtypes,self._femmdl_idgroup_outtypes = generate_OOFEM_input_file(self._oofem_input_filepath, self.femmdl, eltypes,
                                                                       True)
            results = analyse_with_OOFEM(self._oofem_input_filepath, self._oofem_idset_outtypes, self.femmdl,self._femmdl_idgroup_outtypes)
            self.analysis_output = results
            self.add_viewtype_and_results_menus()

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


class DialogShowDeformation(QDialog):
    def __init__(self, parent, fcm: FEMCommand):
        super().__init__(parent)
        self.mainwin = parent
        self.fcm = fcm
        self.femmdl = fcm.femmdl
        self.setWindowTitle("Display deformed model")
        # self.setFixedSize(450, 310)
        self.main_layout = QGridLayout()
        self.setLayout(self.main_layout)

        self.scale_factor = QSpinBox(self)

    def def_scale_label(self):
        label = QLabel(self)
        label.setText("Deformation scale factor")
        self.main_layout.addWidget(label, 0, 0)

    def def_scale_input(self):
        self.scale_factor.setMinimum(1)
        self.scale_factor.setValue(10)
        self.scale_factor.setMaximum(10000)
        self.main_layout.addWidget(self.scale_factor, 1, 0)

    def button_display(self):
        generate_button = QPushButton("Display deformation", self)
        self.main_layout.addWidget(generate_button, 1, 1)
        generate_button.clicked.connect(self.regenerate_deformation)

    def show_deformation_gui(self):
        self.def_scale_label()
        self.def_scale_input()
        self.button_display()

    def regenerate_deformation(self):
        nodal_displacement = self.fcm.get_node_displacement()
        def_scale = self.scale_factor.value()
        self.femmdl.regenerate_deformation(nodal_displacement, def_scale)
        if self.femmdl is not None:
            manager.remove_geometry([self.femmdl])
            manager.add_geometry([self.femmdl])
            manager.show_geometry([self.femmdl])
        self.close()
        pass


def createCommand():
    return FEMCommand()