from typing import List, Dict
from femdir.geofementity import *
from femdir.oofemin import OOFEMAnalysisModel
#d3v imports
from geometry import BBox
from selection import SelectionInfo
#from extendedgeometry import ConnectedModel,ModelBasedGeometry


#helper methods
def get_string_from_list(inlist: List):
    n = len(inlist)
    line = ''
    if n > 0:
        for i in range(n - 1):
            line += str(inlist[i]) + ', '
        line += str(inlist[-1])
        return line
    return line

def str_to_float(in_string):
    try:
        fnum = float(in_string)
        return fnum
    except:
        return np.nan

def get_float_list_from_string(line: str,delim:str = ' '):
    splited = line.split(delim)
    flist = None
    if splited is not None:
        flist = [0.0]*len(splited)
        for i in range(len(splited)):
            flist[i] = str_to_float(splited[i])

    return flist


#class GeoFEM(ConnectedModel):
class GeoFEM():
    def __init__(self):
        super().__init__()
        self.nodes = {}
        self.elements = {}
        self.materials = {}
        self.properties = {}
        self.stiflayouts = {}
        self.groups = {}
        self._boundaryconditions = {}
        self._loads = {}
        self.loadcases: Dict[int,LoadCase] = {}
        self.meshcontrol = 0
        self.units= Units()
        self.mc = MeshControl()
        self.mas = MaestroElementAssociation()
        self.element_results= {}
        self.element_vertex_results = {}
        self.vertex_results = {}
        self.model_results = {}
        self.result_name=""
        self.is_node_result=False
        self.is_element_result = False
        self.attrib_val_functions ={}
        self.populateAtribValFunctionsDictionary()
        self.minValue=0
        self.maxValue=0
        self.numDiffValues =0
        self.valueIndexColor = {}

        self.drawLegend = False
        self.legendValues = []
        self.legendColors = []
        self.legendTitle = ""
        self.fixed_color_list = self.initColorList()
        self.max_fixed_colors= len(self.fixed_color_list)
        self.facetoelement:Dict[int,int] = {}
        self.nodetoelement:Dict[int,int] = {}
        #geometry handling
        #self._allfegeo= ModelBasedGeometry(self)
        #self._allnodgeo = ModelBasedGeometry(self)
        self._allfegeo = None
        self._allnodgeo = None
        self._allnodgeo._show_mesh_wireframe = False
        self.addgeometry(self._allfegeo)
        self.addgeometry(self._allnodgeo)
        self._selected_entitiy = None
        self._selected_node = None
        self._additional_properties:Dict[FEMElementType,Dict[str,float]] = {}
        self._oofem: OOFEMAnalysisModel = None
        pass

    @property
    def oofem(self)->OOFEMAnalysisModel:
        return self._oofem

    @oofem.setter
    def oofem(self, value:OOFEMAnalysisModel):
        self._oofem = value

    def get_input_file_path(self):
        return ""

    def add_additional_property(self,key,value):
        self._additional_properties[key]=value #value is Dictionary

    @property
    def selected_entitiy(self)->GeoEntity:
        return self._selected_entitiy

    def unselect(self):
        self._selected_entitiy = None

    def get_element_for_face(self,fh)->Element:
        return self.facetoelement.get(fh)

    def get_node_for_face(self,fh)->Node:
        return self.nodetoelement.get(fh)

    def onSelected(self, si:SelectionInfo):
        if si.geometry is self._allfegeo:
            sfh = si.getFace()
            el = self.get_element_for_face(sfh)
            if el != None:
                self._selected_entitiy = el
                if len(el.face_handles) > 1:
                    for fh in el.face_handles:
                        if fh != sfh:
                            si.allfaces.append(fh)
        elif si.geometry is self._allnodgeo:
            sfh = si.getFace()
            nod = self.get_node_for_face(sfh)
            if nod != None:
                self._selected_entitiy = nod
                if len(nod.face_handles) > 1:
                    for fh in nod.face_handles:
                        if fh != sfh:
                            si.allfaces.append(fh)
        pass
    def prepareModelForVisualization(self,key):
        if __debug__:
            ts = time.perf_counter()
        self.minValue = float("inf")
        self.maxValue = float("-inf")
        self.numDiffValues = 0
        self.valueIndexColor.clear()

        fatrib= self.attrib_val_functions.get(key)

        if fatrib == None:
            self.doResultValue(key)
        else:
            fatrib(key)
        self.mc.lowertreshold=self.minValue
        self.mc.uppertreshold=self.maxValue
        self.mc.viewtype = ViewType.face_colors
        self.legendTitle=key
        if __debug__:
            dt = time.perf_counter() - ts
            print("GEO FEM Time to prepare data for mesh, s:", dt)
        self.regenerateusingcolor()

    # region Colors handling
    def initColorList(self):
        colors = [[0, 0, 255,255],[128, 0, 128,255],[222, 184, 135,255],[255, 165, 0,255],[0, 255, 0,255],
                  [ 0, 128, 0,255],[128, 0, 0,255],[255, 0, 0,255],[255, 192, 203,255],[222, 184, 135,255],
                  [255, 165, 0,255],[255, 127, 80,255],[128, 128, 0,255],[255, 255, 0,255],[245, 245, 220,255],
                  [0, 255, 0,255],[ 0, 128, 0,255],[245, 255, 250,255],[0, 128, 128,255],[0, 255, 255,255],
                  [0, 0, 128,255],[230, 230, 250,255],[255, 0, 255,255],[205, 133, 63,255]]
        floatColors = []
        for color in colors:
            floatColors.append([x / 255 for x in color])
        return floatColors
    def getContinuousColor(self, v, vmin, vmax):
        color = [1.0, 1.0, 1.0, 1.0]
        if v > self.mc.uppertreshold:
            return self.mc.getUpperTresholdColor()
        elif v < self.mc.lowertreshold:
            return self.mc.getLowerTresholdColor()
        vmin   = max(vmin,self.mc.lowertreshold)
        vmax = min(vmax, self.mc.uppertreshold)

        if v < vmin:
            v = vmin
        if v > vmax:
            v = vmax
        dv = vmax - vmin

        if (v < (vmin + 0.25 * dv)):
            color[0] = 0
            color[1] = 4 * (v - vmin) / dv
        elif (v < (vmin + 0.5 * dv)):
            color[0] = 0
            color[2] = 1 + 4 * (vmin + 0.25 * dv - v) / dv
        elif (v < (vmin + 0.75 * dv)):
            color[0] = 4 * (v - vmin - 0.5 * dv) / dv
            color[2] = 0
        else:
            color[1] = 1 + 4 * (vmin + 0.75 * dv - v) / dv
            color[2] = 0
        return color

    def getColorFromList(self, v, vmin, vmax):
        if v > self.mc.uppertreshold:
            return self.mc.getUpperTresholdColor()
        elif v < self.mc.lowertreshold:
            return self.mc.getLowerTresholdColor()

        index=self.getValueColorIndex(v)
        color = self.fixed_color_list[index]
        return color


    def prepContColorLegend(self,fun_getcolor,minVal, maxVal,nColor):
        self.legendValues.clear()
        self.legendColors.clear()
        self.drawLegend = True
        minVal = max(minVal, self.mc.lowertreshold)
        maxVal = min(maxVal, self.mc.uppertreshold)
        legendValues=np.linspace(minVal,maxVal,nColor)
        for x in legendValues:
            self.legendValues.append(f"{x:.4g}")
        for val in legendValues:
            color = fun_getcolor(val, minVal, maxVal)
            self.legendColors.append(color)

    def prepListColorLegend(self,fun_getcolor):
        self.legendValues.clear()
        self.legendColors.clear()
        self.drawLegend = True
        for key,index in self.valueIndexColor.items():
            if  key < self.mc.lowertreshold or key > self.mc.uppertreshold:
                continue
            self.legendValues.append(f"{key:.4g}")
            color = fun_getcolor(key, 0, self.numDiffValues)
            self.legendColors.append(color)

    # endregion
    # region Attribute Value Functions

    def populateAtribValFunctionsDictionary(self):
        self.addAttValFunc('TPL', self.doTPLValue)
        self.addAttValFunc('Material ID', self.doMaterialIDValue)
        self.addAttValFunc('Property ID', self.doPropertyIDValue)

    def addAttValFunc(self, key, f):
        self.attrib_val_functions[key] = f

    def doTPLValue(self,key):
        for el in self.elements.values():
            val= el.onTPLValue()
            self.checkMinMax(val)


    def doMaterialIDValue(self,key):
        for el in self.elements.values():
            val = el.onMatID()
            self.checkMinMax(val)

    def doPropertyIDValue(self, key):
        for el in self.elements.values():
            val = el.onPropID()
            self.checkMinMax(val)


    def doResultValue(self,key):
        self.setValueToItemResults(key)

    # endregion


    def checkMinMax(self,val):
        if val > self.maxValue:
            self.maxValue = val
        if val < self.minValue:
            self.minValue = val

        if self.numDiffValues <= self.max_fixed_colors:
            index = self.getValueColorIndex(val)
            if index == None:
                self.valueIndexColor[val]=self.numDiffValues
                self.numDiffValues=self.numDiffValues+1


    def getValueColorIndex(self,value):
        return self.valueIndexColor.get(value)


    def isElementFaceResult(self):
        return self.is_element_result and (not self.is_node_result)

    def isElementNodeResult(self):
        return self.is_element_result and self.is_node_result

    def isNodeResult(self):
        return (not self.is_element_result) and self.is_node_result

    def setValueToItemResults(self, resultName):
        result = self.element_results.get(resultName)
        if result != None:
            self.is_element_result= True
            self.is_node_result = False
            for key, el in self.elements.items():
                val = el.setFaceValueUsingElementID(result.feres)
                self.checkMinMax(val)
            return

        result = self.element_vertex_results.get(resultName)
        if result != None:
            self.is_element_result= True
            self.is_node_result = True
            return

        result = self.vertex_results.get(resultName)
        if result != None:
            self.is_element_result = False
            self.is_node_result = True
            return

    @property
    def num_nodes(self):
        return len(self.nodes)

    @property
    def num_elements(self):
        return len(self.elements)

    @property
    def num_properties(self):
        return len(self.properties)

    @property
    def num_materials(self):
        return len(self.materials)

    @property
    def num_groups(self):
        return len(self.groups)

    @property
    def num_lc(self):
        return len(self.loadcases)

    @property
    def num_bc(self):
        return len(self._boundaryconditions)

    @property
    def num_loads(self):
        return len(self._loads)

    @property
    def boundaryconditions(self):
        return self._boundaryconditions

    @property
    def loads(self):
        return self._loads

    def getNode(self,id)->Node:
        return self.nodes.get(id)
    def getElement(self,id)->Element:
        return self.elements.get(id)
    def getMaterial(self,id)->Material:
        return self.materials.get(id)
    def getProperty(self,id)->Property:
        return self.properties.get(id)
    def getStiffLayout(self,id)->StiffLayoutProperty:
        return self.stiflayouts.get(id)
    def getLoadCase(self, id)->LoadCase:
        return self.loadcases.get(id)
    def addLoadCase(self, id,name:str)->LoadCase:
        lc= LoadCase(id,name)
        self.loadcases[id] = lc
        return lc
    def addNode(self,item):
        self.nodes[item.id]=item
    def addElement(self,item):
        self.elements[item.id]=item
    def addProperty(self,item):
        self.properties[item.id]=item
    def addMaterial(self,item):
        self.materials[item.id]=item
    def addStiffLayout(self, item):
        self.stiflayouts[item.id] = item

    def addGroup(self, id:int, group:Group):
        self.groups[id] = group
    def getGroup(self, id:int)->Group:
        return self.groups.get(id)

    def addLoadToLoadcase(self, idLC:int, load:BC_and_Load):
        if not(idLC in self.loadcases):
            self.loadcases[idLC] = LoadCase(idLC,'LC '+str(idLC))
        self.loadcases[idLC].add_load(load)

    def addBoundaryConditionToLoadcase(self, idLC:int, bc:BC_and_Load):
        if not(idLC in self.loadcases):
            self.loadcases[idLC] = LoadCase(idLC,'LC '+str(idLC))
        self.loadcases[idLC].add_boundarycondition(bc)

    def addLoadToAllLoadCasses(self,nlc,load:BC_and_Load):
        idlcmax=0
        for idlc in self.loadcases:
            if idlc > idlcmax:
                idlcmax = idlc
        n=len(self.loadcases)
        # guess idlc of missing loadcases
        while n < nlc:
            idlcmax+=1
            n+=1
            self.loadcases[idlcmax] = LoadCase(idlcmax, 'LC '+ str(idlcmax))
        for idlc in self.loadcases:
            self.loadcases[idlc].add_load(load)

    def addBoundaryCondition(self, bc:BC_and_Load):
        self._boundaryconditions[bc.id] = bc

    def getBoundaryCondition(self,id):
        return self.boundaryconditions.get(id)

    def addLoad(self, load:BC_and_Load):
        self._loads[load.id] = load

    def getLoad(self,id):
        return self.loads.get(id)

    def regenerate(self):
        if __debug__:
            ts = time.perf_counter()
        #mesh= om.TriMesh()
        mesh = om.PolyMesh()
        self.mc.viewtype = ViewType.constant_color
       # mesh.request_face_colors()
        for el in self.elements.values():
            el.updateMesh(mesh, self.facetoelement, self.mc)
        pass
        self._allfegeo.mesh = mesh

        mesh = om.PolyMesh()
        mesh.request_face_colors()
        if self.mc.show_nodes:
            box = self._allfegeo.bbox
            refdim = self.get_box_max_span(box)
            Node.set_sphere(refdim, self.mc)
            for nod in self.nodes.values():
                nod.updateMesh(mesh, self.nodetoelement, self.mc)
            pass
            #colors = mesh.face_colors()
            #colors[:]= [1.0,0.0,0.0,1.0]
        self._allnodgeo.mesh = mesh

        if __debug__:
            dt = time.perf_counter() - ts
            print("GEO FEM Time to regenerate mesh, s:", dt)

    def regenerateusingcolor(self):
        if __debug__:
            ts = time.perf_counter()

        fun_getcolor = self.getColorFromList
        if self.numDiffValues > self.max_fixed_colors:
            fun_getcolor= self.getContinuousColor

        #mesh= om.TriMesh()
        mesh = om.PolyMesh()
        if self.mc.viewtype == ViewType.constant_color or self.mc.viewtype == ViewType.face_colors:
            #mesh.release_vertex_colors()
            mesh.request_face_colors()
        elif self.mc.viewtype == ViewType.face_vertex_colors:
            #mesh.release_face_colors()
            mesh.request_vertex_colors()

        const_color = [0.4, 1.0, 1.0, 1.0]

        for el in self.elements.values():
            el.updateMesh(mesh, self.facetoelement, self.mc, const_color=const_color, fun_getcolor=fun_getcolor,
                          minvalue=self.minValue, maxvalue=self.maxValue)

        self._allfegeo.mesh = mesh


        if __debug__:
            dt = time.perf_counter() - ts
            print("GEO FEM Time to regenerate mesh using color, s:", dt)

        if self.numDiffValues > self.max_fixed_colors:
            self.prepContColorLegend(fun_getcolor,self.minValue, self.maxValue, 12)
        else:
            self.prepListColorLegend(fun_getcolor)


    def get_box_diagonal(self,box:BBox):
        diag = np.linalg.norm(box.maxCoord - box.minCoord)
        return diag

    def get_box_max_span(self,box:BBox):
        span = np.max(box.maxCoord - box.minCoord)
        return span

    def setResultValuesOnElements(self):
        pass








    # def showFaceColorP(self, propDict):
    #     colors = [[0, 0, 255, 255], [128, 0, 128, 255], [222, 184, 135, 255], [255, 165, 0, 255], [0, 255, 0, 255],
    #               [0, 128, 0, 255], [128, 0, 0, 255], [255, 0, 0, 255], [255, 192, 203, 255], [222, 184, 135, 255],
    #               [255, 165, 0, 255], [255, 127, 80, 255], [128, 128, 0, 255], [255, 255, 0, 255], [245, 245, 220, 255],
    #               [0, 255, 0, 255], [0, 128, 0, 255], [245, 255, 250, 255], [0, 128, 128, 255], [0, 255, 255, 255],
    #               [0, 0, 128, 255], [230, 230, 250, 255], [255, 0, 255, 255], [205, 133, 63, 255]]
    #
    #     floatColors = []
    #     for color in colors:
    #         floatColors.append([x / 255 for x in color])
    #     mesh = self.mesh
    #     mesh.request_face_colors()
    #     propColorDict = {}
    #     self.legendValues.clear()
    #     self.legendColors.clear()
    #     self.drawLegend = False
    #     for el in self.element2Face:
    #         idProp=propDict[el]
    #         nuc=len(propColorDict)
    #         indexColor = 0
    #         if idProp in propColorDict:
    #             indexColor=propColorDict[idProp]
    #         else:
    #             propColorDict[idProp]=nuc
    #             indexColor=nuc
    #             self.legendValues.append(str(idProp))
    #             self.legendColors.append(floatColors[indexColor])
    #         for fh in self.element2Face[el]:
    #             mesh.set_color(fh, floatColors[indexColor])
    #         pass
    #     if len(self.legendValues)> 0:
    #         self.drawLegend=True
    #     pass


class Result():
    def __init__(self, name):
        self.name = name
        pass


class GeneralResultsDictionary(Result):
    def __init__(self, name):
        super().__init__(name)
        self.results = {}
        pass

    def appendValue(self, key, value):
        resultList = self.results.setdefault(key, [])
        resultList.append(value)

    def addValues(self, key, values: []):
        self.results[key] = values

    def addValues2(self, key, values: list):
        self.results[key] = values.copy()

    def getValues(self, key):
        return self.results[key]

    def getValue(self, key, index: int):
        return self.results[key][index]

    def appdendListwithResultData(self, x: list, y: list, iresx: int, iresy: int):
        for res in self.results.values():
            x.append(res[iresx])
            y.append(res[iresy])

    def appdendListwithKeyPairedResultData(self, x: list, y: list, iresy: int):
        for key, res in self.results:
            x.append(key)
            y.append(res[iresy])

class GeneralResultsTableModel(Result):
    def __init__(self, name):
        super().__init__(name)
        self.data =[]
        self.column_names =[]
        pass

    def appendName(self, name: str):
        self.column_names.append(name)

    def addRow(self, values: list):
        self.data.append(values)

    def getRowValues(self, row_index):
        return self.data[row_index]

    def getValue(self, row_index, column_index: int):
        return self.data[row_index][column_index]

class ElementResult(Result):
    def __init__(self, name):
        super().__init__(name)
        self.feres = {}
        pass

    def getValue(self, feID):
        return self.feres.get(feID)

    def setValue(self, feID, value):
        self.feres[feID]=value

    def keyExist(self,feID):
        return feID in self.feres


class FEMModelResults:
    def __init__(self, name):
        self.name = name
        pass

    def readOutput(self, path):
        pass

    def setResultsToModel(self, fem: GeoFEM):
        pass


class LusaResults(FEMModelResults):
    def __init__(self, name, mas:MaestroElementAssociation):
        super().__init__(name)
        self.las = LusaElementAssociation()
        self.mas=mas
        self.lers = {} #Lusa element results
        self.modres={}
        pass

    def readOutput(self, path):
        abspath1 = '\\'.join(path.split('\\')[0:-1])
        abspath2 = '/'.join(path.split('/')[0:-1])
        if len(abspath2) > len(abspath1):
            abspath=abspath2 + '/'
        else:
            abspath=abspath1 + '\\'

        abspath_hoggCSD = abspath + 'LUSAhoggCSD.OUT'
        abspath_saggCSD=  abspath + 'LUSAsaggCSD.OUT'

        self.readCSDFile(abspath_hoggCSD,False)
        self.readCSDFile(abspath_saggCSD,True)

        abspath_hogg = abspath + 'LUSAhogg.OUT'
        abspath_sagg = abspath + 'LUSAsagg.OUT'

        iterationResults = GeneralResultsTableModel('Lusa iteration results Sagg')
        self.modres[iterationResults.name]=iterationResults


        iterationResults.appendName('CycleNo Sagg')
        iterationResults.appendName('Moment Sagg, kNm')
        iterationResults.appendName('Curvature Sagg, 1/m')
        iterationResults.appendName('y_NA Sagg, m')

        self.readMainLusaFile(abspath_sagg,True,iterationResults)

        iterationResults = GeneralResultsTableModel('Lusa iteration results Hogg')
        self.modres[iterationResults.name] = iterationResults
        iterationResults.appendName('CycleNo Hogg')
        iterationResults.appendName('Moment Hogg, kNm')
        iterationResults.appendName('Curvature Hogg, 1/m')
        iterationResults.appendName('y_NA Hogg, m')

        self.readMainLusaFile(abspath_hogg, False, iterationResults)

        pass

    def readMainLusaFile(self, path, isSagg, tableResult:GeneralResultsTableModel):
        file = pathlib.Path(path)
        if not file.exists():
            return
        f = open(path, "_r")
        nlines2skip = 0
        isCycleData = False

        for line in f:
            if nlines2skip > 0:
                nlines2skip = nlines2skip - 1
                continue
            line = ' '.join(line.split())
            if line.startswith('*'):
                continue
            if line == "" or line == " ":
                continue

            if 'HULL MODULE RESPONSE DATA' in line:
                nlines2skip = 4
                isCycleData = True
                continue
            if 'ULTIMATE CAPACITY IS' in line:
                f.close()
#                if isSagg:
#                    tableResult.data.reverse()
                return
            sline = line.split(" ")
            if len(sline) == 0:
                continue
            if isCycleData:
                if len(sline) == 4:
                    rowValues=[0]*4
                    rowValues[0]= int(sline[0])
                    rowValues[1] = float(sline[1])
                    rowValues[2] = float(sline[2])
                    rowValues[3] = float(sline[3])
                    tableResult.addRow(rowValues)

        f.close()

    def readCSDFile(self,path,isSagg):
        file=pathlib.Path(path)
        if not file.exists():
            return
        f = open(path, "_r")
        nlines2skip=0
        isSPCdata=False
        isGPCdata = False
        isHCdata = False
        if isSagg:
            collapse_stress = ElementResult('Collapse Stress Sagg')
            collapse_mod = ElementResult('Collapse Mod Sagg')
            collapse_cycle = ElementResult('Collapse Cycle Sagg')
        else:
            collapse_stress = ElementResult('Collapse Stress Hogg')
            collapse_mod = ElementResult('Collapse Mod Hogg')
            collapse_cycle = ElementResult('Collapse Cycle Hogg')

        self.lers[collapse_stress.name]=collapse_stress
        self.lers[collapse_mod.name] = collapse_mod
        self.lers[collapse_cycle.name] = collapse_cycle
        for line in f:
            if  nlines2skip > 0:
                nlines2skip=nlines2skip-1
                continue
            line = ' '.join(line.split())
            if line.startswith('*'):
                continue
            if line == "" or line == " ":
                continue

            if 'Stiffener - Plate Combinations (SPCs)' in line:
                nlines2skip=4
                isSPCdata=True
                continue
            if 'Girder - Plate Combinations (GPCs)' in line:
                nlines2skip=4
                isGPCdata=True
                isSPCdata=False
                continue
            if 'Hard Corners (HCs)' in line:
                nlines2skip = 4
                isGPCdata = False
                isHCdata = True
                continue
            sline = line.split(" ")
            if  len(sline)== 0:
                continue
            el_no_lusa=-1
            if isSPCdata:
                if len(sline) > 5:
                    strakeNo    = int(sline[0])
                    el_no_lusa        = int(sline[1])
                    self.las.addPlate(el_no_lusa, strakeNo)

                    elIDs=self.mas.getPlateElemForStrake(strakeNo)
                    cc_new = int(sline[4])
                    for id_el in elIDs:
                        bAddNew=True
                        cc_old= collapse_cycle.getValue(id_el)
                        if cc_old != None and cc_old < cc_new:
                            pass
                        else:
                            collapse_stress.setValue(id_el, float(sline[2]))
                            collapse_mod.setValue(id_el, float(sline[3]))
                            collapse_cycle.setValue(id_el, cc_new)
            elif isGPCdata:
                if len(sline) > 5:
                    strakeNo    = int(sline[0])
                    el_no_lusa        = int(sline[1])
                    self.las.addSPC(el_no_lusa, strakeNo)

                    elIDs = self.mas.getGirderBeamElemForStrake(strakeNo)
                    elIDsPlate = self.mas.getPlateElemForStrake(strakeNo)
                    if elIDs == None:
                        elIDs = elIDsPlate
                    elif elIDsPlate != None:
                        for el in elIDsPlate:
                            elIDs.append(el)
                    cc_new = int(sline[4])
                    for id_el in elIDs:
                        bAddNew = True
                        cc_old = collapse_cycle.getValue(id_el)
                        if cc_old != None and cc_old < cc_new:
                            pass
                        else:
                            collapse_stress.setValue(id_el, float(sline[2]))
                            collapse_mod.setValue(id_el, float(sline[3]))
                            collapse_cycle.setValue(id_el, cc_new)
            elif isHCdata:
                if len(sline) > 5:
                    endPtNo    = int(sline[0])
                    el_no_lusa        = int(sline[1])
                    self.las.addHC(el_no_lusa, endPtNo)
                    elIDs = self.mas.getPlateElsForEnpoint(endPtNo)
                    cc_new = int(sline[4])
                    for id_el in elIDs:
                        bAddNew = True
                        cc_old = collapse_cycle.getValue(id_el)
                        if cc_old != None and cc_old < cc_new:
                            pass
                        else:
                            collapse_stress.setValue(id_el, float(sline[2]))
                            collapse_mod.setValue(id_el, float(sline[3]))
                            collapse_cycle.setValue(id_el, cc_new)

        f.close()
        pass

    def setResultsToModel(self, fem: GeoFEM):
        pass