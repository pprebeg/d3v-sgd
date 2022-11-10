from femdir.geofem import *
from femdir.oofemin import *

class GeoOOFEM (GeoFEM):
    def __init__(self,fileName,name=''):
        self.filename = fileName
        super().__init__(name)
        self.outfilename= ""
        if fileName != '':
            self.readModel()

    @staticmethod
    def getGeoFEMElementType(keyword)->FEMElementType:
        beams = {'beam3d'}
        rods= {'truss3d'}
        trias = {'tr_shell02', "dktplate", "trplanestress2d", "trplanestressrotallman", "trplanestrrot",'tr_shell01'}
        quads = {'shellqd42','shellqd41', 'mitc4shell', "planestress2d",'linquad3dplanestress'}
        keyword = keyword.lower()
        if keyword in quads:
            return FEMElementType.Quad
        if keyword in beams:
            return FEMElementType.Beam
        if keyword in trias:
            return FEMElementType.Tria
        if keyword in rods:
            return FEMElementType.Rod
        return FEMElementType.Unknown

    def get_input_file_path(self):
        return self.filename

    def readModel(self):
        self.readInputFile()
        self.regenerate()
        return

    @staticmethod
    def isOOFEMFile(filepath):
        return True

    def getoutfilepath(self):
        abspath = '\\'.join(self.filename.split('\\')[0:-1])
        if len(abspath) < 1:
            abspath = '/'.join(self.filename.split('/')[0:-1])
        return abspath + '/' + self.outfilename

    def readInputFile(self):
        if __debug__:
            print('GeoOOFEM read started: '+ self.filename)
        mdl = self
        null = None
        dict_idprop_addcharvalues = {}
        dict_idprop_elkeyword = {}
        drgm = OOFEMInputFile(self.filename) # GeoOOFEM data record group manager
        # Analysis record
        drg = drgm.getRecordGroup(AnalysisRecord.getkeyname())  # data record group
        if drg is not None:
            keyword,nlc = drg.getNeutralFormatRecordData(0)
            if nlc == 'Unknown':
                print('GeoOOFEM input file read stopped due to unhandled Analysis Type')
                return
        propset_idnod = {} # key idset, value [] idelements
        propset_idel = {}  # key idset, value [] idelements
        dict_el_prop = {} #key idel, value idcs

        # Sets
        drg = drgm.getRecordGroup(SetRecords.getkeyname())  # data record group
        if drg is not None:
            n = drg.getNumRecords()
            for i in range(n):
                data = drg.getNeutralFormatRecordData(i)
                if data[0] == 'nodes':
                    propset_idnod[data[1]]=data[2]
                elif data[0] == 'elements':
                    propset_idel[data[1]]=data[2]
        # Materials
        drg = drgm.getRecordGroup(MaterialRecords.getkeyname()) # data record group
        n =  drg.getNumRecords()
        for i in range(n):
            data = drg.getNeutralFormatRecordData(i)
            if data is None:
                print('Unsupported material: ' + str(data))
            elif data[0].lower()== 'isole':
                self.addOOFEMMaterial(data[1:])
            else:
                print('Unsupported material: '+str(data))
        # Properties
        drg = drgm.getRecordGroup(CrossSectionRecords.getkeyname())  # data record group
        if drg is not null:
            n = drg.getNumRecords()
            for i in range(n):
                data = drg.getNeutralFormatRecordData(i)
                idset = data[len(data) - 1]
                if data[0] is null:
                    continue
                if data[0] == 'plate':
                    self.addOOFEMPlateProperty(data[1:])
                elif data[0] == 'plate_add':
                    self.addOOFEMPlateProperty(data[1:5])
                    dict_idprop_addcharvalues[int(data[1])] = data[5]
                    idset = data[len(data) - 2]
                elif data[0] == 'k_beam':
                    self.addOOFEMBeamProperty(data[1:])
                elif data[0] == 'k_rod':
                    self.addOOFEMRodProperty(data[1:])
                if idset is not null: # if set exist
                    idprop= data[1]
                    for idel in propset_idel[idset]:
                        dict_el_prop[idel]=idprop
                    
        # Nodes
        drg = drgm.getRecordGroup(NodeRecords.getkeyname())  # data record group
        if drg is not null:
            n = drg.getNumRecords()
            rigidarmnodes=[]
            for i in range(n):
                data = drg.getNeutralFormatRecordData(i)
                if data[0]=='node':
                    self.addOOFEMNode(data[1],data[2])
                    idata = min(3,len(data))
                    while idata < len(data):
                        if data[idata]== 'lcs':
                            idata += 2
                elif data[0]=='rigidarmnode':
                    rigidarmnodes.append(data)
                else:
                    print('Unknown node type:' + data[0] + ' id = ' + str(data[1]))

            for data in rigidarmnodes:
                self.addOOFEMRigidArmNode(data)


        # Elements
        drg = drgm.getRecordGroup(ElementRecords.getkeyname())  # data record group
        if drg is not null:
            n = drg.getNumRecords()
            for i in range(n):
                data = drg.getNeutralFormatRecordData(i)
                elkeyword =data[0].lower()
                eltype = GeoOOFEM.getGeoFEMElementType(elkeyword)
                idel=data[1]
                if len(data) > 3:
                    idprop = data[3]
                else:
                    idprop = dict_el_prop[idel]
                elkeywordsset = dict_idprop_elkeyword.get(idprop)
                if elkeywordsset is null:
                    elkeywordsset = set()
                    dict_idprop_elkeyword[idprop] = elkeywordsset
                elkeywordsset.add(elkeyword)
                if eltype == FEMElementType.Quad:
                    el=QuadElement()
                    self.addOOFEMElement(el, idel, idprop, data[2])
                elif eltype == FEMElementType.Beam:
                    el=BeamElement()
                    ref_node = mdl.nodes[data[2][-1]]
                    orient = BeamOrientationNode(ref_node) # last node is referent node
                    el.set_beam_orientation(orient)
                    self.addOOFEMElement(el, idel, idprop, data[2][:-1])
                elif eltype == FEMElementType.Tria:
                    el=TriaElement()
                    self.addOOFEMElement(el,idel, idprop, data[2])
                elif eltype == FEMElementType.Rod:
                    el=RodElement()
                    self.addOOFEMElement(el, idel, idprop, data[2])
                else:
                    print('Unknown element type:' + data[0]+' _id = ' + str(data[1]))
        # Add groups
        for id,idnodes in propset_idnod.items():
            self.addOOFEMNodeGroup(id,idnodes)
        for id,idelements in propset_idel.items():
            self.addOOFEMElementGroup(id,idelements)
        # Time function
        timeFunc2Loadcase={}
        timeFuncAllLoadcase = set()
        drg = drgm.getRecordGroup(TimeFunctionRecords.getkeyname())  # data record group
        if drg is not null:
            n = drg.getNumRecords()
            for i in range(n):
                id,idlc = drg.getNeutralFormatRecordData(i)
                if idlc == -999:
                    timeFuncAllLoadcase.add(id)
                elif idlc > 0 :
                    timeFunc2Loadcase[id]=idlc
                else:
                    print('Unknown time function to loadcase data:' + data[0] + ' _id = ' + str(data[1]))

        # Boundary and LoadCase
        drg = drgm.getRecordGroup(BCandLoadRecords.getkeyname())  # data record group
        if drg is not null:
            n = drg.getNumRecords()
            for i in range(n):
                neutral = drg.getNeutralFormatRecordData(i)
                idfun = neutral[3]
                idLC = -1
                if idfun in timeFuncAllLoadcase:
                    idLC = -999
                elif idfun in timeFunc2Loadcase:
                    idLC = timeFunc2Loadcase[idfun]
                else:
                    print('Unknown time function to loadcase data:  ' + drg.getRecordLine(i))
                keyword = neutral[0]
                if keyword == 'nodalload':
                    keyword, id, idset, idfun, ldofs, lvals = neutral
                    self.addOOFEMGroupLoad(id,nlc,idLC,lvals,ldofs,idset)
                elif keyword == 'boundarycondition':
                    keyword, id, idset, idfun, ldofs, lvals = neutral
                    self.addOOFEMGroupBC(id,idLC,lvals,ldofs,idset)
                elif keyword == 'constantsurfaceload':
                    keyword, id, idset, idfun, press = neutral
                    self.addOOFEMGroupPressureLoad(id,nlc,idLC,press,idset)
                elif keyword == 'deadweight':
                    keyword, id, idset, idfun, lvals = neutral
                    self.addOOFEMAccelerationLoad(id,nlc,idLC,lvals,idset)
                else:
                    print('Unhandled Boundary and LoadCase record:  ' + drg.getRecordLine(i))
        #handle additional properties
        for idprop,elkeywordsset in dict_idprop_elkeyword.items():
            for elkeyword in elkeywordsset:
                addcharvalues = dict_idprop_addcharvalues.get(idprop)
                if addcharvalues is not null:
                    mdl.add_additional_property(elkeyword,addcharvalues)


        if __debug__:
            print('GeoOOFEM read ended!')

    def addOOFEMNode(self,id,cords):
        node = Node()
        node.init(id, cords[0], cords[1], cords[2])
        self.addNode(node)

    def addOOFEMRigidArmNode(self,data):
        (keyword, id, coords,idmaster,dofidmask)= data
        masternode = self.getNode(idmaster)
        node = NodeRigidArm()
        node.init(id, coords[0], coords[1], coords[2],masternode,dofidmask)
        self.addNode(node)

    def addOOFEMElement(self, elem:Element, id, idProp, nodeIds):
        elem.init(id)
        elem.property = self.getProperty(idProp)
        for idNod in nodeIds:
            node = self.getNode(idNod)
            elem.addNode(node)
        self.addElement(elem)

    def addOOFEMNodeGroup(self, id:int,idnodes:List[int]):
        group= NodeGroup()
        group.init(id, 'Nodal_group_' + str(id))
        for idnod in idnodes:
            node = self.getNode(idnod)
            group.add_item(node)
        self.addGroup(id,group)

    def addOOFEMElementGroup(self, id:int, idelements:List[int]):
        group= ElementGroup()
        group.init(id, 'Element_group_' + str(id))
        for idel in  idelements:
            el = self.getElement(idel)
            group.add_item(el)
        self.addGroup(id,group)

    def addOOFEMMaterial(self, neudata):
        id = neudata[0]
        E = neudata[1]
        ni = neudata[2]
        rho = neudata[3]
        material = Material()
        material.init(id, 'Mat_' + str(id))
        material.E = E
        material.ni = ni
        material.rho = rho
        self.addMaterial(material)

    def addOOFEMPlateProperty(self, neudata):
        (id, idmat, tp, idset) = neudata
        prop = PlateProperty()
        prop.init(id, 'Plate_prop_' + str(id))
        prop.tp = tp
        prop.material = self.getMaterial(idmat)
        self.addProperty(prop)

    def addOOFEMBeamProperty(self, neudata):
        (id, idmat, sfiff_list, idset)= neudata
        sfiff_list_names = ['area', 'iy', 'iz', 'ik', 'beamshearcoeff', 'shearareay', 'shearareaz']
        prop = StifnessBeamProperty()
        prop.init(id, 'Beam_prop_' + str(id))
        for i in range(len(sfiff_list)):
            if sfiff_list[i] > 0.0:
                prop.add_stiffness_characteristic(sfiff_list_names[i],sfiff_list[i])
        prop.material = self.getMaterial(idmat)
        self.addProperty(prop)

    def addOOFEMRodProperty(self, neudata):
        (id, idmat, area, idset)= neudata
        prop = RodProperty()
        prop.init(id, 'Rod_prop_' + str(id))
        prop.area = area
        prop.material = self.getMaterial(idmat)
        self.addProperty(prop)

    def addOOFEMGroupLoad(self,id,nlc, idLC, values:List[float], dofs:List[int], idset):
        load = GroupDoffBasedLoad(id, self.getGroup(idset), values, dofs)
        self.addLoad(load)
        if idLC > 0:
            self.addLoadToLoadcase(idLC, load)
        elif idLC == -999:
            self.addLoadToAllLoadCasses(nlc,load)

    def addOOFEMGroupPressureLoad(self,id,nlc, idLC, press, idset):
        group=self.getGroup(idset)
        load = GroupPressureLoad(id,group, press)
        self.addLoad(load)
        if idLC > 0:
            self.addLoadToLoadcase(idLC, load)
        elif idLC == -999:
            self.addLoadToAllLoadCasses(nlc,load)

    def addOOFEMAccelerationLoad(self,id,nlc, idLC, values:List[float],idset):
        #assumption is that this type of load is associated to entire model
        load = AccelerationLoad(id,values)
        self.addLoad(load)
        if idLC > 0:
            self.addLoadToLoadcase(idLC, load)
        elif idLC == -999:
            self.addLoadToAllLoadCasses(nlc,load)

    def addOOFEMGroupBC(self,id,idLC:int, values:List[float], dofs:List[int], idset:int):
        bc = GroupNodalBC(id,self.getGroup(idset), values, dofs)
        self.addBoundaryCondition(bc)
        if idLC != -999:
            self.addBoundaryConditionToLoadcase(idLC, bc)



