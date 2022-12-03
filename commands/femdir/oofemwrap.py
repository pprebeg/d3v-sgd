from femdir.oofemin import *
import femdir.geofem
from femdir.geofem import GeoFEM,PlateProperty,BeamProperty,RodProperty,GroupNodalBC,GroupDoffBasedLoad
from femdir.geofem import NodeRigidArm,NodeGroup,ElementGroup,BeamElement,Node,BeamOrientationNode
from femdir.geofem import AccelerationLoad,GroupPressureLoad,ElementResult
from typing import List,Dict
from femdir.oofemenum import CrossSectionProperty as csp


try:
    from femdir.oofemanalysis import OOFEMAnalysisModel
except BaseException as error:
    print('An exception occurred: {}'.format(error))
except:
    print('Unknown exception occurred')



def generate_OOFEM_input_file(file_path:str, mdl:GeoFEM, eltypes:Dict[FEMElementType,str], do_outtype_sets:bool,
                              dict_elkeyword_add_cschhar:Dict[str, Dict[str, float]]=None)->Dict[int, OutputElementType]:
    """
    Generate new GeoOOFEM input file
    If file already exist, it will be overwritten

    :param file_path: file path to the new file
    :param mdl: fem model
    :param eltypes: dictionary where key is element type, value is oofem element keyword for specified element type
    :param do_outtype_sets: do additional sets containing all elements of each type
    :return: Dictionary of generated sets for each element type, key is element type, value is set id
    """
    null = None
    oofemin = OOFEMInputFile('')

    # handle additonal stifness characteristics, if not set in argument list
    # dict_elkeyword_add_cschhar:Dict[FEMElementType,Dict[str,float]]
    # dictionary (key=FEMElementType) holding dictionary of values (value) for each additional cs characteristics (key)
    if dict_elkeyword_add_cschhar is null:
        dict_elkeyword_add_cschhar:Dict[str, Dict[str, float]] = {}

    for key, value in eltypes.items():
        dict_add_chars =dict_elkeyword_add_cschhar.get(value)
        if dict_add_chars is null:
            dict_add_chars = {}
            dict_elkeyword_add_cschhar[value]=dict_add_chars
        addcschar = get_OOFEMelement_additional_cross_section_characteristics(value)
        for cschar in addcschar:
            if dict_add_chars.get(cschar) is null:
                dict_add_chars[cschar] = get_additional_cschar_defaults(cschar)

    # dictionary holding idcs to FEMElementType , only initialized here
    dict_idcs_elkeyword: Dict[int, str] = {}

    nelem =mdl.num_elements
    ncrosssect =mdl.num_properties
    nmat =mdl.num_materials
    nltf = 1 # Constant
    if mdl.num_lc > 1:
        nltf += mdl.num_lc
    # prepare sets
    dict_idset_outtypes:Dict[OutputElementType,int] = {} # output
    d_outtypes_idset_members: Dict[OutputElementType, List[int]] = {}
    idgroup_to_idset:Dict[int,int] = {}
    nset =0
    # set id in oofem is treated as index, not key!!!
    for id in sorted(mdl.groups):
        nset += 1
        idgroup_to_idset[id]=nset

    if do_outtype_sets:
        for key in eltypes:
            d_outtypes_idset_members[get_element_output_type(key,eltypes[key])]=[]
        for id in sorted(mdl.elements):
            el = mdl.getElement(id)
            typekey = el.get_type()  # possible typekeys can be seen in function get_initial_element_types_dict()
            outkey= get_element_output_type(typekey,eltypes[typekey])
            d_outtypes_idset_members[outkey].append(id)
        for key,value in d_outtypes_idset_members.items():
            if len(value)> 0:
                nset += 1
                dict_idset_outtypes[nset] = key
    nset += 1
    id_all_el_set = nset
    nic =0
    # determine number of loads and boundary conditions
    nlbc= mdl.num_bc # boundary conditions present on all load cases
    for id in mdl.loadcases:
        lc = mdl.getLoadCase(id)
        nlbc += lc.num_bcs + lc.num_loads
    #determine additional nodes for beams
    id_next = max(mdl.nodes.keys()) + 1
    nodes = mdl.nodes.copy()
    beam_nodes = {}
    for id in sorted(mdl.elements):
        beam = mdl.getElement(id)
        if(isinstance(beam,BeamElement)):
            bnodes = []
            web_vec = beam.z_vec
            offset = beam.neutral_axis_offset_z
            if offset != 0.0:
                for node in beam.nodes:
                    rigid_arm_pos = node.p + web_vec * offset
                    ra_nod = NodeRigidArm(id_next, rigid_arm_pos, node)
                    nodes[id_next] = ra_nod
                    id_next += 1
                    bnodes.append(ra_nod)

            else:
                for node in beam.nodes:
                    bnodes.append(node)
            ref_node = null
            if isinstance(beam.orientation, BeamOrientationNode):
                ref_node = beam.orientation.ref_node
            if beam.have_neutral_axis_offset or ref_node is None:
                ref_node_pos =get_referent_node_position(bnodes[0].p,bnodes[1].p,web_vec)
                ref_node=Node(id_next,ref_node_pos)
                nodes[id_next] = ref_node
                id_next += 1
            bnodes.append(ref_node)
            beam_nodes[beam.id] = bnodes
    ndofman = len(nodes)
    data = (ndofman, nelem, ncrosssect, nmat, nlbc, nltf, nset, nic)
    oofemin.init_input_records(data)
    #prepare input records
    # prepare  analysis record
    data = (mdl.num_lc)
    rg = oofemin.getRecordGroup(AnalysisRecord.getkeyname())
    rg.setRecordDataNeutralFormat(0, data)
    #prepare  node records
    masteridofs = {}
    for node in nodes.values():
        if isinstance(node, NodeRigidArm):
            if node.master.id not in masteridofs:
                masteridofs[node.master.id] = set()
            for idof in node.dofidmask:
                masteridofs[node.master.id].add(idof)

    i=0
    for id in sorted(nodes):
        node= nodes.get(id)
        xyz_coords = node.p # 1D array with xyz coords
        if isinstance(node,NodeRigidArm):
            data = ('rigidarmnode', id, xyz_coords, node.master.id, node.dofidmask) # dofidmask = [1, 2, 3, 4, 5, 6]
        else:
            if id in masteridofs:
                data = ('node', id, xyz_coords,masteridofs[id])
            else:
                data = ('node', id, xyz_coords)
        rg = oofemin.getRecordGroup(NodeRecords.getkeyname())
        rg.setRecordDataNeutralFormat(i,data)
        i+=1
    # prepare  element records
    nodeids=[]
    i = 0
    for id in sorted(mdl.elements):
        el = mdl.getElement(id)
        typekey = el.get_type() #this is my Enum from geofem
        oofem_el_keyword = eltypes.get(typekey)
        if oofem_el_keyword is null:
            print(str(typekey) + ' unhandled for element id = '+str(id))
            continue
        dict_idcs_elkeyword[el.prop_id] = oofem_el_keyword # necessary for gethering additional data for properties (crosssections)
        nodeids.clear()
        if id in beam_nodes:
            el_nodes=beam_nodes.get(id)
        else:
            el_nodes=el.nodes
        for node in el_nodes:
            nodeids.append(node.id)
        data = (oofem_el_keyword, id, nodeids, el.prop_id)
        rg = oofemin.getRecordGroup(ElementRecords.getkeyname())
        rg.setRecordDataNeutralFormat(i, data)
        i += 1
    #prepare materials
    i = 0
    for id in sorted(mdl.materials):
        mat = mdl.getMaterial(id)
        data = ('isole',id,mat.E,mat.ni,mat.rho)
        rg = oofemin.getRecordGroup(MaterialRecords.getkeyname())
        rg.setRecordDataNeutralFormat(i, data)
        i += 1
    #prepare properties
    i = 0
    for id in sorted(mdl.properties):
        prop = mdl.getProperty(id)
        idmat=prop.material.id
        if id not in dict_idcs_elkeyword:
            continue # unused property
        add_cschar = dict_elkeyword_add_cschhar[dict_idcs_elkeyword[id]]
        if isinstance(prop,PlateProperty):
            if len(add_cschar) > 0:
                #       proptype, id, idmat, tp, idset,add_cschar
                data = ('plate_add', id, idmat, prop.tp, null,add_cschar)
            else:
                #       proptype, id, idmat, tp, idset
                data = ('plate',id,idmat,prop.tp,null)
        elif isinstance(prop,BeamProperty):
            #proptype,id,idmat,area,Iy,Iz,Ik,idset
            #sfiff_list_names = [' area ', ' Iy ', ' Iz ', ' Ik ', ' beamshearcoeff ', ' shearareay ', ' shearareaz ']
            stiff_list={}
            stiff_list['area']=prop.area
            stiff_list['Iy'] = prop.Iy
            stiff_list['Iz'] = prop.Iz
            stiff_list['Ik'] = prop.Ik
            if prop.is_shear_area_set:
                stiff_list['shearareay'] = prop.shear_area_y
                stiff_list['shearareaz'] = prop.shear_area_z
            else:
                stiff_list['beamshearcoeff'] = prop.shear_coeff
            data = ('k_beam', id, idmat, stiff_list, null)
        elif isinstance(prop,RodProperty):
            data = ('k_rod', id, idmat, prop.get_area(), null)
        else:
            print('Unhandled property found in the model: '+ str(prop))
            i += 1
            continue
        rg = oofemin.getRecordGroup(CrossSectionRecords.getkeyname())
        rg.setRecordDataNeutralFormat(i, data)
        i += 1
    # prepare time functions
    nlc= mdl.num_lc
    lcidtotf = {}
    i=0
    idtf = i + 1
    idLCAll = -999
    lcidtotf[idLCAll]=idtf
    data = (idtf, idLCAll,nlc)
    rg = oofemin.getRecordGroup(TimeFunctionRecords.getkeyname())
    rg.setRecordDataNeutralFormat(i, data)
    if mdl.num_lc > 1:
        ilc = 1
        for id in sorted(mdl.loadcases):
            i += 1
            lc = mdl.getLoadCase(id)
            idtf = i+1
            lcidtotf[id] = idtf
            data = (idtf,ilc,nlc)
            rg = oofemin.getRecordGroup(TimeFunctionRecords.getkeyname())
            rg.setRecordDataNeutralFormat(i, data)
            ilc+=1
    else:
        for id in sorted(mdl.loadcases):
            lcidtotf[id] = idtf

    # prepare boundary conditions and loads
    i=0
    keyword = 'boundarycondition'
    for idbc,bc in mdl.boundaryconditions.items():
        idLCAll = -999
        if isinstance(bc, GroupNodalBC):
            idset = idgroup_to_idset[bc.id_group]
            idfun = lcidtotf[idLCAll]
            dofs = bc.dofs
            vals = bc.values
            data = (keyword, i+1, idset, idfun, dofs, vals)
            rg = oofemin.getRecordGroup(BCandLoadRecords.getkeyname())
            rg.setRecordDataNeutralFormat(i, data)
            i+=1
        else:
            print('Unhandled boundary condition type found in the model: ' + str(bc))

    # for id in sorted(mdl.loadcases):
    #     lc = mdl.getLoadCase(id)
    #     keyword = 'boundarycondition'
    #     for bc in lc.boundaryconditions:
    #         if isinstance(bc, GroupNodalBC):
    #             idset = idgroup_to_idset[bc.id_group]
    #             idfun = lcidtotf[id]
    #             dofs = bc.dofs
    #             vals = bc.values
    #             data = (keyword, i + 1, idset, idfun, dofs, vals)
    #             rg = oofemin.getRecordGroup(BCandLoadRecords.getkeyname())
    #             rg.setRecordDataNeutralFormat(i, data)
    #             i += 1
    #         else:
    #             print('Unhandled boundary condition type found in the model: '+ str(bc))
    for id in sorted(mdl.loadcases):
        lc = mdl.getLoadCase(id)
        for load in lc.loads:
            idload= i + 1 # old way
            #idload = load.id # new way
            if isinstance(load, GroupDoffBasedLoad):
                keyword = 'nodalload'
                idset = idgroup_to_idset[load.id_group]
                idfun = lcidtotf[id]
                dofs = load.dofs
                vals = load.values
                data = (keyword, idload, idset, idfun, dofs, vals)
                rg = oofemin.getRecordGroup(BCandLoadRecords.getkeyname())
                rg.setRecordDataNeutralFormat(i, data)
                i += 1
            elif isinstance(load, AccelerationLoad):
                keyword = 'accelerationload'
                idset = id_all_el_set
                idfun = lcidtotf[id]
                vals = load.get_6dof_values()
                data = (keyword, idload, idset, idfun, vals)
                rg = oofemin.getRecordGroup(BCandLoadRecords.getkeyname())
                rg.setRecordDataNeutralFormat(i, data)
                i += 1
            elif isinstance(load, GroupPressureLoad):
                keyword = 'pressureload'
                idset = idgroup_to_idset[load.id_group]
                idfun = lcidtotf[id]
                press = load.signed_pressure
                data = (keyword, idload, idset, idfun, press)
                rg = oofemin.getRecordGroup(BCandLoadRecords.getkeyname())
                rg.setRecordDataNeutralFormat(i, data)
                i += 1
            else:
                print('Unhandled load type found in the model: '+ str(load))
    # prepare sets (groups)
    iditems = []
    i = 0
    rg = oofemin.getRecordGroup(SetRecords.getkeyname())
    for id in sorted(mdl.groups):
        idset = idgroup_to_idset[id]
        group= mdl.getGroup(id)
        iditems.clear()
        for item in group.items:
            iditems.append(item.id)
        if isinstance(group,NodeGroup):
            data = ('nodes',idset,iditems)
        elif isinstance(group,ElementGroup):
            data = ('elements',idset,iditems)
        else:
            print('Unhandled group type found in the model: ' + str(group))
            continue
        rg.setRecordDataNeutralFormat(i, data)
        i+=1
    #handle  generated sets for element outputs
    for idset,outtype in dict_idset_outtypes.items():
        data = ('elements',idset,d_outtypes_idset_members[outtype])
        rg.setRecordDataNeutralFormat(i, data)
        group = ElementGroup()
        group.init(idset, 'Element_group_' + str(idset))
        idelements = d_outtypes_idset_members[outtype]
        for idel in  idelements:
            el = mdl.getElement(idel)
            group.add_item(el)
        mdl.addGroup(id, group)
        i += 1
    # add the set with all elements at the end
    data = ('allelements', id_all_el_set,[])
    rg.setRecordDataNeutralFormat(i, data)

    #write file from records
    oofemin.write_from_records(file_path)

    return dict_idset_outtypes

    def addOOFEMElementGroup(self, id:int, idelements:List[int]):
        group= ElementGroup()
        group.init(id, 'Element_group_' + str(id))
        for idel in  idelements:
            el = self.getElement(idel)
            group.add_item(el)
        self.addGroup(id,group)

def assign_oofem_analysis_model(mdl:GeoFEM, eltypes:Dict[FEMElementType,str],
                              dict_elkeyword_add_cschhar:Dict[str, Dict[str, float]]=None):
    dir_path = os.path.dirname(mdl.get_input_file_path())
    file_path = os.path.abspath(dir_path) + '\\' + 'oofem_auto_input.in'
    dict_idset_outtypes = generate_OOFEM_input_file(file_path,mdl, eltypes, True,dict_elkeyword_add_cschhar)
    mdl.oofem = OOFEMAnalysisModel(file_path,dict_idset_outtypes)
    pass

def analyse_with_OOFEM(file_path,dict_idset_outtypes, mdl:GeoFEM=None):
    if file_path == '':
        print('Error, file path not set!')
    if not os.path.exists(file_path):
        print('Error, file path does not exist: '+file_path)
        return
    oofem = OOFEMAnalysisModel(file_path)
    #test cs change
    cs = oofem.get_cs(1)
    if False: # Fast test of
        cs.setPropertyValue(csp.CS_Thickness.value,0.15)

    node_out, shell_out, beam_out, react_out, node_id_out= oofem.analyse_all_loadcases(dict_idset_outtypes)
    if mdl is not None:
        shell_group=None
        for id_gr, outtype in dict_idset_outtypes.items():
            if outtype== OutputElementType.Shell:
                shell_group= mdl.getGroup(id_gr)
        if shell_group is not None:
            sxx_top = ElementResult('Shell_CG_Sxx_Top')
            sxx_bottom = ElementResult('Shell_CG_Sxx_Bottom')
            syy_top = ElementResult('Shell_CG_Syy_Top')
            syy_bottom = ElementResult('Shell_CG_Syy_Bottom')
            txy_top = ElementResult('Shell_CG_Tau_xy_Top')
            txy_bottom = ElementResult('Shell_CG_Tau_xy_Bottom')
            mdl.add_element_results(sxx_top)
            mdl.add_element_results(syy_top)
            mdl.add_element_results(txy_top)
            mdl.add_element_results(sxx_bottom)
            mdl.add_element_results(syy_bottom)
            mdl.add_element_results(txy_bottom)
            ilc=0
            for lcID in sorted(mdl.loadcases):
                shell_out_lc = shell_out[ilc]
                sxx_top_lc=sxx_top.add_feres(lcID)
                syy_top_lc = syy_top.add_feres(lcID)
                txy_top_lc = txy_top.add_feres(lcID)
                sxx_bottom_lc = sxx_bottom.add_feres(lcID)
                syy_bottom_lc = syy_bottom.add_feres(lcID)
                txy_bottom_lc = txy_bottom.add_feres(lcID)
                lcshellres = shell_out[ilc]
                iel_res=0
                for ent in shell_group.items:
                    shell_res_el = lcshellres[iel_res]
                    id=ent.id
                    sxx_top_lc[id]=shell_res_el[0]
                    syy_top_lc[id]=shell_res_el[1]
                    txy_top_lc[id]=shell_res_el[2]
                    sxx_bottom_lc[id]=shell_res_el[3]
                    syy_bottom_lc[id]=shell_res_el[4]
                    txy_bottom_lc[id]=shell_res_el[5]
                    iel_res +=1

                ilc+=1
    # print (node_out)
    # print(shell_out)
    # print(beam_out)
    # print(react_out)
    # print(node_id_out)
    return node_out, shell_out, beam_out, react_out, node_id_out
    if mdl is not None:
        id_nodes = [15]
        for id_nod in id_nodes:
            try:
                index_nod = list(mdl.nodes).index(id_nod)
                print('node id = ', id_nod)
                for ilc in range(0, mdl.num_lc):
                    print(node_out[ilc][index_nod])
            except ValueError:
                pass
        id_shell_els = [47, 70]
        for id_set, out_type in dict_idset_outtypes.items():
            if out_type == OutputElementType.Shell:
                group = mdl.getGroup(id_set)

                for id_el in id_shell_els:
                    el=mdl.getElement(id_el)
                    try:
                        index_el=group.items.index(el)
                        print('el id = ',id_el)
                        for ilc in range(0,mdl.num_lc):
                            print(shell_out[ilc][index_el])
                    except ValueError:
                        pass
    return
