import numpy as np
from typing import List, Dict, Set,Tuple
import os
import errno
oofempydir=''
oofempydir='..\\fembin'
if(oofempydir !=''):
    import sys
    sys.path.append(oofempydir)
is_loaded = False
try:
    import fembin.oofempy as oofempy
    is_loaded = True
    print('Used oofempy library: ' + oofempy.__file__)
except ImportError as e:
    print('Warning: oofempy library not found, OOFEM analysis will be disabled! ')
    oofempy = None
from femdir.oofemin import *

class OOFEMAnalysisModel():
    _is_loaded = False
    def __init__(self, fileName,dict_idset_outtypes:Dict[int,OutputElementType]=None):
        if not OOFEMAnalysisModel._is_loaded:
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), 'oofempy')
        null = None
        if fileName is null:
            self._file_path = ''
        self._eng_model = null               # oofempy.EngngModel
        self._vtkxml_nodes = null            # oofempy.VTKXMLExportModule
        self._vtkxml_shell = null            # oofempy.VTKXMLExportModule
        self._vtkxml_beam = null             # oofempy.VTKXMLExportModule
        self._nlc:int = 0
        self._file_path,self._dr,self._eng_model,self._nlc = self.init_problem(fileName)
        self._wdir = os.getcwd()
        self._dict_idset_outtypes:Dict[int,OutputElementType] = dict_idset_outtypes

    @property
    def num_lc(self):
        return self._nlc

    def get_cs(self,id_cs:int):
        domain = self._eng_model.giveDomain(1)
        cs = domain.giveCrossSection(id_cs)
        return cs

    #def init_problem(self,file_path)->(oofempy.EngngModel,oofempy.VTKXMLExportModule,oofempy.VTKXMLExportModule,oofempy.VTKXMLExportModule):
    def init_problem(self, file_path):
        null = None
        try:
            if not os.path.exists(file_path):
                print('Error, file path does not exist: ' + file_path)
                return file_path,null, null,0
            dr = oofempy.OOFEMTXTDataReader(file_path)
            eng_mdl = oofempy.InstanciateProblem(dr, oofempy.problemMode.processor, False, None, False) #: oofempy.EngngModel
            nlc = eng_mdl.giveNumberOfSteps()
        except BaseException as err:
            print(f"Unexpected {err=}, {type(err)=}")
            file_path=None
            dr=None
            eng_mdl = None
            nlc =0
        return file_path, dr, eng_mdl,nlc

    def _check_num_elements_in_set(self,domain,id_set)->int:
        set = domain.giveSet(id_set)
        lst = set.giveElementList()
        #have_elements = lst.isEmpty()
        num_elements = lst.giveSize()
        return num_elements


    def init_analyisis(self, dict_idset_outtypes:Dict[int,OutputElementType]):
        list_shell=[]
        list_beam = []
        domain = self._eng_model.giveDomain(1)
        nset = domain.giveNumberOfSets()  # sets can not have unused ids!
        if dict_idset_outtypes is not None:
            for idset, outtype in dict_idset_outtypes.items():
                if idset > nset:
                    continue
                if idset < 0: # additional option to handle the case where output sets are added at the end of existing sets
                    idset = nset + 1 - idset
                if self._check_num_elements_in_set(domain,idset) == 0:
                    continue
                if outtype is OutputElementType.Shell:
                    list_shell.append(idset)
                elif outtype is OutputElementType.Beam:
                    list_beam.append(idset)

        self._wdir = os.getcwd()
        os.chdir(os.path.dirname(self._file_path))

        self._vtkxml_nodes = oofempy.vtkxml(1, self._eng_model, domain_all=True, tstep_all=True, dofman_all=True,
                                            element_all=True, primvars=(1,), pythonExport=1)

        if len(list_shell) == 0:
            self._vtkxml_shell = None
        else:
            tuple_shell = tuple(list_shell)
            self._vtkxml_shell = oofempy.vtkxml(2, self._eng_model, domain_all=True, tstep_all=True, dofman_all=True,
                                                element_all=False, cellvars=(200,), stype=1,
                                                regionsets=tuple_shell, pythonExport=1)

        if len(list_beam) == 0:
            self._vtkxml_beam = None
        else:
            tuple_beam = tuple(list_beam)
            self._vtkxml_beam = oofempy.vtkxml(3, self._eng_model, domain_all=True, tstep_all=True, dofman_all=True,
                                               element_all=False, cellvars=(201,),
                                               regionsets=tuple_beam, pythonExport=1)

    def finalize_analysis(self):
        os.chdir(self._wdir)

    def analyse_next_loadcase(self):
        null = None
        self._eng_model.preInitializeNextStep()
        self._eng_model.giveNextStep()
        currLC = self._eng_model.giveCurrentStep()
        self._eng_model.initializeYourself(currLC)
        self._eng_model.solveYourselfAt(currLC)
        self._eng_model.updateYourself(currLC)
        self._eng_model.terminate(currLC)
        disp_node_out = self._vtkxml_nodes.getPrimaryVars()['DisplacementVector']
        if self._vtkxml_shell is null:
            shell_el_out = null
        else:
            shell_el_out = self._vtkxml_shell.getCellVars()['IST_Shell_SxSyTxy_Top_Bottom']
        if self._vtkxml_beam is null:
            beam_el_out = null
        else:
            beam_el_out = self._vtkxml_beam.getCellVars()['IST_Beam_FxMyMz_Start_End']
        reactions_dict = self._vtkxml_nodes.getReactionForces()
        reactions = [reactions_dict['nodeid'],
                     reactions_dict['dofid'],
                     reactions_dict['reaction']]
        nodes_def_data=self._vtkxml_nodes.getNodes()['1']
        node_ids = []
        for nod_data in nodes_def_data:
            node_ids.append(nod_data[2])
        return disp_node_out, shell_el_out, beam_el_out,reactions,node_ids


    def analyse_all_loadcases(self, dict_idset_outtypes:Dict[int,OutputElementType]):
        node_out=[]
        shell_out = []
        beam_out = []
        react_out = []
        node_id_out = []
        try:
            self.init_analyisis(dict_idset_outtypes)
            for ilc in range(self.num_lc):
                disp_node_out, shell_el_out, beam_el_out, reactions,node_ids =  self.analyse_next_loadcase()
                node_out.append(disp_node_out)
                shell_out.append(shell_el_out)
                beam_out.append(beam_el_out)
                react_out.append(reactions)
                node_id_out.append(node_ids)
            self.finalize_analysis()
        except BaseException as err:
            print(f"Unexpected {err=}, {type(err)=}")

        return node_out, shell_out, beam_out, react_out, node_id_out



    def analyse(self):
        return self.analyse_all_loadcases(self._dict_idset_outtypes)
OOFEMAnalysisModel._is_loaded = is_loaded

def is_oofempy_loaded():
   return OOFEMAnalysisModel._is_loaded

