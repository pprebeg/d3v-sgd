from enum import Enum
from typing import List, Dict, Set


class FEMElementType(Enum):
    Unknown = 0
    Rod = 1
    Beam = 10
    Tria = 20
    Quad = 30
    StiffQuad = 35

class OutputElementType(Enum):
    Unknown = 0
    Rod = 1
    Beam = 10
    Shell = 20

class CrossSectionProperty(Enum): # Have to be the same as in oofem crossection.h
    CS_Thickness        =400
    CS_Width            =400+1 #< Width
    CS_BeamShearCoeff   =400+2 #< Shear coefficient of beam
    CS_Area             =400+3 #< Area
    CS_InertiaMomentY   =400+4 #< Moment of inertia around y-axis
    CS_InertiaMomentZ   =400+5 #< Moment of inertia around z-axis
    CS_TorsionMomentX   =400+6 #< Moment of inertia around x-axis
    CS_ShearAreaY       =400+7 #< Shear area in y direction
    CS_ShearAreaZ       =400+8 #< Shear area in z direction
    CS_DrillingStiffness=400+9 #< Penalty stiffness for drilling DOFs.
    CS_RelDrillingStiffness=400+10 #< Relative penalty stiffness for drilling DOFs.
    CS_DrillingType     =400+11 #< Type of artificially added drilling stiffness for drilling DOFs.
    CS_TopZCoord        =400+12 #< Top z coordinate
    CS_BottomZCoord     =400+13 #< Bottom z coordinate
    CS_NumLayers        =400+14 #< Number of layers that makes up the cross section
    CS_DirectorVectorX  =400+15 #< Director vector component in x-axis
    CS_DirectorVectorY  =400+16 #< Director vector component in y-axis
    CS_DirectorVectorZ  =400+17 #< Director vector component in z-axis

def get_element_output_type(el_type:FEMElementType, oofem_element_keyname:str='')->OutputElementType:
    #oofem_element_keyname - reserved for future more precise determination
    if el_type is FEMElementType.Quad:
        return OutputElementType.Shell
    if el_type is FEMElementType.Beam:
        return OutputElementType.Beam
    if el_type is FEMElementType.Tria:
        return OutputElementType.Shell
    if el_type is FEMElementType.Rod:
        return OutputElementType.Rod

def get_initial_element_types_dict()->Dict[FEMElementType,str]:
    eltypes  = get_available_oofem_element_types_dict()
    eltypesfirst= {}
    for key,value in eltypes.items():
        eltypesfirst[key]=value[0]
    return eltypesfirst

def get_available_oofem_element_types_dict()->Dict[FEMElementType,List[str]]:
    eltypes  = {
    FEMElementType.Rod:  ["Bar3D"],
    FEMElementType.Beam: ["Beam3D"],
    FEMElementType.Tria: ["tr_shell02"],
    FEMElementType.Quad: ['shellqd42','shellqd41', "mitc4shell"]
    }
    return eltypes


def get_available_cross_section_keywords()->Set[str]:
    cschar= {'thick','area','iy','iz','ik','beamshearcoeff','shearareay','shearareaz','drilltype','reldrillstiffness',
            'material','set'}
    return cschar

def get_OOFEMelement_necessary_cross_section_characteristics(element_keyword:str)->Set[str]:

    if element_keyword == 'shellqd41':
        return {'thick','reldrillstiffness'}
    if element_keyword == 'tr_shell02':
        return {'thick'}
    if element_keyword == 'Beam3D':
        return {'area','iy','beamshearcoeff'}
    if element_keyword == 'mitc4shell':
        return {'thick','drilltype','reldrillstiffness'}
    return {}

def get_OOFEMelement_additional_cross_section_characteristics(element_keyword:str)->Set[str]:
    if element_keyword == 'mitc4shell':
        return {'drilltype','reldrillstiffness'}
    if element_keyword == 'shellqd41':
        return {'reldrillstiffness'}
    return {}

def get_OOFEMelement_supported_cross_section_characteristics_(element_keyword:str)->Set[str]:
    if element_keyword == 'mitc4shell':
        return {'thick','drilltype','reldrillstiffness'}
    if element_keyword == 'tr_shell02':
        return {'thick'}
    if element_keyword == 'Beam3D':
        return {'area','iy','iz','ik','beamshearcoeff','shearareay','shearareaz'}
    return {}

def get_additional_cschar_defaults(cschar_keyword:str):
    if cschar_keyword == 'drilltype':
        return 1
    if cschar_keyword == 'reldrillstiffness':
        return 100.0
    return -999

def get_nip_str_for_element(element_keyword:str):
    if element_keyword == 'tr_shell02':
        return ' nip 4'
    if element_keyword == 'tr_shell01':
        return ' nip 1'
    return ''