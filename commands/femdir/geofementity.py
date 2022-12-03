import time
import openmesh as om
from enum import Enum
import numpy as np
import pathlib
from typing import List, Dict
from femdir.oofemenum import FEMElementType
from femdir.primitives import Icosphere,Cylinder

class ViewType(Enum):
    """!
    ViewType
    """

    constant_color = 0
    face_colors = 1
    face_vertex_colors = 2

class MeshControl():
    def __init__(self):
        self.viewtype=ViewType.constant_color
        self.useviewtreshold=False
        self.uppertreshold=0
        self.lowertreshold = 0
        self.show_beams_with_rigid_arm_as_quad = False
        self.show_nodes = True
        self.show_elements = True
        self.node_complexity = 2
        self.relative_node_radius = 0.005
        self.relative_cylinder_height = 0.05
        self.relative_cylinder_radius = self.relative_node_radius
        self.relative_cylinder_arrow = 2*self.relative_node_radius
        pass
    def getUpperTresholdColor(self):
        return [0.1,0.1,0.1,1.0]
    def getLowerTresholdColor(self):
        return [0.3,0.3,0.3,1.0]
class Property ():
    def __init__(self):
        self.id = 0
        self.name= ''
        pass
    def init(self,id,name):
        self.id = id
        self.name = name
    def get_info(self)->str:
        msg = self.__class__.__name__
        msg +='\n'+self.name+'; id='+str(self.id)
        msg += "\n" + str(self.descriptors)
        return msg

    # optimization property interface
    def is_opt_prop(self):
        return False

    def num_opt_desc(self):
        return 0

    def get_opt_desc_value(self,index):
        return 0

    def set_opt_desc_value(self,index,value):
        pass

    def get_opt_desc_bounds(self,index):
        return 0.0,0.0

    def set_opt_desc_upper_bound(self, index, ub):
        pass

    def set_opt_desc_lower_bound(self, index, lb):
        pass

    @property
    def descriptors (self):
        return None
    @property
    def upperbounds(self):
        return None

    @property
    def lowerbounds(self):
        return  None

class GeoEntity():
    def __init__(self,id:int=0):
        self.id=id
        self.face_value = 0
        self.vertex_based_values = []
        self.face_handles = []
        pass

    def updateMesh(self, mesh: om.PolyMesh, face2entity: Dict[int, int], mc: MeshControl,
                   const_color=[0.4, 1.0, 1.0, 1.0], fun_getcolor=None, minvalue=0, maxvalue=1):
        pass

    def updateDeformedMesh(self, mesh: om.PolyMesh, face2entity: Dict[int, int], mc: MeshControl, deformation,
                           def_scale, const_color=[0.4, 1.0, 1.0, 1.0], fun_getcolor=None, minvalue=0, maxvalue=1):
        pass

    def setFaceValueUsingElementID(self, result: dict):
        self.face_value = result.get(self.id)
        return self.face_value

    def getColorForFaceResult(self, fun_getcolor, minvalue, maxvalue):
        # define color
        if fun_getcolor != None:
            color = fun_getcolor(self.face_value, minvalue, maxvalue)
        else:
            color = []
        return color

    def getColorForFaceVertexResult(self, fun_getcolor, minvalue, maxvalue):
        # define color
        colors = []
        if fun_getcolor != None:
            for i in range(len(self.vertex_based_values)):
                colors.append(fun_getcolor(self.vertex_based_values[i], minvalue, maxvalue))
        return colors

    def init(self,id):
        self.id=id

    def get_info(self)->str:
        msg = self.__class__.__name__
        msg += '\nid=' + str(self.id)
        return msg

class Group(GeoEntity):
    def __init__(self,id:int=0,name:str =''):
        super().__init__(id)
        self._items:List[GeoEntity] = []
        self._name = name
        pass

    def init(self, id, name):
        super().init(id)
        self._name = name

    def extend_items(self,items:List[GeoEntity]):
        self._items.extend(items)

    def add_item(self,item:GeoEntity):
        self._items.append(item)

    def insert_item_to_index(self,index,item:GeoEntity):
            self._items.insert(index,item)

    def insert_item_after_referent_item(self,referent_item:GeoEntity,item:GeoEntity):
        if referent_item in self._items:
            index = self._items.index(referent_item)
            self._items.insert(index+1,item)

    def insert_item_before_referent_item(self,referent_item:GeoEntity,item:GeoEntity):
        if referent_item in self._items:
            index = self._items.index(referent_item)
            self._items.insert(index,item)

    def get_index_of_consecutive_items(self, firstiteme, seconditem):
        for i in range(self.num_items-1):
            if self._items[i] is firstiteme and self._items[i+1] is seconditem:
                return i
        return -1

    def clear(self):
        self._items.clear()

    @property
    def items (self):
        return self._items

    @property
    def name(self):
        return self._name

    @property
    def num_items(self):
        return len(self._items)

class NodeGroup(Group):
    def __init__(self,id:int=0,name:str =''):
        super().__init__(id,name)

class ElementGroup(Group):
    def __init__(self,id:int=0,name:str =''):
        super().__init__(id,name)

class Material (Property):
    def __init__(self):
        super().__init__()
        self._values = [0.0]*4
        pass

    @property
    def ReH(self):
        return self._values[3]

    @ReH.setter
    def ReH(self, value):
        self._values[3] = value

    @property
    def E(self):
        return self._values[0]

    @E.setter
    def E(self, value):
        self._values[0] = value

    @property
    def ni(self):
        return self._values[1]

    @ni.setter
    def ni(self, value):
        self._values[1] = value

    @property
    def rho(self):
        return self._values[2]

    @rho.setter
    def rho(self, value):
        self._values[2] = value

    def init(self, id, name):
        super().init(id,name)

class RodProperty (Property):
    def __init__(self):
        super().__init__()
        self.area = 0.0
        self.diameter = 0.0
        pass

    def get_area(self):
        if self.area > 0.0:
            return  self.area
        else:
            return np.pi*self.diameter**2*0.25

    def init(self, id, name):
        super().init(id,name)

class BeamProperty (Property):
    """
    Beam Property
    GeoOOFEM Beam3d coordinate system definition is used
    y axix is in reference plane defined by reference node and two beam definition nodes
    z axis is perpendicular to reference plane
    z axis coincides with web direction for beams used in ship structures
    """
    def __init__(self,material:Material = None):
        super().__init__()
        self.material=material
        pass
    def init(self, id, name):
        super().init(id,name)
    #region Stifness properties
    @property
    def area(self):
        return  0.0
    @property
    def Iy(self):
        return  0.0
    @property
    def Iz(self):
        return  0.0
    @property
    def Ik(self):
        return 0.0
    @property
    def shear_coeff(self):
        return 0.0
    @property
    def shear_area_y(self):
        return 0.0
    @property
    def shear_area_z(self):
        return 0.0
    @property
    def z_na(self):
        return 0.0
    @property
    def y_na(self):
        return 0.0
    @property
    def is_shear_area_set(self):
        return (self.shear_area_y+self.shear_area_z) > 0.0
    #endregion
    def get_yz_point_pairs_for_plane_visualization(self)->List[np.ndarray]:
        pp = []
        h= np.sqrt(self.Iy*12.0/self.area)
        zna =self.z_na
        pp.append(np.array([0, -zna]))
        pp.append(np.array([0, -zna+h]))
        return pp

class DescriptorBeamProperty (BeamProperty):
    def __init__(self,values:List[float]=None,material:Material = None):
        super().__init__(material)
        null = None
        if values is None:
            self._descriptors = [0.0] * self.get_num_vals()
        else:
            self._descriptors = values.copy()
        self.secType=''
        self._lowerbounds = null
        self._upperbounds = null

    def get_desc_names(self):
        return None

    @property
    def descriptors(self):
        return self._descriptors

    @property
    def upperbounds(self):
        return self._upperbounds

    @property
    def lowerbounds(self):
        return  self._lowerbounds


    def get_num_vals(self):
        return 0

    def set_upperbounds(self,bounds:List[float]):
        self._upperbounds = bounds.copy()

    def set_lowerbounds(self, bounds:List[float]):
        self._lowerbounds = bounds.copy()

    # optimization property interface
    def is_opt_prop(self):
        return (self._lowerbounds is not None) and (self._upperbounds is not None)

    def num_opt_desc(self):
        if self.is_opt_prop():
            return self.get_num_vals()
        return 0

    def get_opt_desc_value(self, index):
        return 0

    def set_opt_desc_value(self, index, value):
        pass

    def get_opt_desc_bounds(self, index):
        return 0.0, 0.0

    def set_opt_desc_upper_bound(self, index, ub):
        if self._upperbounds is None:
            self._upperbounds = [0.0]*self.get_num_vals()
        self._upperbounds[index] = ub

    def set_opt_desc_lower_bound(self, index, lb):
        if self._lowerbounds is None:
            self._lowerbounds = [0.0]*self.get_num_vals()
        self._lowerbounds[index] = lb


class T_Profile_BeamProperty (DescriptorBeamProperty):
    def __init__(self, values: List[float] = None, material: Material = None):
        super().__init__(values, material)

    def get_num_vals(self):
        return 4

    def get_desc_names(self):
        return ['hw', 'tw', 'bf', 'tf']

    @property
    def hw(self):
        return self._descriptors[0]

    @hw.setter
    def hw(self, value):
        self._descriptors[0] = value

    @property
    def tw(self):
        return self._descriptors[1]

    @tw.setter
    def tw(self, value):
        self._descriptors[1] = value

    @property
    def bf(self):
        return self._descriptors[2]

    @bf.setter
    def bf(self, value):
        self._descriptors[2] = value

    @property
    def tf(self):
        return self._descriptors[3]

    @tf.setter
    def tf(self, value):
        self._descriptors[3] = value

    def init(self, id, name):
        super().init(id, name)

    @property
    def z_na(self):
        hw = self.hw
        tf = self.tf
        aw = self.aw
        af = self.af
        area = self.area
        z_na = ((af * (hw + tf * 0.5)) + (aw * (hw * 0.5))) / area
        return z_na

    @property
    def y_na(self):
        return 0.0

    @property
    def aw(self):
        return self.hw * self.tw

    @property
    def af(self):
        return self.bf * self.tf

    @property
    def area(self):
        return self.aw + self.af

    @property
    def Iy(self):
        hw = self.hw
        tf = self.tf
        aw = self.aw
        af = self.af
        zna = self.z_na
        iy = (af * ((tf ** 2) / 12 + (tf * 0.5 + hw - zna) ** 2))
        iy += aw * ((hw ** 2) / 12 + (hw * 0.5 - zna) ** 2)
        return iy

    @property
    def Iz(self):
        return 0.0

    @property
    def Ik(self):
        return 0.0

    @property
    def shear_coeff(self):
        return 0.0

    @property
    def shear_area_y(self):
        return 0.0

    @property
    def shear_area_z(self):
        return 0.0

    def get_yz_point_pairs_for_plane_visualization(self) -> List[np.ndarray]:
        pp = []
        hl = - self.z_na
        hu = hl + self.hw
        pp.append(np.array([0, hl]))
        pp.append(np.array([0, hu]))
        pp.append(np.array([- self.bf / 2.0, hu]))
        pp.append(np.array([self.bf / 2.0, hu]))
        return pp

    def get_info(self)->str:
        msg= super().get_info()
        msg+='\n'+'još nešto'
        return msg


class Half_T_Profile_BeamProperty (T_Profile_BeamProperty):
    def __init__(self, values: List[float] = None, material: Material = None):
        super().__init__(values, material)

    @property
    def area(self):
        return (self.aw + self.af) / 2

    @property
    def z_na(self):
        hw = self.hw
        tf = self.tf
        aw = self.aw / 2
        af = self.af / 2
        area = self.area
        z_na = ((af * (hw + tf * 0.5)) + (aw * (hw * 0.5))) / area
        return z_na

    @property
    def Iy(self):
        hw = self.hw
        tf = self.tf
        aw = self.aw / 2
        af = self.af / 2
        zna = self.z_na
        iy = (af * ((tf ** 2) / 12 + (tf * 0.5 + hw - zna) ** 2))
        iy += aw * ((hw ** 2) / 12 + (hw * 0.5 - zna) ** 2)
        return iy

    def get_yz_point_pairs_for_plane_visualization(self) -> List[np.ndarray]:
        pp = []
        hl = - self.z_na
        hu = hl + self.hw
        pp.append(np.array([0, hl]))
        pp.append(np.array([0, hu]))
        pp.append(np.array([- self.bf / 2, hu]))
        pp.append(np.array([0, hu]))
        return pp


class L_Profile_BeamProperty (T_Profile_BeamProperty):
    def __init__(self, values: List[float] = None, material: Material = None):
        super().__init__(values, material)

    def get_yz_point_pairs_for_plane_visualization(self) -> List[np.ndarray]:
        pp = []
        hl = - self.z_na
        hu = hl + self.hw
        pp.append(np.array([0, hl]))
        pp.append(np.array([0, hu]))
        pp.append(np.array([- self.bf, hu]))
        pp.append(np.array([0, hu]))
        return pp

    def get_info(self)->str:
        msg= super().get_info()
        # msg+='\n'+'još nešto'
        return msg


class Half_L_Profile_BeamProperty (L_Profile_BeamProperty):
    def __init__(self, values: List[float] = None, material: Material = None):
        super().__init__(values, material)

    @property
    def area(self):
        return (self.aw + self.af) / 2

    @property
    def z_na(self):
        hw = self.hw
        tf = self.tf
        aw = self.aw / 2
        af = self.af / 2
        area = self.area
        z_na = ((af * (hw + tf * 0.5)) + (aw * (hw * 0.5))) / area
        return z_na

    @property
    def Iy(self):
        hw = self.hw
        tf = self.tf
        aw = self.aw / 2
        af = self.af / 2
        zna = self.z_na
        iy = (af * ((tf ** 2) / 12 + (tf * 0.5 + hw - zna) ** 2))
        iy += aw * ((hw ** 2) / 12 + (hw * 0.5 - zna) ** 2)
        return iy

    def get_yz_point_pairs_for_plane_visualization(self) -> List[np.ndarray]:
        pp = []
        hl = - self.z_na
        hu = hl + self.hw
        pp.append(np.array([0, hl]))
        pp.append(np.array([0, hu]))
        pp.append(np.array([- self.bf / 2, hu]))
        pp.append(np.array([0, hu]))
        return pp


class Hat_Profile_BeamProperty (DescriptorBeamProperty):
    def __init__(self, values: List[float] = None, material: Material = None):
        super().__init__(values, material)

    def get_num_vals(self):
        return 4

    def get_desc_names(self):
        return ['h', 't', 'bf', 'fi']

    @property
    def h(self):
        return self._descriptors[0]

    @h.setter
    def h(self, value):
        self._descriptors[0] = value

    @property
    def t(self):
        return self._descriptors[1]

    @t.setter
    def t(self, value):
        self._descriptors[1] = value

    @property
    def bf(self):
        return self._descriptors[2]

    @bf.setter
    def bf(self, value):
        self._descriptors[2] = value

    @property
    def fi(self):
        return self._descriptors[3]

    @fi.setter
    def fi(self, value):
        self._descriptors[3] = value

    def init(self, id, name):
        super().init(id, name)

    @property
    def tan_fi(self):
        fi = self.fi
        return np.tan(np.radians(fi))

    @property
    def sin_fi(self):
        fi = self.fi
        return np.sin(np.radians(fi))

    @property
    def hat_s1(self):
        h = self.h
        t = self.t
        bf = self.bf
        fi = self.fi
        tan_fi = self.tan_fi
        s1 = bf + ((2 * h + t) / tan_fi) - (t * np.tan(np.radians(fi / 2)))
        return s1

    @property
    def area(self):
        h = self.h
        t = self.t
        bf = self.bf
        tan_fi = self.tan_fi
        sin_fi = self.sin_fi
        area = ((2 * h * t) / sin_fi) + (t ** 2 / tan_fi) + (bf * t)
        return area

    @property
    def z_na(self):
        h = self.h
        t = self.t
        tan_fi = self.tan_fi
        sin_fi = self.sin_fi
        area = self.area
        z_na = ((h * t / sin_fi) * (h + t) + t ** 3 /
                (6 * tan_fi)) / area
        return z_na

    @property
    def y_na(self):
        return 0.0

    @property
    def Iy(self):
        h = self.h
        t = self.t
        bf = self.bf
        tan_fi = self.tan_fi
        sin_fi = self.sin_fi
        zna = self.z_na
        Iy = 2 * (t * h ** 3 / (12 * sin_fi))
        Iy += 2 * (((t + h) / 2 - zna) ** 2) * (h * t / sin_fi)
        Iy += 2 * (t ** 4 / (36 * tan_fi))
        Iy += 2 * (((zna - t / 6) ** 2) * (t ** 2 / (2 * tan_fi)))
        Iy += ((bf * t ** 3 / 12) + zna ** 2 * bf * t)
        return Iy

    @property
    def Iz(self):
        return 0.0

    @property
    def Ik(self):
        return 0.0

    @property
    def shear_coeff(self):
        return 0.0

    @property
    def shear_area_y(self):
        return 0.0

    @property
    def shear_area_z(self):
        return 0.0

    def get_yz_point_pairs_for_plane_visualization(self) -> List[np.ndarray]:
        # hl is set 5mm below plating to not show web lines on the other side
        pp = []
        hl = - self.z_na + 5
        hu = hl + self.h
        s1 = self.hat_s1
        bf = self.bf
        pp.append(np.array([- s1 / 2, hl]))
        pp.append(np.array([- bf / 2, hu]))
        pp.append(np.array([- bf / 2, hu]))
        pp.append(np.array([bf / 2, hu]))
        pp.append(np.array([bf / 2, hu]))
        pp.append(np.array([s1 / 2, hl]))
        return pp

    def get_info(self) -> str:
        msg= super().get_info()
        # msg+='\n'+'još nešto'
        return msg


class Half_Hat_Profile_BeamProperty (Hat_Profile_BeamProperty):
    def __init__(self, values: List[float] = None, material: Material = None):
        super().__init__(values, material)

    @property
    def area(self):
        h = self.h
        t = self.t
        bf = self.bf
        tan_fi = self.tan_fi
        sin_fi = self.sin_fi
        area = ((2 * h * t) / sin_fi) + (t ** 2 / tan_fi) + (bf * t)
        return area / 2

    @property
    def z_na(self):
        h = self.h
        t = self.t
        tan_fi = self.tan_fi
        sin_fi = self.sin_fi
        area = self.area
        z_na = ((h * t / sin_fi) * (h + t) + t ** 3 /
                (6 * tan_fi)) / area
        return z_na / 2

    @property
    def Iy(self):
        h = self.h
        t = self.t
        bf = self.bf
        tan_fi = self.tan_fi
        sin_fi = self.sin_fi
        zna = self.z_na
        Iy = t * h ** 3 / (12 * sin_fi)
        Iy += (((t + h) / 2 - zna) ** 2) * (h * t / sin_fi)
        Iy += t ** 4 / (36 * tan_fi)
        Iy += ((zna - t / 6) ** 2) * (t ** 2 / (2 * tan_fi))
        Iy += ((bf / 2) * t ** 3 / 12) + zna ** 2 * (bf / 2) * t
        return Iy

    def get_yz_point_pairs_for_plane_visualization(self) -> List[np.ndarray]:
        # hl is set 5mm below plating to not show web lines on the other side
        pp = []
        hl = - self.z_na + 5
        hu = hl + self.h
        s1 = self.hat_s1
        bf = self.bf
        pp.append(np.array([- s1 / 2, hl]))
        pp.append(np.array([- bf / 2, hu]))
        pp.append(np.array([- bf / 2, hu]))
        pp.append(np.array([0, hu]))
        return pp


class Bulb_Profile_BeamProperty (DescriptorBeamProperty):
    def __init__(self, values: List[float] = None, material: Material = None):
        super().__init__(values, material)

    def get_num_vals(self):
        return 4

    def get_desc_names(self):
        return ['hw_ekv', 'tw_ekv', 'bf_ekv', 'tf_ekv']

    @property
    def hw_ekv(self):
        return self._descriptors[0]

    @hw_ekv.setter
    def hw_ekv(self, value):
        self._descriptors[0] = value

    @property
    def tw_ekv(self):
        return self._descriptors[1]

    @tw_ekv.setter
    def tw_ekv(self, value):
        self._descriptors[1] = value

    @property
    def bf_ekv(self):
        return self._descriptors[2]

    @bf_ekv.setter
    def bf_ekv(self, value):
        self._descriptors[2] = value

    @property
    def tf_ekv(self):
        return self._descriptors[3]

    @tf_ekv.setter
    def tf_ekv(self, value):
        self._descriptors[3] = value

    def init(self, id, name):
        super().init(id, name)

    @property
    def z_na(self):
        hw = self.hw_ekv
        tf = self.tf_ekv
        aw = self.aw
        af = self.af
        area = self.area
        z_na = ((af * (hw + tf * 0.5)) + (aw * (hw * 0.5))) / area
        return z_na

    @property
    def y_na(self):
        return 0.0

    @property
    def aw(self):
        return self.hw_ekv * self.tw_ekv

    @property
    def af(self):
        return self.bf_ekv * self.tf_ekv

    @property
    def area(self):
        return self.aw + self.af

    @property
    def Iy(self):
        hw = self.hw_ekv
        tf = self.tf_ekv
        aw = self.aw
        af = self.af
        zna = self.z_na
        iy = (af * ((tf ** 2) / 12 + (tf * 0.5 + hw - zna) ** 2))
        iy += aw * ((hw ** 2) / 12 + (hw * 0.5 - zna) ** 2)
        return iy

    @property
    def Iz(self):
        return 0.0

    @property
    def Ik(self):
        return 0.0

    @property
    def shear_coeff(self):
        return 0.0

    @property
    def shear_area_y(self):
        return 0.0

    @property
    def shear_area_z(self):
        return 0.0

    def get_yz_point_pairs_for_plane_visualization(self) -> List[np.ndarray]:
        pp = []
        hl = - self.z_na
        hu = hl + self.hw_ekv
        pp.append(np.array([0, hl]))
        pp.append(np.array([0, hu]))
        pp.append(np.array([- self.bf_ekv, hu]))
        pp.append(np.array([0, hu]))
        return pp


class Half_Bulb_Profile_BeamProperty (Bulb_Profile_BeamProperty):
    def __init__(self, values: List[float] = None, material: Material = None):
        super().__init__(values, material)

    @property
    def area(self):
        return (self.aw + self.af) / 2

    @property
    def z_na(self):
        hw = self.hw_ekv
        tf = self.tf_ekv
        aw = self.aw / 2
        af = self.af / 2
        area = self.area
        z_na = ((af * (hw + tf * 0.5)) + (aw * (hw * 0.5))) / area
        return z_na

    @property
    def Iy(self):
        hw = self.hw_ekv
        tf = self.tf_ekv
        aw = self.aw / 2
        af = self.af / 2
        zna = self.z_na
        iy = (af * ((tf ** 2) / 12 + (tf * 0.5 + hw - zna) ** 2))
        iy += aw * ((hw ** 2) / 12 + (hw * 0.5 - zna) ** 2)
        return iy

    def get_yz_point_pairs_for_plane_visualization(self) -> List[np.ndarray]:
        pp = []
        hl = - self.z_na
        hu = hl + self.hw_ekv
        pp.append(np.array([0, hl]))
        pp.append(np.array([0, hu]))
        pp.append(np.array([- self.bf_ekv / 2, hu]))
        pp.append(np.array([0, hu]))
        return pp


class FB_Profile_BeamProperty (DescriptorBeamProperty):
    def __init__(self, values: List[float] = None, material: Material = None):
        super().__init__(values, material)

    def get_num_vals(self):
        return 2

    def get_desc_names(self):
        return ['hw', 'tw']

    @property
    def hw(self):
        return self._descriptors[0]

    @hw.setter
    def hw(self, value):
        self._descriptors[0] = value

    @property
    def tw(self):
        return self._descriptors[1]

    @tw.setter
    def tw(self, value):
        self._descriptors[1] = value

    def init(self, id, name):
        super().init(id, name)

    @property
    def z_na(self):
        hw = self.hw
        return hw * 0.5

    @property
    def y_na(self):
        return 0.0

    @property
    def area(self):
        return self.hw * self.tw

    @property
    def Iy(self):
        hw = self.hw
        tw = self.tw
        return (tw * hw ** 3) / 12

    @property
    def Iz(self):
        return 0.0

    @property
    def Ik(self):
        return 0.0

    @property
    def shear_coeff(self):
        return 0.0

    @property
    def shear_area_y(self):
        return 0.0

    @property
    def shear_area_z(self):
        return 0.0

    def get_yz_point_pairs_for_plane_visualization(self) -> List[np.ndarray]:
        pp = []
        hl = - self.z_na
        hu = hl + self.hw
        pp.append(np.array([0, hl]))
        pp.append(np.array([0, hu]))
        return pp


class Half_FB_Profile_BeamProperty (FB_Profile_BeamProperty):
    def __init__(self, values: List[float] = None, material: Material = None):
        super().__init__(values, material)

    @property
    def area(self):
        return self.hw * self.tw / 2

    @property
    def Iy(self):
        hw = self.hw
        tw = self.tw / 2
        return (tw * hw ** 3) / 12


class T_ProfileAttachPlate_BeamProperty (T_Profile_BeamProperty):
    def __init__(self,values:List[float]=None,material:Material = None):
        super().__init__(values,material)


    def get_num_vals(self):
        return 6

    def get_desc_names(self):
        return ['hw','tw','bf','tf','bp','tp']

    def init(self, id, name):
        super().init(id,name)

    @property
    def bp(self):
        return self._descriptors[4]
    @bp.setter
    def bp(self, value):
        self._descriptors[4] = value

    @property
    def tp(self):
        return self._descriptors[5]

    @tp.setter
    def tp(self, value):
        self._descriptors[5] = value

    @property
    def area(self):
        return self.aw + self.af + self.ap
    @property
    def Iy(self):
        yna = self.y_na
        aw= self.aw
        af = self.af
        ap = self.ap
        dw = yna - (self.hw+self.tp)/2.0
        df = yna - self.hw+self.tf/2.0
        dp = yna - self.tp / 2.0
        iy = (self.hw**2.0*aw) / 12.0 + self.aw * dw ** 2.0
        iy += (self.tf**2.0*af) / 12.0 + self.af * df ** 2.0
        iy += (self.tp ** 2.0 * ap) / 12.0 + self.ap * dp ** 2.0
        return iy
    @property
    def Iz(self):
        return  0.0
    @property
    def Ik(self):
        return 0.0
    @property
    def shear_coeff(self):
        return 0.0
    @property
    def shear_area_y(self):
        return 0.0
    @property
    def shear_area_z(self):
        return 0.0
    @property
    def z_na(self):
        return (self.aw * (self.hw + self.tp) / 2.0 + self.af * (+self.tp / 2.0 + self.hw + self.bf / 2.0)) / self.area
    @property
    def ap(self):
        return self.bp * self.tp

    def get_yz_point_pairs_for_plane_visualization(self)->List[np.ndarray]:
        pp = []
        hl=-self.z_na
        hu=hl+self.hw
        pp.append(np.array([-self.bp / 2.0, hl]))
        pp.append(np.array([self.bp / 2.0, hl]))
        pp.append(np.array([0,hl]))
        pp.append(np.array([0,hu]))
        pp.append(np.array([-self.bf / 2.0, hu]))
        pp.append(np.array([self.bf / 2.0, hu]))
        return pp


class StifnessBeamProperty (BeamProperty):
    def __init__(self,material:Material = None):
        super().__init__(material)
        self._dstifnesses:Dict[str,float] = {}

    def add_stiffness_characteristic(self,key:str,value:float):
        self._dstifnesses[key]=value
    def init(self, id, name):
        super().init(id,name)

    def _get_stiff_value(self,key):
        value= self._dstifnesses.get(key)
        if value is None:
            return 0.0
        else:
            return value
    @property
    def area(self):
        return  self._get_stiff_value('area')
    @property
    def Iy(self):
        return  self._get_stiff_value('iy')
    @property
    def Iz(self):
        return  self._get_stiff_value('iz')
    @property
    def Ik(self):
        return  self._get_stiff_value('ik')
    @property
    def shear_coeff(self):
        return  self._get_stiff_value('beamshearcoeff')
    @property
    def shear_area_y(self):
        return  self._get_stiff_value('shearareay')
    @property
    def shear_area_z(self):
        return  self._get_stiff_value('shearareaz')
    @property
    def z_na(self):
        return 0.0
    @property
    def y_na(self):
        return 0.0
    def get_yz_point_pairs_for_plane_visualization(self)->List[np.ndarray]:
        pp = []
        h= np.sqrt(self.Iy*12.0/self.area)
        zna = h/2
        pp.append(np.array([0, -zna]))
        pp.append(np.array([0, -zna+h]))
        return pp

class StiffLayoutProperty (Property):
    def __init__(self):
        self.beam = 0
        pass
    def init(self, id, name):
        super().init(id,name)
class PlateProperty (Property):
    def __init__(self):
        self._tp=0
        self.material =0
        pass
    def init(self, id, name):
        super().init(id,name)


    @property
    def tp(self):
        return self._tp

    @tp.setter
    def tp(self, value):
        self._tp = value

    @property
    def descriptors (self):
        return [self._tp]




class Node(GeoEntity):
    def __init__(self,id:int=0,p:np.ndarray = None):
        super().__init__(id)
        if p is None:
            self.p =np.zeros((3))
        else:
            self.p = p
        pass

    _sphere:Icosphere = None

    @property
    def sphere(self)->Icosphere:
        return Node._sphere

    @staticmethod
    def set_sphere(referent_dimension:float,mc:MeshControl):
        radius =mc.relative_node_radius*referent_dimension
        Node._sphere = Icosphere(radius,mc.node_complexity-1)

    @staticmethod
    def get_distance_betweeen_nodes(node1, node2):
        p1 = node1.p
        p2 = node2.p
        p3= p2-p1
        return np.linalg.norm(p3)

    @staticmethod
    def get_point_betweeen_nodes(node1, node2,relative_position):
        p1 = node1.p
        p2 = node2.p
        p3 = p2 - p1
        l= np.linalg.norm(p3)
        p3=p3/l
        pbwn= p1+relative_position*l*p3
        return pbwn

    @staticmethod
    def get_area_between_nodes(nodes:List):
        n = len(nodes)
        area=0.0
        if n >= 3:
            for i in range(1,n-1):
                p1p2 = np.subtract(nodes[i], nodes[0])
                p1p3 = np.subtract(nodes[i+1], nodes[0])
                u =  np.cross(p1p2, p1p3)
                area += 0.5* np.linalg.norm(u)
        return area

    def x(self):
        return self.p[0]
    def y(self):
        return self.p[1]
    def z(self):
        return self.p[2]
    def init(self,id,x,y,z):
        super().init(id)
        self.p[0] = x
        self.p[1] = y
        self.p[2] = z

    def updateMesh(self, mesh: om.PolyMesh, face2entity: Dict[int, int], mc: MeshControl,
                   const_color=[1.0, 0.0, 0.0, 1.0], fun_getcolor=None, minvalue=0, maxvalue=1):
        color = const_color
        # define color
        if fun_getcolor != None:
            color = self.getColorForFaceResult(fun_getcolor, minvalue, maxvalue)
        vhandle = []
        handleColorIndex = []

        self.face_handles =self.sphere.append_mesh_at_center(mesh,self.p)
        f2e = dict(zip(self.face_handles,[self]*len(self.face_handles)))
        face2entity.update(f2e)
        colors = mesh.face_colors()
        colors[self.face_handles]=color
        return

    def get_info(self)->str:
        msg = self.__class__.__name__
        msg+='\nid='+str(self.id)
        msg += '\ncoords' + str(self.p)
        return msg



class NodeRigidArm(Node):
    def __init__(self,id:int=0,p:np.ndarray=None,
                 master:Node=None,dofidmask:List[int]=None):
        super().__init__(id,p)
        self._master = master
        if dofidmask is None:
            self._dofidmask = [1, 2, 3, 4, 5, 6]
        else:
            self._dofidmask = dofidmask
        pass
    def init(self,id,x,y,z,master:Node,dofidmask:List[int]=[1,2,3,4,5,6]):
        super().init(id,x,y,z)
        self._master = master
        self._dofidmask = dofidmask


    @property
    def master(self)->Node:
        return self._master

    @property
    def dofidmask(self)->List[int]:
        return self._dofidmask
    def updateMesh(self, mesh: om.PolyMesh, face2entity: Dict[int, int], mc: MeshControl,
                   const_color=[0.0, 1.0, 0.0, 1.0], fun_getcolor=None, minvalue=0, maxvalue=1):
        color = const_color
        # define color
        if fun_getcolor != None:
            color = self.getColorForFaceResult(fun_getcolor, minvalue, maxvalue)
        vhandle = []
        handleColorIndex = []

        self.face_handles =self.sphere.append_mesh_at_center(mesh,self.p)
        f2e = dict(zip(self.face_handles,[self]*len(self.face_handles)))
        face2entity.update(f2e)
        colors = mesh.face_colors()
        colors[self.face_handles]=color
        return
class Element(GeoEntity):
    def __init__(self):
        super().__init__()
        self.property:Property = 0
        self.nodes = []

        self.value_name=""

        pass

    def get_type(self)->FEMElementType:
        return FEMElementType.Unknown

    @property
    def prop_id(self):
        return self.property.id



    def addNode(self,node):
        self.nodes.append(node)
        self.vertex_based_values.append(0)

    def init(self,id):
        super().init(id)
    def updateMesh(self, mesh: om.PolyMesh, face2entity: Dict[int, int], mc: MeshControl,
                   const_color=[0.4, 1.0, 1.0, 1.0], fun_getcolor=None, minvalue=0, maxvalue=1):
        pass

    def updateDeformedMesh(self, mesh: om.PolyMesh, face2entity: Dict[int, int], mc: MeshControl, deformation,
                           def_scale, const_color=[0.4, 1.0, 1.0, 1.0], fun_getcolor=None, minvalue=0, maxvalue=1):
        pass

    def onTPLValue(self):
        self.face_value = 0
        return self.face_value

    def onMatID(self):
        self.face_value = self.property.material.id
        return self.face_value

    def onPropID(self):
        self.face_value = self.property.id
        return self.face_value

    def get_info(self)->str:
        msg = self.__class__.__name__
        msg+='\nid='+str(self.id)
        msg += '\nnodes'
        for node in self.nodes:
            msg += ' ' + str(node.id)
        msg+= '\n'+self.property.get_info()
        return msg

    @property
    def num_nodes(self):
        return len(self.nodes)

    def get_element_mass_ycg(self):
        pass


class RodElement(Element):
    def __init__(self):
        super().__init__()
        pass
    def get_type(self):
        return FEMElementType.Rod

class BeamOrientation():
    def __init__(self):
        pass
    def get_z_axis_unit_vector(self, force:bool, n1:np.ndarray=None, n2:np.ndarray=None)->np.ndarray:
        pass
    def get_referent_node_position(self, force:bool, n1:np.ndarray=None, n2:np.ndarray=None)->np.ndarray:
        pass

class BeamOrientationVector(BeamOrientation) :
    def __init__(self, direction:np.ndarray):
        self._z_axis_unit_vector = direction / np.linalg.norm(direction)

    def set_z_axis_unit_vector(self, direction:np.ndarray):
        self._z_axis_unit_vector = direction / np.linalg.norm(direction)

    def get_z_axis_unit_vector(self, force:bool, n1:np.ndarray=None, n2:np.ndarray=None)->np.ndarray:
        return self._z_axis_unit_vector

    def get_referent_node_position(self, force:bool, n1:np.ndarray=None, n2:np.ndarray=None) ->np.ndarray:
        if force:
            xvec = n2-n1
            y_vec = - np.cross(xvec,self._z_axis_unit_vector)
            y_vec = y_vec/np.linalg.norm(y_vec)
            ref_vec = (xvec + y_vec * np.linalg.norm(xvec)) / 2.0
            ref_node_pos = n1 + ref_vec
            return ref_node_pos
        return None


class BeamOrientationNode(BeamOrientation):
    def __init__(self,node: Node = None):
        self.set_referent_node(node)

    def set_referent_node(self, node: Node):
        self._ref_node = node

    def get_z_axis_unit_vector(self, force:bool, n1:np.ndarray=None, n2:np.ndarray=None) -> np.ndarray:
        if force:
            xvec = n2-n1
            kvec = self._ref_node.p - n1
            web_vec = np.cross(xvec,kvec)
            return web_vec/np.linalg.norm(web_vec)
        return None
    def get_referent_node_position(self, force:bool, n1:np.ndarray=None, n2:np.ndarray=None) ->np.ndarray:
        return self._ref_node.p
    @property
    def ref_node(self):
        return self._ref_node

class BeamElement(Element):
    """
    Beam Element
    GeoOOFEM Beam3d coordinate system definition is used
    x axis from first node to second node
    y axix is in reference plane defined by reference node and two beam definition nodes
    z axis is perpendicular to reference plane
    z axis coincides with web direction for beams used in ship structures
    """
    def __init__(self):
        super().__init__()
        self._orientation:BeamOrientation = None
        self._length=0.0

    @property
    def length(self):
        if self._length == 0.0:
            if self.num_nodes >=2:
                self._length = Node.get_distance_betweeen_nodes(self.nodes[1],self.nodes[0])
            pass
        return  self._length

    @property
    def beam_prop(self)->BeamProperty:
        return self.property
    @property
    def z_vec(self):
        return self._orientation.get_z_axis_unit_vector(True,self.nodes[0].p, self.nodes[1].p)
    @property
    def have_neutral_axis_offset(self):
        return (self.neutral_axis_offset_z != 0.0 or self.neutral_axis_offset_y != 0.0)
    @property
    def neutral_axis_offset_z(self):
        return 0.0
    @property
    def neutral_axis_offset_y(self):
        return 0.0
    @property
    def orientation(self):
        return self._orientation

    def set_beam_orientation(self,orientation:BeamOrientation):
        self._orientation = orientation

    def get_type(self):
        return FEMElementType.Beam
    def get_element_mass_ycg(self):
        mass=self.beam_prop.area*self.beam_prop.material.rho
        mid_point = Node.get_point_betweeen_nodes(self.nodes[0],self.nodes[1],0.5)
        cg_point = mid_point + self.orientation.get_z_axis_unit_vector(True)*self.beam_prop.z_na
        return mass,cg_point[1]

    def onTPLValue(self):
        self.face_value = self.property.tw
        return self.face_value

    def updateTriMesh(self,mesh:om.TriMesh,mc:MeshControl, const_color = [0.4, 1.0, 1.0, 1.0],
                   fun_getcolor=None,minvalue = 0,maxvalue = 1 ):
        color = const_color
        vhandle = []
        handleColorIndex = []
        fhs = []

        hw=self.property.hw
        tw=self.property.tw
        bf=self.property.bf
        tf = self.property.tf
        x = self.nodes[0].p - self.nodes[1].p
        y = self.z_vec
        v = np.cross(x, y)
        # z = self.nodes[0].p + self.wo*hw - (v*bf*0.5)

        vhandle.append(mesh.add_vertex(self.nodes[0].p)) #0 točka 1
        handleColorIndex.append(0)
        vhandle.append(mesh.add_vertex(self.nodes[1].p)) #1 točka 2
        handleColorIndex.append(1)

        data = self.nodes[0].p + self.z_vec * hw  #2 točka 3
        vhandle.append(mesh.add_vertex(data))
        handleColorIndex.append(0)
        data = self.nodes[1].p + self.z_vec * hw  # 3 točka 4
        vhandle.append(mesh.add_vertex(data))
        handleColorIndex.append(1)
        data = self.nodes[0].p + self.z_vec * hw - (v * bf * 0.5)  #4 točka 5
        vhandle.append(mesh.add_vertex(data))
        handleColorIndex.append(0)
        data = self.nodes[1].p + self.z_vec * hw - (v * bf * 0.5)  #5 točka 6
        vhandle.append(mesh.add_vertex(data))
        handleColorIndex.append(1)
        data = self.nodes[1].p + self.z_vec * hw + (v * bf * 0.5)  #6 točka 7
        vhandle.append(mesh.add_vertex(data))
        handleColorIndex.append(1)
        data = self.nodes[0].p + self.z_vec * hw + (v * bf * 0.5)  #7 točka 8
        vhandle.append(mesh.add_vertex(data))
        handleColorIndex.append(0)

        # data = np.array([self.nodes[0].p]) + np.array(self.wo)*self.hw - np.cross(np.array([self.nodes[1].p]),np.array(self.wo))*(self.bf*0.5) #4 točak 5
        # vhandle.append(mesh.add_vertex(data))

        fhs.append(mesh.add_face(vhandle[0], vhandle[1], vhandle[2])) #1-2-3
        fhs.append(mesh.add_face(vhandle[2], vhandle[1], vhandle[3])) #3-2-4
        fhs.append(mesh.add_face(vhandle[4], vhandle[5], vhandle[6])) #5-6-7
        fhs.append(mesh.add_face(vhandle[4], vhandle[6], vhandle[7])) # 5-7-8

        # define color
        if fun_getcolor != None:
            if mc.viewtype == ViewType.face_colors:
                color = self.getColorForFaceResult(fun_getcolor, minvalue, maxvalue)
                for fh in fhs:
                    mesh.set_color(fh, color)
            if mc.viewtype == ViewType.face_vertex_colors:
                colors = self.getColorForFaceVertexResult(fun_getcolor, minvalue, maxvalue)
                for ivh in range(len(vhandle)):
                    mesh.set_color(vhandle[ivh], color[handleColorIndex[ivh]])

    def updateMesh_old(self,mesh:om.PolyMesh,face2element:Dict[int,int],mc:MeshControl, const_color = [0.4, 1.0, 1.0, 1.0],
                   fun_getcolor=None,minvalue = 0,maxvalue = 1 ):
        color = const_color
        vhandle = []
        handleColorIndex = []

        hw=self.property.hw
        tw=self.property.tw
        bf=self.property.bf
        tf = self.property.tf
        x = self.nodes[0].p - self.nodes[1].p
        y = self.z_vec
        v = np.cross(x, y)
        # z = self.nodes[0].p + self.wo*hw - (v*bf*0.5)

        vhandle.append(mesh.add_vertex(self.nodes[0].p)) #0 točka 1
        handleColorIndex.append(0)
        vhandle.append(mesh.add_vertex(self.nodes[1].p)) #1 točka 2
        handleColorIndex.append(1)

        data = self.nodes[0].p + self.z_vec * hw * 0.99  #2 točka 4
        vhandle.append(mesh.add_vertex(data))
        handleColorIndex.append(0)
        data = self.nodes[1].p + self.z_vec * hw * 0.99  # 3 točka 3
        vhandle.append(mesh.add_vertex(data))
        handleColorIndex.append(1)
        data = self.nodes[0].p + self.z_vec * hw - (v * bf * 0.5)  #4 točka 5
        vhandle.append(mesh.add_vertex(data))
        handleColorIndex.append(0)
        data = self.nodes[1].p + self.z_vec * hw - (v * bf * 0.5)  #5 točka 6
        vhandle.append(mesh.add_vertex(data))
        handleColorIndex.append(1)
        data = self.nodes[1].p + self.z_vec * hw + (v * bf * 0.5)  #6 točka 7
        vhandle.append(mesh.add_vertex(data))
        handleColorIndex.append(1)
        data = self.nodes[0].p + self.z_vec * hw + (v * bf * 0.5)  #7 točka 8
        vhandle.append(mesh.add_vertex(data))
        handleColorIndex.append(0)

        # data = np.array([self.nodes[0].p]) + np.array(self.wo)*self.hw - np.cross(np.array([self.nodes[1].p]),np.array(self.wo))*(self.bf*0.5) #4 točak 5
        # vhandle.append(mesh.add_vertex(data))
        fh = mesh.add_face(vhandle[0], vhandle[1], vhandle[3], vhandle[2]).idx() #1-2-3-4
        self.face_handles.append(fh)
        face2element[fh] = self
        fh = mesh.add_face(vhandle[4], vhandle[5], vhandle[6], vhandle[7]).idx() #5-6-7-8
        self.face_handles.append(fh)
        face2element[fh] = self

        # define color
        if fun_getcolor != None:
            if mc.viewtype == ViewType.face_colors:
                color = self.getColorForFaceResult(fun_getcolor, minvalue, maxvalue)
                for fhidx in self.face_handles:
                    fh = mesh.face_handle(fhidx)
                    mesh.set_color(fh, color)
            if mc.viewtype == ViewType.face_vertex_colors:
                colors = self.getColorForFaceVertexResult(fun_getcolor, minvalue, maxvalue)
                for ivh in range(len(vhandle)):
                    mesh.set_color(vhandle[ivh], color[handleColorIndex[ivh]])

    def updateMesh(self, mesh: om.PolyMesh, face2entity: Dict[int, int], mc: MeshControl,
                   const_color=[0.4, 1.0, 1.0, 1.0], fun_getcolor=None, minvalue=0, maxvalue=1):
        color = const_color
        vhandle = []
        handleColorIndex = []
        n1 = self.nodes[0]
        n2 = self.nodes[1]
        if (isinstance(self.beam_prop,StifnessBeamProperty) and mc.show_beams_with_rigid_arm_as_quad
                and isinstance(n1,NodeRigidArm) and isinstance(n2,NodeRigidArm)):
            vhandle.append(mesh.add_vertex(n1.p))
            handleColorIndex.append(0)
            vhandle.append(mesh.add_vertex(n2.p))
            handleColorIndex.append(1)
            vhandle.append(mesh.add_vertex(n2.master.p))
            handleColorIndex.append(1)
            vhandle.append(mesh.add_vertex(n1.master.p))
            handleColorIndex.append(0)
            fh = mesh.add_face(vhandle[0], vhandle[1], vhandle[2], vhandle[3])
            self.face_handles.append(fh.idx())
            face2entity[fh.idx()] = self
            if fun_getcolor is None:
                mesh.set_color(fh, [1.0, 0.6, 0.0, 1.0])
        else:
            e_x = self.nodes[0].p - self.nodes[1].p
            e_z = self.z_vec
            e_y = - np.cross(e_x, e_z)
            e_y = e_y/np.linalg.norm(e_y)
            i = 0
            pp = self.beam_prop.get_yz_point_pairs_for_plane_visualization()
            n = len(pp)/2
            yo = self.neutral_axis_offset_y
            zo = self.neutral_axis_offset_z
            while i < n:
                for j in range(2):
                    pj_y = pp[i * 2 + j][0]
                    pj_z = pp[i * 2 + j][1]
                    for k in range(2):
                        v_kj = self.nodes[k].p+(pj_y+yo)*e_y+(pj_z+zo)*e_z
                        vhandle.append(mesh.add_vertex(v_kj))
                        handleColorIndex.append(k)
                fh = mesh.add_face(vhandle[i*4], vhandle[i*4+1], vhandle[i*4+3], vhandle[i*4+2])  # 0-1-3-2
                self.face_handles.append(fh.idx())
                face2entity[fh.idx()] = self
                if fun_getcolor is None:
                    mesh.set_color(fh,[1.0, 0.8, 0.4, 1.0])
                i+=1

        # define color
        if fun_getcolor is not None:
            if mc.viewtype == ViewType.face_colors:
                color = self.getColorForFaceResult(fun_getcolor, minvalue, maxvalue)
                for fhidx in self.face_handles:
                    fh = mesh.face_handle(fhidx)
                    mesh.set_color(fh, color)
            if mc.viewtype == ViewType.face_vertex_colors:
                colors = self.getColorForFaceVertexResult(fun_getcolor, minvalue, maxvalue)
                for ivh in range(len(vhandle)):
                    mesh.set_color(vhandle[ivh], color[handleColorIndex[ivh]])

    def updateDeformedMesh(self, mesh: om.PolyMesh, face2entity: Dict[int, int], mc: MeshControl, deformation,
                           def_scale, const_color=[0.4, 1.0, 1.0, 1.0], fun_getcolor=None, minvalue=0, maxvalue=1):
        color = const_color
        vhandle = []
        handleColorIndex = []
        n1 = self.nodes[0]
        n2 = self.nodes[1]

        node1_pos = np.array([n1.p[0], n1.p[1], n1.p[2]])
        node2_pos = np.array([n2.p[0], n2.p[1], n2.p[2]])

        node1_def = np.array([deformation[n1.id][0], deformation[n1.id][1], deformation[n1.id][2]])
        node2_def = np.array([deformation[n2.id][0], deformation[n2.id][1], deformation[n2.id][2]])

        node1_disp = node1_def * def_scale + node1_pos
        node2_disp = node2_def * def_scale + node2_pos

        def_nodes = [node1_disp, node2_disp]

        if (isinstance(self.beam_prop,StifnessBeamProperty) and mc.show_beams_with_rigid_arm_as_quad
                and isinstance(n1,NodeRigidArm) and isinstance(n2,NodeRigidArm)):
            vhandle.append(mesh.add_vertex(node1_disp.p))
            handleColorIndex.append(0)
            vhandle.append(mesh.add_vertex(node2_disp.p))
            handleColorIndex.append(1)
            vhandle.append(mesh.add_vertex(node2_disp.master.p))
            handleColorIndex.append(1)
            vhandle.append(mesh.add_vertex(node1_disp.master.p))
            handleColorIndex.append(0)
            fh = mesh.add_face(vhandle[0], vhandle[1], vhandle[2], vhandle[3])
            self.face_handles.append(fh.idx())
            face2entity[fh.idx()] = self
            if fun_getcolor is None:
                mesh.set_color(fh, [1.0, 0.6, 0.0, 1.0])
        else:
            e_x = node1_disp - node2_disp
            e_z = self.z_vec
            e_y = - np.cross(e_x, e_z)
            e_y = e_y/np.linalg.norm(e_y)
            i = 0
            pp = self.beam_prop.get_yz_point_pairs_for_plane_visualization()
            n = len(pp) / 2
            yo = self.neutral_axis_offset_y
            zo = self.neutral_axis_offset_z
            while i < n:
                for j in range(2):
                    pj_y = pp[i * 2 + j][0]
                    pj_z = pp[i * 2 + j][1]
                    for k in range(2):
                        v_kj = def_nodes[k]+(pj_y+yo)*e_y+(pj_z+zo)*e_z
                        vhandle.append(mesh.add_vertex(v_kj))
                        handleColorIndex.append(k)
                fh = mesh.add_face(vhandle[i*4], vhandle[i*4+1], vhandle[i*4+3], vhandle[i*4+2])  # 0-1-3-2
                self.face_handles.append(fh.idx())
                face2entity[fh.idx()] = self
                if fun_getcolor is None:
                    mesh.set_color(fh, [1.0, 0.8, 0.4, 1.0])
                i += 1

        # define color
        if fun_getcolor is not None:
            if mc.viewtype == ViewType.face_colors:
                color = self.getColorForFaceResult(fun_getcolor, minvalue, maxvalue)
                for fhidx in self.face_handles:
                    fh = mesh.face_handle(fhidx)
                    mesh.set_color(fh, color)
            if mc.viewtype == ViewType.face_vertex_colors:
                colors = self.getColorForFaceVertexResult(fun_getcolor, minvalue, maxvalue)
                for ivh in range(len(vhandle)):
                    mesh.set_color(vhandle[ivh], color[handleColorIndex[ivh]])

class BeamElementShipStructure(BeamElement):
    def __init__(self):
        super().__init__()

    @property
    def neutral_axis_offset_z(self):
        return self.beam_prop.z_na


class TriaElement(Element):
    def __init__(self):
        super().__init__()
        pass

    def get_type(self):
        return FEMElementType.Tria

    def onTPLValue(self):
        self.face_value = self.property.tp
        return self.face_value
    def updateMesh(self, mesh: om.PolyMesh, face2entity: Dict[int, int], mc: MeshControl,
                   const_color=[0.4, 1.0, 1.0, 1.0], fun_getcolor=None, minvalue=0, maxvalue=1):

        color = const_color
        vhandle = []
        data = np.array([self.nodes[0].p[0], self.nodes[0].p[1], self.nodes[0].p[2]])
        vhandle.append(mesh.add_vertex(data))
        data = np.array([self.nodes[1].p[0], self.nodes[1].p[1], self.nodes[1].p[2]])
        vhandle.append(mesh.add_vertex(data))
        data = np.array([self.nodes[2].p[0], self.nodes[2].p[1], self.nodes[2].p[2]])
        vhandle.append(mesh.add_vertex(data))

        fh = mesh.add_face(vhandle[0], vhandle[1], vhandle[2])
        self.face_handles.append(fh.idx())
        face2entity[fh.idx()] = self
        if fun_getcolor is None:
            mesh.set_color(fh, [0.0, 0.9, 0.0, 1.0])

        # define color
        if fun_getcolor != None:
            if mc.viewtype == ViewType.face_colors:
                color =self.getColorForFaceResult(fun_getcolor, minvalue, maxvalue)
                for fhidx in self.face_handles:
                    fh = mesh.face_handle(fhidx)
                    mesh.set_color(fh, color)
            elif mc.viewtype == ViewType.face_vertex_colors:
                colors =self.getColorForFaceVertexResult(fun_getcolor, minvalue, maxvalue)
                for ivh in range(len(vhandle)):
                    mesh.set_color(vhandle[ivh], color[ivh])

    def updateDeformedMesh(self, mesh: om.PolyMesh, face2entity: Dict[int, int], mc: MeshControl, deformation,
                           def_scale, const_color=[0.4, 1.0, 1.0, 1.0], fun_getcolor=None, minvalue=0, maxvalue=1):
        node1 = self.nodes[0]
        node2 = self.nodes[1]
        node3 = self.nodes[2]

        node1_pos = np.array([node1.p[0], node1.p[1], node1.p[2]])
        node2_pos = np.array([node2.p[0], node2.p[1], node2.p[2]])
        node3_pos = np.array([node3.p[0], node3.p[1], node3.p[2]])

        node1_def = np.array([deformation[node1.id][0], deformation[node1.id][1], deformation[node1.id][2]])
        node2_def = np.array([deformation[node2.id][0], deformation[node2.id][1], deformation[node2.id][2]])
        node3_def = np.array([deformation[node3.id][0], deformation[node3.id][1], deformation[node3.id][2]])

        color = const_color
        vhandle = []
        data = node1_def * def_scale + node1_pos
        vhandle.append(mesh.add_vertex(data))
        data = node2_def * def_scale + node2_pos
        vhandle.append(mesh.add_vertex(data))
        data = node3_def * def_scale + node3_pos
        vhandle.append(mesh.add_vertex(data))


        fh = mesh.add_face(vhandle[0], vhandle[1], vhandle[2])
        self.face_handles.append(fh.idx())
        face2entity[fh.idx()] = self
        if fun_getcolor is None:
            mesh.set_color(fh, [0.0, 0.9, 0.0, 1.0])

        # define color
        if fun_getcolor != None:
            if mc.viewtype == ViewType.face_colors:
                color =self.getColorForFaceResult(fun_getcolor, minvalue, maxvalue)
                for fhidx in self.face_handles:
                    fh = mesh.face_handle(fhidx)
                    mesh.set_color(fh, color)
            elif mc.viewtype == ViewType.face_vertex_colors:
                colors =self.getColorForFaceVertexResult(fun_getcolor, minvalue, maxvalue)
                for ivh in range(len(vhandle)):
                    mesh.set_color(vhandle[ivh], color[ivh])

class StiffTriaElement(TriaElement):
    def __init__(self):
        super().__init__()
        self.layout=0
        pass
    def init(self,id):
        self.id=id


class QuadElement(Element):
    def __init__(self):
        super().__init__()
        pass

    def init(self,id):
        self.id=id

    def get_type(self):
        return FEMElementType.Quad

    def onTPLValue(self):
        self.face_value = self.property.tp
        return self.face_value
    def updateMesh(self, mesh: om.PolyMesh, face2entity: Dict[int, int], mc: MeshControl,
                   const_color=[0.4, 1.0, 1.0, 1.0], fun_getcolor=None, minvalue=0, maxvalue=1):

        color = const_color
        vhandle = []
        data = np.array([self.nodes[0].p[0], self.nodes[0].p[1], self.nodes[0].p[2]])
        vhandle.append(mesh.add_vertex(data))
        data = np.array([self.nodes[1].p[0], self.nodes[1].p[1], self.nodes[1].p[2]])
        vhandle.append(mesh.add_vertex(data))
        data = np.array([self.nodes[2].p[0], self.nodes[2].p[1], self.nodes[2].p[2]])
        vhandle.append(mesh.add_vertex(data))
        data = np.array([self.nodes[3].p[0], self.nodes[3].p[1], self.nodes[3].p[2]])
        vhandle.append(mesh.add_vertex(data))

        fh = mesh.add_face(vhandle[0], vhandle[1], vhandle[2], vhandle[3])
        self.face_handles.append(fh.idx())
        face2entity[fh.idx()] = self
        if fun_getcolor is None:
            mesh.set_color(fh, [0.0, 0.0, 0.9, 1.0])

        # define color
        if fun_getcolor != None:
            if mc.viewtype == ViewType.face_colors:
                color =self.getColorForFaceResult(fun_getcolor, minvalue, maxvalue)
                for fhidx in self.face_handles:
                    fh = mesh.face_handle(fhidx)
                    mesh.set_color(fh, color)
            elif mc.viewtype == ViewType.face_vertex_colors:
                colors =self.getColorForFaceVertexResult(fun_getcolor, minvalue, maxvalue)
                for ivh in range(len(vhandle)):
                    mesh.set_color(vhandle[ivh], color[ivh])

    def updateDeformedMesh(self, mesh: om.PolyMesh, face2entity: Dict[int, int], mc: MeshControl, deformation,
                           def_scale, const_color=[0.4, 1.0, 1.0, 1.0], fun_getcolor=None, minvalue=0, maxvalue=1):
        node1 = self.nodes[0]
        node2 = self.nodes[1]
        node3 = self.nodes[2]
        node4 = self.nodes[3]

        node1_pos = np.array([node1.p[0], node1.p[1], node1.p[2]])
        node2_pos = np.array([node2.p[0], node2.p[1], node2.p[2]])
        node3_pos = np.array([node3.p[0], node3.p[1], node3.p[2]])
        node4_pos = np.array([node4.p[0], node4.p[1], node4.p[2]])

        node1_def = np.array([deformation[node1.id][0], deformation[node1.id][1], deformation[node1.id][2]])
        node2_def = np.array([deformation[node2.id][0], deformation[node2.id][1], deformation[node2.id][2]])
        node3_def = np.array([deformation[node3.id][0], deformation[node3.id][1], deformation[node3.id][2]])
        node4_def = np.array([deformation[node4.id][0], deformation[node4.id][1], deformation[node4.id][2]])

        color = const_color
        vhandle = []
        data = node1_def * def_scale + node1_pos
        vhandle.append(mesh.add_vertex(data))
        data = node2_def * def_scale + node2_pos
        vhandle.append(mesh.add_vertex(data))
        data = node3_def * def_scale + node3_pos
        vhandle.append(mesh.add_vertex(data))
        data = node4_def * def_scale + node4_pos
        vhandle.append(mesh.add_vertex(data))
        
        fh = mesh.add_face(vhandle[0], vhandle[1], vhandle[2], vhandle[3])
        self.face_handles.append(fh.idx())
        face2entity[fh.idx()] = self
        if fun_getcolor is None:
            mesh.set_color(fh, [0.0, 0.0, 0.9, 1.0])

        # define color
        if fun_getcolor != None:
            if mc.viewtype == ViewType.face_colors:
                color =self.getColorForFaceResult(fun_getcolor, minvalue, maxvalue)
                for fhidx in self.face_handles:
                    fh = mesh.face_handle(fhidx)
                    mesh.set_color(fh, color)
            elif mc.viewtype == ViewType.face_vertex_colors:
                colors =self.getColorForFaceVertexResult(fun_getcolor, minvalue, maxvalue)
                for ivh in range(len(vhandle)):
                    mesh.set_color(vhandle[ivh], color[ivh])

class StiffQuadElement(QuadElement):
    def __init__(self):
        super().__init__()
        self.layout=0
        pass

    def get_type(self):
        return FEMElementType.StiffQuad

    def init(self,id):
        self.id=id
    def updateTriMesh(self,mesh:om.TriMesh,mc:MeshControl, const_color = [0.4, 1.0, 1.0, 1.0],
                   fun_getcolor=None,minvalue = 0,maxvalue = 1 ):

        color = const_color
        vhandle = []
        fhs = []
        data = np.array([self.nodes[0].p[0], self.nodes[0].p[1], self.nodes[0].p[2]])
        vhandle.append(mesh.add_vertex(data))
        data = np.array([self.nodes[1].p[0], self.nodes[1].p[1], self.nodes[1].p[2]])
        vhandle.append(mesh.add_vertex(data))
        data = np.array([self.nodes[2].p[0], self.nodes[2].p[1], self.nodes[2].p[2]])
        vhandle.append(mesh.add_vertex(data))
        data = np.array([self.nodes[3].p[0], self.nodes[3].p[1], self.nodes[3].p[2]])
        vhandle.append(mesh.add_vertex(data))

        fhs.append(mesh.add_face(vhandle[0], vhandle[1], vhandle[2]))
        fhs.append(mesh.add_face(vhandle[0], vhandle[2], vhandle[3]))

        # define color
        if fun_getcolor != None:
            if mc.viewtype == ViewType.face_colors:
                color =self.getColorForFaceResult(fun_getcolor, minvalue, maxvalue)
                for fh in fhs:
                    mesh.set_color(fh, color)
            elif mc.viewtype == ViewType.face_vertex_colors:
                colors =self.getColorForFaceVertexResult(fun_getcolor, minvalue, maxvalue)
                for ivh in range(len(vhandle)):
                    mesh.set_color(vhandle[ivh], color[ivh])

    def updateMesh(self, mesh: om.PolyMesh, face2entity: Dict[int, int], mc: MeshControl,
                   const_color=[0.4, 1.0, 1.0, 1.0], fun_getcolor=None, minvalue=0, maxvalue=1):

        color = const_color
        vhandle = []
        data = np.array([self.nodes[0].p[0], self.nodes[0].p[1], self.nodes[0].p[2]])
        vhandle.append(mesh.add_vertex(data))
        data = np.array([self.nodes[1].p[0], self.nodes[1].p[1], self.nodes[1].p[2]])
        vhandle.append(mesh.add_vertex(data))
        data = np.array([self.nodes[2].p[0], self.nodes[2].p[1], self.nodes[2].p[2]])
        vhandle.append(mesh.add_vertex(data))
        data = np.array([self.nodes[3].p[0], self.nodes[3].p[1], self.nodes[3].p[2]])
        vhandle.append(mesh.add_vertex(data))

        fh= mesh.add_face(vhandle[0], vhandle[1], vhandle[2], vhandle[3])
        self.face_handles.append(fh.idx())
        face2entity[fh.idx()] = self
        if fun_getcolor is None:
            mesh.set_color(fh, [0.0, 0.0, 0.9, 1.0])

        # define color
        if fun_getcolor != None:
            if mc.viewtype == ViewType.face_colors:
                color =self.getColorForFaceResult(fun_getcolor, minvalue, maxvalue)
                for fhidx in self.face_handles:
                    fh = mesh.face_handle(fhidx)
                    mesh.set_color(fh, color)
            elif mc.viewtype == ViewType.face_vertex_colors:
                colors =self.getColorForFaceVertexResult(fun_getcolor, minvalue, maxvalue)
                for ivh in range(len(vhandle)):
                    mesh.set_color(vhandle[ivh], color[ivh])


class BC_and_Load():
    def __init__(self,id = -1):
        self._id= id
        pass

    @property
    def id(self):
        return self._id

    @id.setter
    def id(self, value):
        self._id = value

    def updateMesh(self):
        pass


class DoffBasedBCandLoad(BC_and_Load):
    def __init__(self, id, values:List[float], dofs:List[int]=None):
        super().__init__(id)
        if dofs is None:
            self._dofs = np.array(list(range(1, len(values) + 1)))
        else:
            self._dofs = np.array(dofs) # indice values are 1 based (from 1 to 6)
        self._values = np.array(values)

    @property
    def values(self):
        return self._values

    @property
    def dofs(self):
        return self._dofs

    def get_num_data(self):
        return self._dofs.size

    def get_value(self,idata):
        return self._values[idata]

    def get_6dof_values(self):
        if self.get_num_data() == 6:
            return self.values.copy()
        else:
            vals = [0.0]*6
            for i in range(self._dofs):
                vals[self._dofs[i]]=self.values[i]
            return vals

    def get_position(self,idata):
        return self._dofs[idata]

    def have_force(self):
        for i in self._dofs:
            if i < 4:
                return True
        return False
    def have_moment(self):
        for i in self._dofs:
            if i > 3:
                return True
        return False

    def have_translation(self):
        return self.have_force()

    def have_rotation(self):
        return self.have_moment()

    def get_force_vector(self):
        vec = np.zeros(3)
        for i in range(len(self._dofs)):
            ix=self._dofs[i]
            if ix < 4:
                vec[ix-1]=self._values[i]
        return vec

    def get_moment_vector(self):
        vec = np.zeros(3)
        for i in range(len(self._dofs)):
            ix=self._dofs[i]
            if ix > 3:
                vec[ix-4]=self._values[i]
        return vec

    def updateMesh(self):
        fmsum = self._getComponetsSum()
        if fmsum[0] > 0: # force exist
            pass
        elif fmsum[1] > 0: # moment exist
            pass

        pass

class PointDoffBasedLoad(DoffBasedBCandLoad):
    def __init__(self, id, node, values:List[float], dofs:List[int]):
        super().__init__(id, values, dofs)
        self._node:Node = node
        pass
    def updateMesh(self):
        pass

class PointNodalBC(PointDoffBasedLoad):
    def __init__(self, id, node, values:List[float], dofs:List[int]):
        super().__init__(id, node, values, dofs)
        pass
    def updateMesh(self):
        pass

class GroupDoffBasedLoad(DoffBasedBCandLoad):
    def __init__(self, id, group, values:List[float], dofs:List[int]):
        super().__init__(id, values, dofs)
        self._group:Group = group
        pass
    def updateMesh(self):
        pass
    @property
    def id_group(self):
        return self._group.id

    @property
    def group(self):
        return self._group

class NodalLineLoad(GroupDoffBasedLoad):
    def __init__(self, id, group, values1:List[float], values2:List[float], dofs:List[int]):
        super().__init__(id, group, values1, dofs)
        self._values2 = np.array(values2)
        pass

    @property
    def values1(self):
        return self._values

    @property
    def values2(self):
        return self._values2

    def get_eqivalent_point_loads(self,id_nextload,id_nextgroup):
        # if the line load is constant for each dof, only two nodal groups are necessary
        # otherwise one group per node is necessary groups are necessary
        # l = length, q1 load on first, q2 load on second
        # F1 = l*(2*q1+q2)/6, F2 = l*(2*q2+q1)/6
        llinenodevalues = [self.values1]
        num_internal_nodes = self.group.num_items-2
        for i in range(1,self.group.num_items-1):
            nodevals = self.values1+(self.values2-self.values1)*i/(num_internal_nodes+1)
            llinenodevalues.append(nodevals)
        llinenodevalues.append(self.values2)
        dict_node_values={}
        for node in self.group.items:
            dict_node_values[node] = np.zeros((3))
        for i in range(self.group.num_items-1):
            node1 = self.group.items[i]
            node2 = self.group.items[i+1]
            l =Node.get_distance_betweeen_nodes(node1,node2)
            q1 = llinenodevalues[i]
            q2 = llinenodevalues[i+1]
            F1 = l * (2.0 * q1 + q2) / 6.0
            F2 = l * (2.0 * q2 + q1) / 6.0
            dict_node_values[node1] = dict_node_values[node1] + F1
            dict_node_values[node2] = dict_node_values[node2] + F2
        # determine minimal number of groups
        dict_previously_assigned={}
        nnodes = self.group.num_items
        for i1 in range(nnodes):
            cur_node = self.group.items[i1]
            if cur_node in dict_previously_assigned:
                continue
            dict_previously_assigned[cur_node] = cur_node
            cur_value = dict_node_values[cur_node]
            for i2 in range(i1+1,nnodes):
                test_node = self.group.items[i2]
                test_value = dict_node_values[test_node]
                isclose = np.allclose(cur_value,test_value, rtol=1e-05, atol=1e-08)
                if isclose:
                    dict_previously_assigned[test_node]=cur_node
        dict_nodal_load_group_items = {}
        for slave,master in dict_previously_assigned.items():
            if master not in dict_nodal_load_group_items:
                grlist = []
                dict_nodal_load_group_items[master]=grlist
            else:
                grlist = dict_nodal_load_group_items[master]
            grlist.append(slave)
        # create groups and loads
        i=0
        nn = len(dict_nodal_load_group_items)
        nodal_groups = [None]*nn
        nodal_loads = [None]*nn
        for key,grlist in dict_nodal_load_group_items.items():
            new_group = NodeGroup(id_nextgroup,'gr_'+str(i+1)+'_from_NLL_id_'+str(self.id))
            new_group.extend_items(grlist)
            nodal_groups[i]=new_group
            id_nextgroup+=1
            new_load = GroupDoffBasedLoad(id_nextload,new_group,dict_node_values[key],self.dofs)
            id_nextload+=1
            nodal_loads[i] = new_load
            i+=1
        return nodal_loads,nodal_groups


    def updateMesh(self):
        pass

class GroupNodalBC(GroupDoffBasedLoad):
    def __init__(self, id, group, values:List[float], dofs:List[int]):
        super().__init__(id, group, values, dofs)
    def updateMesh(self):
        pass

class PressureLoad(BC_and_Load):
    def __init__(self,id,pressure:float,flip:bool=False):
        super().__init__(id)
        self.pressure = pressure
        if pressure > 0.0 and flip:
            self._flip = flip

    @property
    def pressure(self):
        return self._press

    @property
    def signed_pressure(self):
        if self.flip:
            return -self._press
        else:
            return self._press

    @pressure.setter
    def pressure(self, value):
        if value < 0.0:
            self._flip = True
            self._press= - value
        else:
            self._press = value
            self._flip = False
        pass

    @property
    def flip(self):
        return self._flip

    @flip.setter
    def flip(self, value):
        self._flip = value
        pass

class GroupPressureLoad(PressureLoad):
    def __init__(self,id,group,pressure:float,flip:bool=False):
        super().__init__(id,pressure)
        if pressure > 0.0 and flip:
            self._flip = flip
        self._group:Group = group
        pass
    def updateMesh(self):
        pass
    @property
    def id_group(self):
        return self._group.id

class AccelerationLoad(DoffBasedBCandLoad):
    def __init__(self, id, values:List[float], dofs:List[int]=None):
        super().__init__(id, values, dofs)


    def updateMesh(self):
        pass

class LoadCaseData():
    def __init__(self):
        pass

    def updateMesh(self):
        pass

class LoadCase():
    def __init__(self,id,name):
        self._loads: List[BC_and_Load] = []
        self._boundaryconditions: List[BC_and_Load] = []
        self._id = id
        self._name = name
        pass

    def add_load(self, load:BC_and_Load):
        self._loads.append(load)

    def add_boundarycondition(self, bc:BC_and_Load):
        self._boundaryconditions.append(bc)

    @property
    def loads(self):
        return self._loads

    @property
    def boundaryconditions(self):
        return self._boundaryconditions

    @property
    def num_loads(self):
        return len(self._loads)

    @property
    def num_bcs(self):
        return len(self._boundaryconditions)

    @property
    def name(self):
        return self._name


    @property
    def id(self):
        return self._id

    def updateMesh(self):
        pass


class Units():
    def __init__(self):
        self.user2si_length=1
        self.user2si_force=1
        self.name_length='m'
        self.name_force = 'N'
        pass
class MaestroElementAssociation():
    def __init__(self):
        self.strakeGirder2fe = {}
        self.strakePlate2fe = {}
        self.strakeFrame2fe = {}
        self.endPoint2node = {}
        self.endPointStrakes = {}

    def addStrakePlate(self,key:int,elList:[]):
        self.strakePlate2fe[key]=elList

    def addStrakeGirder(self,key:int,elList:[]):
        self.strakeGirder2fe[key]=elList

    def addStrakeFrame(self,key:int,elList:[]):
        self.strakeFrame2fe[key]=elList

    def addEndPointNode(self,key:int,nodeFeTag:int):
        feTagList = self.endPoint2node.setdefault(key,[])
        feTagList.append(nodeFeTag)

    def addEndPointStrake(self, key:int, strakeID:int):
        strakeList = self.endPointStrakes.setdefault(key,[])
        strakeList.append(strakeID)

    def getPlateElsForEnpoint(self, endpointID):
        connectedElements=[]
        strakeList = self.endPointStrakes.get(endpointID)
        if strakeList is not None:
            for idStrake in strakeList:
                for idEl in self.strakePlate2fe[idStrake]:
                    connectedElements.append(idEl)
        return connectedElements

    def getPlateElemForStrake(self,strakeID):
        return self.strakePlate2fe.get(strakeID)

    def getGirderBeamElemForStrake(self,strakeID):
        return self.strakeGirder2fe.get(strakeID)


class LusaElementAssociation():
    def __init__(self):
        self.spc2fe = {}
        self.plate2fe = {}
        self.hc2fe = {}

    def addPlate(self, key:int, strakeID):
        self.plate2fe[key]=strakeID
    def addSPC(self,key:int,strakeID):
        self.spc2fe[key]=strakeID
    def addHC(self,key:int,strakeID):
        self.hc2fe[key]=strakeID
        
