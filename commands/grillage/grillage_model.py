"""
Tool for Grillage Structure Analysis
University of Zagreb, Faculty of Mechanical Engineering and Naval Architecture
Department of Naval Architecture and Ocean Engineering


MODULE FOR GRILLAGE STRUCTURE DEFINITION


Minimum input data for grillgeo definition with initial uniform spacing of all elements:
    Main dimensions: L, B
    Number of elements: longitudinal and transverse primary supporting members, stiffeners between primary supporting members
    Stiffener orientation: LONGITUDINAL / TRANSVERSE
    Corrosion addition value: tc
    Primary supporting member beam property: hw, tw, bf, tf, material
    Stiffener HP property: hw, tw, material
    Stiffener layout: beam property, definition type, definition value
    Plate property: t, material


Measurement units:
    * Grillage main dimensions: [m]
    * Primary supporting members relative distance values: [m]
    * All beam property, plate thickness, corrosion addition input dimensions: [mm]
    * Material density: [kg/m3]
    * Modulus of elasticity, yield strength: [N/mm2]


Global coordinate system (csy):
    x - longitudinal axis, in the direction of overall length L
    y - transverse axis, in the direction of overall width B
    z - vertical axis, oriented from the flange to the plating

"""

import numpy as np
import itertools
from enum import Enum


class BeamDirection(Enum):
    TRANSVERSE = 0
    LONGITUDINAL = 1


class FlangeDirection(Enum):   # L beam flange orientation
    INWARD = 1
    OUTWARD = 2


class Ref(Enum):    # Reference edge
    EDGE1 = 1
    EDGE2 = 2


class AOS(Enum):    # Axis Of Symmetry
    TRANSVERSE = 0
    LONGITUDINAL = 1
    BOTH = 2
    NONE = 3


class Node:
    def __init__(self, id_: int, x=0.0, y=0.0, z=0.0):
        self._id = id_
        self._coords = np.array([x, y, z])

    @property
    def id(self):
        return self._id

    @property
    def coords(self):
        return self._coords


class MaterialProperty:
    def __init__(self, id_, E, v, ro, Reh, name):
        self._id_ = id_
        self._E = float(E)         # Young's modulus
        self._v = float(v)         # Poisson's ratio
        self._ro = float(ro)       # Density
        self._Reh = float(Reh)     # Yield strength
        self._name = name          # Material name

    @property
    def id(self):
        return self._id_

    @property
    def E(self):
        return self._E

    @E.setter
    def E(self, value):
        self._E = value

    @property
    def v(self):
        return self._v

    @v.setter
    def v(self, value):
        self._v = value

    @property
    def ro(self):
        return self._ro

    @ro.setter
    def ro(self, value):
        self._ro = value

    @property
    def Reh(self):
        return self._Reh

    @Reh.setter
    def Reh(self, value):
        self._Reh = value

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        self._name = value


class CorrosionAddition:
    def __init__(self, id_, tc):
        self._id_ = id_
        self._tc = float(tc)

    @property
    def id(self):
        return self._id_

    @property
    def tc(self):
        return self._tc

    @tc.setter
    def tc(self, value):
        self._tc = value


class PlateProperty:
    def __init__(self, id_, tp, plate_mat):
        self._id_ = id_
        self._tp = float(tp)
        self._plate_mat = plate_mat

    @property
    def id(self):
        return self._id_

    @property
    def tp(self):
        return self._tp

    @tp.setter
    def tp(self, value):
        self._tp = value

    @property
    def plate_mat(self):
        return self._plate_mat

    @plate_mat.setter
    def plate_mat(self, value):
        self._plate_mat = value

    @staticmethod
    def tp_net(corr_add: CorrosionAddition, tp):
        # Net attached plating thickness
        if tp > 0 and tp > corr_add.tc:
            tp_net = tp - corr_add.tc
        else:
            tp_net = 0
        return tp_net


class BeamProperty:
    def __init__(self, id_):
        self._id_ = id_

    @property
    def id(self):
        return self._id_

    def get_z_na_I(self, bp, tp, corr_add: CorrosionAddition):
        return 0.0

    def get_z_max_I(self, bp, tp, corr_add: CorrosionAddition):
        return 0.0

    def get_Iy_I(self, bp, tp, corr_add: CorrosionAddition):
        return 0.0

    def get_Iw(self, corr_add: CorrosionAddition):
        return 0.0

    def get_Ip(self, corr_add: CorrosionAddition):
        return 0.0

    def get_It(self, corr_add: CorrosionAddition):
        return 0.0

    def getArea(self):
        return 0.0

    def getWmin(self, bp, tp, corr_add: CorrosionAddition):  # Minimum net section modulus
        zna = self.get_z_na_I(bp, tp, corr_add)
        zmax = self.get_z_max_I(bp, tp, corr_add)
        iy = self.get_Iy_I(bp, tp, corr_add)
        wp = iy / (zna * 0.1)                       # Net section modulus at the plating - Wp, [cm3]
        wf = iy / ((zmax - zna) * 0.1)              # Net section modulus at the flange - Wf, [cm3]
        if wp < wf:
            wmin = wp
            return wmin
        elif wf <= wp:
            wmin = wf
            return wmin

    @property
    def beam_type(self):
        if isinstance(self, TBeamProperty) and not \
                (isinstance(self, FBBeamProperty) or isinstance(self, LBeamProperty)):
            return "T"

        elif isinstance(self, LBeamProperty):
            return "L"

        elif isinstance(self, FBBeamProperty):
            return "FB"


class TBeamProperty(BeamProperty):
    def __init__(self, id_, hw, tw, bf, tf, mat: MaterialProperty):
        super().__init__(id_)
        self._hw = float(hw)  # Web height
        self._tw = float(tw)  # Web thickness
        self._bf = float(bf)  # Flange width
        self._tf = float(tf)  # Flange thickness
        self._mat = mat       # Material

    @property
    def hw(self):
        return self._hw

    @hw.setter
    def hw(self, value):
        self._hw = value

    @property
    def tw(self):
        return self._tw

    @tw.setter
    def tw(self, value):
        self._tw = value

    @property
    def bf(self):
        return self._bf

    @bf.setter
    def bf(self, value):
        self._bf = value

    @property
    def tf(self):
        return self._tf

    @tf.setter
    def tf(self, value):
        self._tf = value

    @property
    def mat(self):
        return self._mat

    @mat.setter
    def mat(self, value):
        self._mat = value

    def hw_net(self, corr_add: CorrosionAddition, tp):
        # Net web height of a T profile with attached plating
        if self.tf > 0 and tp > 0:              # If both flange and attached plating exist, net web height increases by tc
            hw_net = self.hw + corr_add.tc
        elif self.tf == 0 and tp == 0:          # If there is no flange or attached plating, beam becomes a FB profile
            hw_net = self.hw - corr_add.tc
        else:                                   # Takes into account combinations of only flange and only attached plating
            hw_net = self.hw                    # web height increases by tc/2 and decreases by tc/2, so there is no difference
        return hw_net

    def tw_net(self, corr_add: CorrosionAddition):
        tw_net = self.tw - corr_add.tc
        return tw_net

    def bf_net(self, corr_add: CorrosionAddition):
        # Net flange width of a T profile
        if self.bf > 0 and self.bf > corr_add.tc:
            bf_net = self.bf - corr_add.tc
        else:
            bf_net = 0
        return bf_net

    def tf_net(self, corr_add: CorrosionAddition):
        # Net flange thickness of a T profile
        if self.tf > 0 and self.tf > corr_add.tc:
            tf_net = self.tf - corr_add.tc
        else:
            tf_net = 0
        return tf_net

    def getShArea_T(self, tp_gross, corr_add: CorrosionAddition):  # Net shear area of T profile - Ash, [cm2]
        hw = self.hw_net(corr_add, tp_gross)
        tw = self.tw_net(corr_add)
        Ash = (hw * tw) * 0.01
        return Ash

    @property
    def getArea(self):  # Gross cross sectional area of T profile, without attached plating - A, [cm2]
        hw = self.hw
        tw = self.tw
        bf = self.bf
        tf = self.tf
        A = (bf * tf + hw * tw) * 0.01
        return A

    def getArea_I(self, bp, tp_gross, corr_add: CorrosionAddition):  # Net area of I beam with attached plating - A, [cm2]
        hw = self.hw_net(corr_add, tp_gross)
        tw = self.tw_net(corr_add)
        bf = self.bf_net(corr_add)
        tf = self.tf_net(corr_add)
        tp_net = PlateProperty.tp_net(corr_add, tp_gross)
        A = (bf * tf + hw * tw + bp * tp_net) * 0.01
        return A

    def get_z_na_I(self, bp, tp_gross, corr_add: CorrosionAddition):  # Neutral axis measured from the top of the attached plating - zna, [mm]
        hw = self.hw_net(corr_add, tp_gross)
        tw = self.tw_net(corr_add)
        bf = self.bf_net(corr_add)
        tf = self.tf_net(corr_add)
        tp = PlateProperty.tp_net(corr_add, tp_gross)
        A = self.getArea_I(bp, tp_gross, corr_add)
        zna = ((bf * tf * (tp + hw + tf * 0.5)) + (hw * tw * (tp + hw * 0.5) + (bp * (tp ** 2) * 0.5))) / (100 * A)
        return zna

    def get_z_max_I(self, bp, tp, corr_add: CorrosionAddition):  # Maximum z coordinate, from the bottom of the flange to the top of the plating
        zmax = self._tf + self._hw + tp - corr_add.tc
        return zmax

    def get_Iy_I(self, bp, tp_gross, corr_add: CorrosionAddition):  # Net moment of inertia around the local y axis - Iy, [cm4]
        hw = self.hw_net(corr_add, tp_gross)
        tw = self.tw_net(corr_add)
        bf = self.bf_net(corr_add)
        tf = self.tf_net(corr_add)
        tp = PlateProperty.tp_net(corr_add, tp_gross)
        zna = self.get_z_na_I(bp, tp_gross, corr_add)
        Iy = (bf * tf * ((tf ** 2) / 12 + (tf * 0.5 + hw + tp - zna) ** 2) + hw * tw * (
                (hw ** 2) / 12 + (hw * 0.5 + tp - zna) ** 2) + bp * tp * (
                (tp ** 2) / 12 + (zna - tp * 0.5) ** 2))
        return Iy * (10 ** (-4))

    def get_Iw(self, corr_add: CorrosionAddition):  # Sectorial moment of inertia in [cm6]
        hw = self._hw + corr_add.tc
        bf = self.bf_net(corr_add)
        tf = self.tf_net(corr_add)
        Iw = (tf * (bf ** 3) * (hw ** 2)) * 10 ** (-6)
        return Iw

    def get_Ip(self, corr_add: CorrosionAddition):  # Polar moment of inertia in [cm4]
        hw = self._hw + corr_add.tc
        tw = self.tw_net(corr_add)
        bf = self.bf_net(corr_add)
        tf = self.tf_net(corr_add)
        Ip = ((((hw ** 3) * tw) / 3) + ((hw ** 2) * bf * tf)) * 10 ** (-4)
        return Ip

    def get_It(self, corr_add: CorrosionAddition):  # St Venant's moment of inertia in [cm4]
        hw = self._hw + corr_add.tc
        tw = self.tw_net(corr_add)
        bf = self.bf_net(corr_add)
        tf = self.tf_net(corr_add)
        It = (1 / 3) * ((hw * tw ** 3) + (bf * tf ** 3) * (1 - 0.63 * (tf / bf))) * 10 ** (-4)
        return It


class FBBeamProperty(TBeamProperty):
    def __init__(self, id_, hw, tw, mat, bf=0.0, tf=0.0):
        super().__init__(id_, hw, tw, bf, tf, mat)
        self._hw = float(hw)  # Flat Bar height
        self._tw = float(tw)  # Flat Bar thickness
        self._bf = float(bf)  # Second plate attached plating width, default bf=0.0 as ordinary FB stiffner
        self._tf = float(tf)  # Second plate attached plating thickness, default tf=0.0 as ordinary FB stiffner
        self._mat = mat       # Material

    def get_Iw(self, corr_add: CorrosionAddition):  # Sectorial moment of inertia in [cm6]
        hw = self._hw + corr_add.tc
        tw = self.tw_net(corr_add)
        Iw = (((hw ** 3) * (tw ** 3)) / 36) * 10 ** (-6)
        return Iw

    def get_Ip(self, corr_add: CorrosionAddition):  # Polar moment of inertia in [cm4]
        hw = self._hw + corr_add.tc
        tw = self.tw_net(corr_add)
        Ip = (((hw ** 3) * tw) / 3) * 10 ** (-4)
        return Ip

    def get_It(self, corr_add: CorrosionAddition):  # St Venant's moment of inertia in [cm4]
        hw = self._hw + corr_add.tc
        tw = self.tw_net(corr_add)
        It = (((hw ** 3) * tw) / 3) * 10 ** (-4)
        return It


class LBeamProperty(TBeamProperty):
    def __init__(self, id_, hw, tw, bf, tf, mat: MaterialProperty):
        super().__init__(id_, hw, tw, bf, tf, mat)
        self._hw = float(hw)  # Web height
        self._tw = float(tw)  # Web thickness
        self._bf = float(bf)  # Flange width
        self._tf = float(tf)  # Flange thickness
        self._mat = mat       # Material

    def get_Iw(self, corr_add: CorrosionAddition):  # Sectorial moment of inertia in [cm6]
        hw = self.hw + corr_add.tc
        tw = self.tw - corr_add.tc
        bf = self.bf - corr_add.tc
        tf = self.tf - corr_add.tc
        Iw = ((bf ** 3 * hw ** 2) / (12 * (bf + hw) ** 2)) *\
             (tf * (bf ** 2 + 2 * bf * hw + 4 * hw ** 2) + 3 * tw * bf * hw) * 10 ** (-6)
        return Iw


class BulbBeamProperty(BeamProperty):
    def __init__(self, id_, hw_HP, tw_HP, mat: MaterialProperty):
        super().__init__(id_)
        self._id_ = id_
        self._hw_HP = float(hw_HP)  # HP profile height
        self._tw_HP = float(tw_HP)  # Web thickness (net)
        self._mat = mat  # Material

    @property
    def hw_HP(self):
        return self._hw_HP

    @hw_HP.setter
    def hw_HP(self, value):
        self.hw_HP = value

    @property
    def tw_HP(self):
        return self._tw_HP

    @tw_HP.setter
    def tw_HP(self, value):
        self._tw_HP = value

    @property
    def mat(self):
        return self._mat

    @mat.setter
    def mat(self, value):
        self._mat = value

    # Equivalent angle calculation according to:
    # IACS Common Structural Rules, July 2012, Chapter 3, Section 6, 4.1.1 Stiffener profile with a bulb section
    def alpha(self):
        hw_HP = self._hw_HP
        if hw_HP > 120:
            alpha = 1
            return alpha
        if hw_HP <= 120:
            alpha = 1.1 + ((120 - hw_HP) ** 2 / 3000)
            return alpha

    @property
    def hw_ekv(self):  # Web height of equivalent L profile
        hw_HP = self._hw_HP
        hw = hw_HP - (hw_HP / 9.2) + 2
        return hw

    @property
    def tw_ekv(self):  # Web thickness of equivalent L profile
        tw = self._tw_HP
        return tw

    @property
    def bf_ekv(self):  # Flange width of equivalent L profile
        hw_HP = self._hw_HP
        tw_HP = self._tw_HP
        alpha = self.alpha()
        bf = alpha * (tw_HP + (hw_HP / 6.7) - 2)
        return bf

    @property
    def tf_ekv(self):  # Flange thickness of equivalent L profile
        hw_HP = self._hw_HP
        tf = (hw_HP / 9.2) - 2
        return tf

    def getShArea_HP(self, corr_add: CorrosionAddition):  # Net shear area of HP profile - Ash, [cm2]
        hw = self.hw_ekv + corr_add.tc
        tw = self.tw_ekv - corr_add.tc
        Ash = (hw * tw) * 0.01
        return Ash

    @property
    def getArea(self):  # Gross cross sectional area of HP stiffener, without attached plating - A, [cm2]
        hw = self.hw_ekv
        tw = self.tw_ekv
        bf = self.bf_ekv
        tf = self.tf_ekv
        A = (bf * tf + hw * tw) * 0.01
        return A

    def getArea_I(self, bp, tp_gross, corr_add: CorrosionAddition):  # Net area of HP stiffener with attached plating - A, [cm2]
        hw = self.hw_ekv + corr_add.tc
        tw = self.tw_ekv - corr_add.tc
        bf = self.bf_ekv - corr_add.tc
        tf = self.tf_ekv - corr_add.tc
        tp_net = tp_gross - corr_add.tc
        A = (bf * tf + hw * tw + bp * tp_net) * 0.01
        return A

    def get_z_na_I(self, bp, tp_gross, corr_add: CorrosionAddition):  # Neutral axis measured from the top of the attached plating - zna, [mm]
        hw = self.hw_ekv + corr_add.tc
        tw = self.tw_ekv - corr_add.tc
        bf = self.bf_ekv - corr_add.tc
        tf = self.tf_ekv - corr_add.tc
        tp = tp_gross - corr_add.tc
        A = self.getArea_I(bp, tp_gross, corr_add)
        zna = ((bf * tf * (tp + hw + tf * 0.5)) + (hw * tw * (tp + hw * 0.5) + (bp * (tp ** 2) * 0.5))) / (100 * A)
        return zna

    def get_z_max_I(self, bp, tp, corr_add: CorrosionAddition):  # Maximum z coordinate, from the bottom of the flange to the top of the plating
        zmax = self.tf_ekv + self.hw_ekv + tp - corr_add.tc
        return zmax

    def get_Iy_I(self, bp, tp_gross, corr_add: CorrosionAddition):  # Net moment of inertia around the local y axis - Iy, [cm4]
        hw = self.hw_ekv + corr_add.tc
        tw = self.tw_ekv - corr_add.tc
        bf = self.bf_ekv - corr_add.tc
        tf = self.tf_ekv - corr_add.tc
        tp = tp_gross - corr_add.tc
        zna = self.get_z_na_I(bp, tp_gross, corr_add)
        Iy = (bf * tf * ((tf ** 2) / 12 + (tf * 0.5 + hw + tp - zna) ** 2) + hw * tw * (
                (hw ** 2) / 12 + (hw * 0.5 + tp - zna) ** 2) + bp * tp * (
                (tp ** 2) / 12 + (zna - tp * 0.5) ** 2))
        return Iy * (10 ** (-4))

    def get_Iw(self, corr_add: CorrosionAddition):  # Sectorial moment of inertia in [cm6]
        hw = self.hw_ekv + corr_add.tc
        tw = self.tw_ekv - corr_add.tc
        bf = self.bf_ekv - corr_add.tc
        tf = self.tf_ekv - corr_add.tc
        Iw = ((bf ** 3 * hw ** 2) / (12 * (bf + hw) ** 2)) *\
             (tf * (bf ** 2 + 2 * bf * hw + 4 * hw ** 2) + 3 * tw * bf * hw) * 10 ** (-6)
        return Iw

    def get_Ip(self, corr_add: CorrosionAddition):  # Polar moment of inertia in [cm4]
        hw = self.hw_ekv + corr_add.tc
        tw = self.tw_ekv - corr_add.tc
        bf = self.bf_ekv - corr_add.tc
        tf = self.tf_ekv - corr_add.tc
        Ip = ((((hw ** 3) * tw) / 3) + ((hw ** 2) * bf * tf)) * 10 ** (-4)
        return Ip

    def get_It(self, corr_add: CorrosionAddition):  # St Venant's moment of inertia in [cm4]
        hw = self.hw_ekv + corr_add.tc
        tw = self.tw_ekv - corr_add.tc
        bf = self.bf_ekv - corr_add.tc
        tf = self.tf_ekv - corr_add.tc
        It = (1 / 3) * ((hw * tw ** 3) + (bf * tf ** 3) * (1 - 0.63 * (tf / bf))) * 10 ** (-4)
        return It


class HatBeamProperty(BeamProperty):
    def __init__(self, id_, h, t, bf, fi, mat):
        super().__init__(id_)
        self._id_ = id_
        self._h = float(h)      # Profile height - distance between inner flange surface and closest plating surface
        self._t = float(t)      # Web and flange thickness
        self._bf = float(bf)    # Flange width between outer corners of the Hat profile
        self._fi = float(fi)    # Web angle between plating and web
        self._mat = mat         # Material

    @property
    def h(self):
        return self._h

    @h.setter
    def h(self, value):
        self._h = value

    @property
    def t(self):
        return self._t

    @t.setter
    def t(self, value):
        self._t = value

    @property
    def bf(self):
        return self._bf

    @bf.setter
    def bf(self, value):
        self._bf = value

    @property
    def fi(self):
        return self._fi

    @fi.setter
    def fi(self, value):
        self._fi = value

    @property
    def mat(self):
        return self._mat

    @mat.setter
    def mat(self, value):
        self._mat = value

    def getS1_Hat(self, corr_add: CorrosionAddition):  # Net width of Hat profile webs, outside corners at the connection with attached plating - S1, [mm]
        h = self._h + corr_add.tc       # Net profile height
        t = self._t - corr_add.tc       # Net flange and web thickness
        alpha = (180 - self._fi) / 2
        bf = self._bf - 2 * corr_add.tc / (2 * np.tan(np.radians(alpha)))        # Net flange width
        fi = self._fi
        S1 = bf + ((2 * (h + t)) / np.tan(np.radians(fi)))
        return S1

    def getShArea_Hat(self, corr_add: CorrosionAddition):  # Net shear area of Hat profile - Ash, [cm2]
        h = self.h + corr_add.tc
        t = self.t - corr_add.tc
        fi = self.fi
        Ash = 2 * h * (t / np.sin(np.radians(fi))) * 0.01
        return Ash

    @property
    def getArea(self):  # Gross cross sectional area of Hat profile, without attached plating - A, [cm2]
        h = self.h
        t = self.t
        bf = self.bf
        fi = self.fi
        A = (((2 * h * t) / np.sin(np.radians(fi))) + (t ** 2 / np.tan(np.radians(fi))) + (bf * t)) * 0.01
        return A

    def getArea_I(self, bp, tp_gross, corr_add: CorrosionAddition):  # Net area of Hat stiffener with attached plating - A, [cm2]
        h = self._h + corr_add.tc       # Net profile height
        t = self._t - corr_add.tc       # Net flange and web thickness
        alpha = (180 - self._fi) / 2
        bf = self._bf - 2 * corr_add.tc / (2 * np.tan(np.radians(alpha)))        # Net flange width
        tp_net = tp_gross - corr_add.tc
        fi = self._fi
        A = ((2 * h * t) / np.sin(np.radians(fi)) + t ** 2 / np.tan(np.radians(fi)) + bf * t + bp * tp_net) * 0.01
        return A

    def get_z_na_Hat(self, bp, tp_gross, corr_add: CorrosionAddition):
        h = self._h + corr_add.tc       # Net profile height
        t = self._t - corr_add.tc       # Net flange and web thickness
        fi = self._fi
        tp = tp_gross - corr_add.tc
        A = self.getArea_I(bp, tp_gross, corr_add)
        zna = ((h * t / np.sin(np.radians(fi))) * (h + t) + t ** 3 / (6 * np.tan(np.radians(fi))) +
               (h + (tp + t) / 2) * bp * tp) / (A * 100)
        return zna

    def get_Iy_I(self, bp, tp_gross, corr_add: CorrosionAddition):  # Moment of inertia around the local y axis - Iy, [cm4]
        h = self._h + corr_add.tc       # Net profile height
        t = self._t - corr_add.tc       # Net flange and web thickness
        alpha = (180 - self._fi) / 2
        bf = self._bf - 2 * corr_add.tc / (2 * np.tan(np.radians(alpha)))        # Net flange width
        fi = self._fi
        tp = tp_gross - corr_add.tc
        zna = self.get_z_na_Hat(bp, tp_gross, corr_add)
        Iy = (2 * ((t * h ** 3 / (12 * np.sin(np.radians(fi)))) + (((t + h) / 2 - zna) ** 2) * (h * t / np.sin(np.radians(fi)))) +
              2 * ((t ** 4 / (36 * np.tan(np.radians(fi)))) + (((zna - t / 6) ** 2) * (t ** 2 / (2 * np.tan(np.radians(fi)))))) +
              ((bf * t ** 3 / 12) + zna ** 2 * bf * t) + ((bp * tp ** 3 / 12) + ((h + ((tp + t) / 2 - zna)) ** 2 * bp * tp)))
        return Iy * (10 ** (-4))

    def get_Ip(self, corr_add: CorrosionAddition):  # Polar moment of inertia in [cm4]
        h = self._h + corr_add.tc       # Net profile height
        t = self._t - corr_add.tc       # Net flange and web thickness
        alpha = (180 - self._fi) / 2
        bf = self._bf - 2 * corr_add.tc / (2 * np.tan(np.radians(alpha)))        # Net flange width
        Ip = ((((h ** 3) * t) / 3) + ((h ** 2) * bf * t)) * 10 ** (-4)
        return Ip

    def get_It(self, corr_add: CorrosionAddition):  # St Venant's moment of inertia in [cm4]
        h = self._h + corr_add.tc       # Net profile height
        t = self._t - corr_add.tc       # Net flange and web thickness
        alpha = (180 - self._fi) / 2
        bf = self._bf - 2 * corr_add.tc / (2 * np.tan(np.radians(alpha)))        # Net flange width
        It = (1 / 3) * ((h * t ** 3) + (bf * t ** 3) * (1 - 0.63 * (t / bf))) * 10 ** (-4)
        return It

    def getWmin(self, bp, tp, corr_add: CorrosionAddition):
        Iy = self.get_Iy_I(bp, tp, corr_add)
        zna = self.get_z_na_Hat(bp, tp, corr_add)
        t = self._t - corr_add.tc                   # Net flange and web thickness
        h = self._h + corr_add.tc                   # Net profile height
        wf = Iy / ((zna + (t / 2)) * 0.1)           # Section modulus of the flange - Wf, [cm3]
        wp = Iy / ((tp + h + t * 0.5 - zna) * 0.1)  # Section modulus of the plating - Wp, [cm3]
        if wp < wf:
            wmin = wp
            return wmin
        elif wf <= wp:
            wmin = wf
            return wmin


class VariableSection(BeamProperty):
    def __init__(self, id_, hw_min, hw_max, hw_l1, tw, bf_min, bf_max, bf_l1, tf, mat: MaterialProperty):
        super().__init__(id_)
        self._hw_min = float(hw_min)        # Minimum web height
        self._hw_max = float(hw_max)        # Maximum web height
        self._hw_l1 = float(hw_l1)          # Distance of web height change from the end with dimension hw_min
        self._tw = float(tw)                # Web thickness
        self._bf_min = float(bf_min)        # Minimum flange width
        self._bf_max = float(bf_max)        # Maximum flange width
        self._bf_l1 = float(bf_l1)          # Distance of flange width change from the end with dimension bf_min
        self._tf = float(tf)                # Flange thickness
        self._mat = mat                     # Material

    # dodati:
    #   property i setteri
    #   metode za izračun točaka na bilo kojem mjestu duž prirubnice i struka - za položaj čvorova


class PrimarySuppMem:
    def __init__(self, idbeam, direction: BeamDirection, rel_dist, grillage):
        self._id = idbeam
        self._segments = []                             # List of segments for each primary supporting member
        self._direction: BeamDirection = direction      # Direction of the primary supporting member
        self._flange_direction: FlangeDirection = grillage.flange_direction
        self._rel_dist = float(rel_dist)                # Relative distance - position of the primary supporting member
        self._grillage = grillage
        self._symmetric_member = None

    def generate_segment(self, beamprop):
        id_segment = 1
        cross_beams = self.grillage.get_oposite_dir_beams(self._direction)
        for i_cross in range(1, len(cross_beams)):
            cross1 = cross_beams[i_cross]
            cross2 = cross_beams[i_cross + 1]
            current_segment = Segment(id_segment, beamprop, self, cross1, cross2)
            self._segments.append(current_segment)
            id_segment += 1

    @property
    def flange_direction(self):
        """
        :return: Flange direction unit vector, based on Primary supporting member direction and flange_direction variable
            set for the entire grillage model as INWARD (towards the centerline), or OUTWARD (opposite of centerline).
            Default value of flange_direction is set as INWARD to ensure edge L beams are oriented correctly.
        """
        flange_dir = self._flange_direction

        if self.direction == BeamDirection.LONGITUDINAL:
            if (self.rel_dist <= 0.5 and flange_dir == FlangeDirection.INWARD) \
                    or (self.rel_dist > 0.5 and flange_dir == FlangeDirection.OUTWARD):
                return np.array((0, 1, 0))

            elif (self.rel_dist <= 0.5 and flange_dir == FlangeDirection.OUTWARD) \
                    or (self.rel_dist > 0.5 and flange_dir == FlangeDirection.INWARD):
                return np.array((0, -1, 0))

        if self.direction == BeamDirection.TRANSVERSE:
            if (self.rel_dist <= 0.5 and flange_dir == FlangeDirection.INWARD) \
                    or (self.rel_dist > 0.5 and flange_dir == FlangeDirection.OUTWARD):
                return np.array((1, 0, 0))

            elif (self.rel_dist <= 0.5 and flange_dir == FlangeDirection.OUTWARD) \
                    or (self.rel_dist > 0.5 and flange_dir == FlangeDirection.INWARD):
                return np.array((-1, 0, 0))

    @flange_direction.setter
    def flange_direction(self, value):
        self._flange_direction = value

    @property
    def symmetric_member(self):
        return self._symmetric_member

    @symmetric_member.setter
    def symmetric_member(self, value):
        if self._symmetric_member is None:
            self._symmetric_member = value
            self._symmetric_member.symmetric_member = self

    @property
    def has_symmetric_memb(self):
        return isinstance(self._symmetric_member, PrimarySuppMem)

    def set_rel_dist_symmetric(self, value):
        self._rel_dist = value
        self._symmetric_member._rel_dist = 1.0 - value

    @property
    def end_nodes(self):
        return Segment.end_nodes(self)

    @property
    def grillage(self):
        return self._grillage

    @property
    def id(self):
        return self._id

    @property
    def segments(self):
        return self._segments

    @segments.setter
    def segments(self, value):
        self._segments = value

    @property
    def direction(self):
        return self._direction

    @property
    def rel_dist(self):
        return self._rel_dist

    @rel_dist.setter
    def rel_dist(self, value):
        self._rel_dist = value


class Segment:
    def __init__(self, idsegment, beam_prop, primary_supp_mem, cross_member1, cross_member2):
        self._id = idsegment
        self._beam_prop: BeamProperty = beam_prop
        self._primary_supp_mem = primary_supp_mem
        self._cross_member1 = cross_member1
        self._cross_member2 = cross_member2
        self._symmetric_segment = None

    @property
    def id(self):
        return self._id

    @property
    def beam_prop(self):
        return self._beam_prop

    @beam_prop.setter
    def beam_prop(self, value):
        self._beam_prop = value

    @property
    def primary_supp_mem(self):
        return self._primary_supp_mem

    @primary_supp_mem.setter
    def primary_supp_mem(self, value):
        self._primary_supp_mem = value

    @property
    def cross_member1(self):
        return self._cross_member1

    @cross_member1.setter
    def cross_member1(self, value):
        self._cross_member1 = value

    @property
    def cross_member2(self):
        return self._cross_member2

    @cross_member2.setter
    def cross_member2(self, value):
        self._cross_member2 = value

    @property
    def symmetric_segment(self):
        return self._symmetric_segment

    @symmetric_segment.setter
    def symmetric_segment(self, value):
        if self._symmetric_segment is None:
            self._symmetric_segment = value
            self._symmetric_segment.symmetric_segment = self

    @staticmethod
    def end_nodes(member: PrimarySuppMem):
        """
        :param member: Primary supporting member.
        :return: Primary supporting member end node coordinates in [mm], at the point of connection with plating.
                Flange of the primary supporting member is at vertical coordinate z = 0.
        """
        grillage = member.grillage
        hw_end1 = member.segments[0].beam_prop.hw                           # Primary supporting member web height at x or y = 0
        hw_end2 = member.segments[len(member.segments) - 1].beam_prop.hw    # Primary supporting member web height at x = L or y = B

        if member.direction == BeamDirection.LONGITUDINAL:  # Longitudinal primary supporting members
            node1 = Node(1, 0, member.rel_dist * grillage.B_overall * 1000, hw_end1)
            node2 = Node(2, grillage.L_overall * 1000, member.rel_dist * grillage.B_overall * 1000, hw_end2)
            return node1.coords, node2.coords

        if member.direction == BeamDirection.TRANSVERSE:  # Transverse primary supporting members
            node1 = Node(1, member.rel_dist * grillage.L_overall * 1000, 0, hw_end1)
            node2 = Node(2, member.rel_dist * grillage.L_overall * 1000, grillage.B_overall * 1000, hw_end2)
            return node1.coords, node2.coords

    def get_segment_node1(self):
        return self.get_segment_node(self._cross_member1)

    def get_segment_node2(self):
        return self.get_segment_node(self._cross_member2)

    def get_segment_node(self, cross_member):
        (n1, n2) = self._primary_supp_mem.end_nodes
        (n3, n4) = cross_member.end_nodes
        node = self._primary_supp_mem.grillage.get_intersection(n1, n2, n3, n4)
        return node

    def segment_len(self):
        # Segment length
        rel_dist1 = self._cross_member1.rel_dist
        rel_dist2 = self._cross_member2.rel_dist
        if self._primary_supp_mem.direction == BeamDirection.LONGITUDINAL:
            return (rel_dist2 - rel_dist1) * self._primary_supp_mem.grillage.L_overall
        if self._primary_supp_mem.direction == BeamDirection.TRANSVERSE:
            return (rel_dist2 - rel_dist1) * self._primary_supp_mem.grillage.B_overall

    def get_attplate(self):
        # Segment attached plating
        tps = []  # List of plating thickness of attached plate (1 or 2 values)
        spn = []  # List of half distances between between adjacent primary supporting members (1 or 2 values)
        att_plating_n = 0

        for plate in self._primary_supp_mem.grillage.plating().values():
            if plate.test_plate_segment(self):
                tps.append(plate.plate_prop.tp)

                if self._primary_supp_mem.direction == BeamDirection.TRANSVERSE:
                    spn.append(0.5 * Plate.plate_longitudinal_dim(plate))
                if self._primary_supp_mem.direction == BeamDirection.LONGITUDINAL:
                    spn.append(0.5 * Plate.plate_transverse_dim(plate))

                att_plating_n += 1

                # Stop searching through plating zones if all associated plate zones are identified
                if self._primary_supp_mem.rel_dist == 1 or self._primary_supp_mem.rel_dist == 0 and att_plating_n == 1:
                    break  # Maximum number of associated plate zones for edge primary supporting members is 1
                elif att_plating_n == 2:
                    break  # Maximum number of associated plate zones is 2

        tpsr_seg = 0
        for tp in tps:
            tpsr_seg += tp                  # Average value of attached plating thickness
        tpsr_seg = tpsr_seg / len(tps)      # Attached plate thickness

        bp_seg = 0
        bp = 0
        for sp in spn:
            if self._primary_supp_mem.direction == BeamDirection.LONGITUDINAL:
                bp = np.minimum(0.165 * self._primary_supp_mem.grillage.L_overall, sp)
            bp_seg += bp                # Attached plate width for longitudinal segments
            if self._primary_supp_mem.direction == BeamDirection.TRANSVERSE:
                bp = np.minimum(0.165 * self._primary_supp_mem.grillage.B_overall, sp)
                bp_seg += bp            # Attached plate width for transverse segments
        return bp_seg * 1000, tpsr_seg  # Returns attached plate width and thickness, both in [mm]

    @property
    def Wmin(self):
        """
        :return: Minimum net section modulus of the given Segment. Corrosion addition is stored as
                the first input value (ID = 1) in dictionary grillgeo.corrosion_addition
        """
        (bp, tp) = self.get_attplate()
        tc = self._primary_supp_mem.grillage.corrosion_addition()[1]
        return self._beam_prop.getWmin(bp, tp, tc)

    @property
    def Iy(self):
        (bp, tp) = self.get_attplate()
        tc = self._primary_supp_mem.grillage.corrosion_addition()[1]
        return self._beam_prop.get_Iy_I(bp, tp, tc)


class StiffenerLayout:
    def __init__(self, id_layout, beam_prop, definition_type, definition_value):
        self._id = id_layout
        self._beam_prop = beam_prop                 # Stiffener beam property
        self._definition_type = definition_type     # Type of stiffener definition: spacing or number of stiffeners between primary supporting members
        self._definition_value = float(definition_value)   # Value for chosen stiffener definition: number or spacing

    @property
    def id(self):
        return self._id

    @property
    def beam_prop(self):
        return self._beam_prop

    @beam_prop.setter
    def beam_prop(self, value):
        self._beam_prop = value

    @property
    def definition_type(self):
        return self._definition_type

    @definition_type.setter
    def definition_type(self, value):
        self._definition_type = value

    @property
    def definition_value(self):
        return self._definition_value

    @definition_value.setter
    def definition_value(self, value):
        self._definition_value = value


class Plate:
    def __init__(self, idplate, plate_prop, long_seg1, trans_seg1, long_seg2, trans_seg2, stiff_layout, stiff_dir, ref_edge):
        self._id = idplate
        self._plate_prop = plate_prop
        self._stiff_layout = stiff_layout
        self._stiff_dir: BeamDirection = stiff_dir
        self._ref_edge: Ref = ref_edge      # Starting segment for layouts defined by spacing between stiffeners
        self._segments = [long_seg1, trans_seg1, long_seg2, trans_seg2]
        self._symmetric_plate_zones = []

    def test_plate_segment(self, test_segment: Segment):
        # Segment association test
        for seg in range(4):
            if self._segments[seg] is test_segment:
                return True     # Returns true if test_segment defines the plating zone
        return False

    def test_plate_psm(self, test_member: PrimarySuppMem):
        # Primary supporting member association test
        for seg in range(4):
            if self._segments[seg].primary_supp_mem is test_member:
                return True     # Returns true if the plate zone is located along the test_member primary supporting member
        return False

    def test_plate_between_psm(self, test_member1: PrimarySuppMem, test_member2: PrimarySuppMem):
        # Primary supporting member association test
        if test_member1.direction and test_member2.direction == BeamDirection.LONGITUDINAL:
            if self.long_seg1.primary_supp_mem is test_member1 and self.long_seg2.primary_supp_mem is test_member2:
                return True
            return False
        elif test_member1.direction and test_member2.direction == BeamDirection.TRANSVERSE:
            if self.trans_seg1.primary_supp_mem is test_member1 and self.trans_seg2.primary_supp_mem is test_member2:
                return True
            return False    # Returns true if the plate zone is located between two adjacent primary supporting members

    def plate_longitudinal_dim(self):
        longitudinal_dim = Segment.segment_len(self._segments[0])
        return longitudinal_dim

    def plate_transverse_dim(self):
        transverse_dim = Segment.segment_len(self._segments[1])
        return transverse_dim

    def plate_dim_parallel_to_stiffeners(self):
        if self.stiff_dir == BeamDirection.LONGITUDINAL:
            return self.plate_longitudinal_dim()
        elif self.stiff_dir == BeamDirection.TRANSVERSE:
            return self.plate_transverse_dim()

    def get_reference_segment(self):
        # Reference segment for stiffener placement
        if self._stiff_dir == BeamDirection.LONGITUDINAL:
            if self._ref_edge == Ref.EDGE1:
                return self._segments[0]
            elif self._ref_edge == Ref.EDGE2:
                return self._segments[2]
        if self._stiff_dir == BeamDirection.TRANSVERSE:
            if self._ref_edge == Ref.EDGE1:
                return self._segments[1]
            elif self._ref_edge == Ref.EDGE2:
                return self._segments[3]

    def get_path_length(self):
        # Returns the stiffener path length based on orientation - dimension of the plating zone perpendicular to stiffener direction
        if self._stiff_dir == BeamDirection.LONGITUDINAL:
            return self.plate_transverse_dim()
        if self._stiff_dir == BeamDirection.TRANSVERSE:
            return self.plate_longitudinal_dim()

    def get_stiffener_spacing(self):  # Returns stiffener spacing on any plating zone regardless of definition type
        if self.stiff_layout.definition_type == "number":
            return self.get_path_length() / (self.stiff_layout.definition_value + 1)
        elif self.stiff_layout.definition_type == "spacing":
            return self.stiff_layout.definition_value

    def get_stiffener_number(self):  # Returns number of stiffeners on any plating zone regardless of definition type
        if self.stiff_layout.definition_type == "number":
            return self.stiff_layout.definition_value
        elif self.stiff_layout.definition_type == "spacing":
            return np.ceil((np.round(self.get_path_length(), decimals=6) / self.stiff_layout.definition_value) - 1)

    def get_equal_stiffener_offset(self):
        stiff_num = self.get_stiffener_number()
        stiff_spacing = self.get_stiffener_spacing()
        path_len = self.get_path_length()
        equal_offset = (path_len - ((stiff_num - 1) * stiff_spacing)) / 2
        return equal_offset

    def get_stiff_coords(self, stiffener_n: int):
        """
        :param stiffener_n: n-th stiffener on the plating zone. Enter a integer value from 1 to number of stiffeners on the plating zone.
        :return: Coordinates of the n-th stiffener on the plating zone in [mm], at the point of connection with plating.
        """
        perpendicular_vector = np.array((0, 0, 0))
        ref_node1 = Segment.get_segment_node1(self.get_reference_segment())     # Reference node 1 coordinates in [mm]
        ref_node2 = Segment.get_segment_node2(self.get_reference_segment())     # Reference node 2 coordinates in [mm]
        ref_vector = ref_node2 - ref_node1                                      # Reference vector in the direction of the reference segment
        ref_vector_magnitude = np.linalg.norm(ref_vector)                       # Reference vector length
        unit_ref_vector = ref_vector / ref_vector_magnitude                     # Unit reference vector
        normal_vector = np.array((0, 0, 1))                                     # Vector normal to the plating surface
        # U slucaju oplate sa prelukom, vektor normale nije (0,0,1) !!!
        if self.stiff_dir == BeamDirection.LONGITUDINAL:
            perpendicular_vector = np.cross(normal_vector, unit_ref_vector)
        elif self.stiff_dir == BeamDirection.TRANSVERSE:
            perpendicular_vector = np.cross(unit_ref_vector, normal_vector)

        stiff_spacing = self.get_stiffener_spacing() * 1000                     # Stiffener spacing in [mm]
        spacing_vector = perpendicular_vector * stiff_spacing                   # Vector with magnitude of stiffener spacing
        stiff_offset = self.get_equal_stiffener_offset() * 1000                 # First and last stiffener offset in [mm]
        offset_vector = perpendicular_vector * stiff_offset                     # Vector with magnitude of stiffener offset

        if stiffener_n == 1:
            stiff_n_node1 = ref_node1 + offset_vector
        else:
            stiff_n_node1 = ref_node1 + offset_vector + spacing_vector * (stiffener_n - 1)
        stiff_n_node2 = stiff_n_node1 + ref_vector

        return stiff_n_node1, stiff_n_node2

    def Wmin(self):
        """
        :return: Minimum net section modulus of stiffeners on the given plating zone. Corrosion addition is stored as
                the first input value (ID = 1) in dictionary grillgeo.corrosion_addition
        """
        stiff_property = self.stiff_layout.beam_prop
        bp = self.get_stiffener_spacing() * 1000                                        # Stiffener spacing in [mm]
        tp = self.plate_prop.tp                                                         # Plating thickness in [mm]
        corr_add = self.long_seg1.primary_supp_mem.grillage.corrosion_addition()[1]     # CorrosionAddition object
        return stiff_property.getWmin(bp, tp, corr_add)

    def Iy(self):
        stiff_property = self.stiff_layout.beam_prop
        bp = self.get_stiffener_spacing() * 1000                                        # Stiffener spacing in [mm]
        tp = self.plate_prop.tp                                                         # Plating thickness in [mm]
        corr_add = self.long_seg1.primary_supp_mem.grillage.corrosion_addition()[1]     # CorrosionAddition object
        return stiff_property.get_Iy_I(bp, tp, corr_add)

    def tp_net(self):
        tp_gross = self.plate_prop.tp
        corr_add = self.long_seg1.primary_supp_mem.grillage.corrosion_addition()[1]
        tp_net = PlateProperty.tp_net(corr_add, tp_gross)
        return tp_net

    # Symmetric plating zones
    @property
    def symmetric_plate_zones(self):
        return self._symmetric_plate_zones

    @symmetric_plate_zones.setter
    def symmetric_plate_zones(self, value):
        if len(self._symmetric_plate_zones) == 0:
            self._symmetric_plate_zones.append(value)
        elif len(self._symmetric_plate_zones) < 3 and value is not self and value not in self._symmetric_plate_zones:
            self._symmetric_plate_zones.append(value)

    @property
    def long_seg1(self):
        return self._segments[0]

    @property
    def trans_seg1(self):
        return self._segments[1]

    @property
    def long_seg2(self):
        return self._segments[2]

    @property
    def trans_seg2(self):
        return self._segments[3]

    @property
    def id(self):
        return self._id

    @property
    def plate_prop(self):
        return self._plate_prop

    @plate_prop.setter
    def plate_prop(self, value):
        self._plate_prop = value

    @property
    def stiff_layout(self):
        return self._stiff_layout

    @stiff_layout.setter
    def stiff_layout(self, value):
        self._stiff_layout = value

    @property
    def stiff_dir(self):
        return self._stiff_dir

    @stiff_dir.setter
    def stiff_dir(self, value):
        self._stiff_dir = value

    @property
    def ref_edge(self):
        return self._ref_edge

    @ref_edge.setter
    def ref_edge(self, value):
        self._ref_edge = value

    @property
    def plate_segments(self):
        return self._segments


class Grillage():
    def __init__(self, L_overall, B_overall, N_longitudinal, N_transverse):
        self._L_overall = L_overall                 # Overall length, m
        self._B_overall = B_overall                 # Overall width, m
        self._N_longitudinal = N_longitudinal       # Number of longitudinal primary supporting members
        self._N_transverse = N_transverse           # Number of transverse primary supporting members
        self._flange_direction = FlangeDirection.INWARD   # Orientation of L primary supporting member flange; default = INWARD

        self._longitudinal_memb = {}
        self._transverse_memb = {}
        self._material_props = {}
        self._beam_props = {}
        self._plate_props = {}
        self._stiffener_layouts = {}
        self._plating = {}
        self._corrosion_add = {}

    @property
    def L_overall(self):
        return self._L_overall

    @L_overall.setter
    def L_overall(self, value):
        self._L_overall = value

    @property
    def B_overall(self):
        return self._B_overall

    @B_overall.setter
    def B_overall(self, value):
        self._B_overall = value

    @property
    def N_longitudinal(self):
        return self._N_longitudinal

    @N_longitudinal.setter
    def N_longitudinal(self, value):
        self._N_longitudinal = value

    @property
    def N_transverse(self):
        return self._N_transverse

    @N_transverse.setter
    def N_transverse(self, value):
        self._N_transverse = value

    @property
    def flange_direction(self):
        return self._flange_direction

    @flange_direction.setter
    def flange_direction(self, value):
        self._flange_direction = value

    def longitudinal_members(self):
        return self._longitudinal_memb

    def transverse_members(self):
        return self._transverse_memb

    def material_props(self):
        return self._material_props

    def add_material(self, mat: MaterialProperty):
        self._material_props[mat.id] = mat

    def beam_props(self):
        return self._beam_props

    def add_beam_prop(self, beam: BeamProperty):
        self._beam_props[beam.id] = beam

    def plate_props(self):
        return self._plate_props

    def add_plate_prop(self, plate: PlateProperty):
        self._plate_props[plate.id] = plate

    def stiffener_layouts(self):
        return self._stiffener_layouts

    def add_stiffener_layout(self, layout: StiffenerLayout):
        self._stiffener_layouts[layout.id] = layout

    def plating(self):
        return self._plating

    def add_plating(self, plate: Plate):
        self._plating[plate.id] = plate

    def corrosion_addition(self):
        return self._corrosion_add

    def add_corrosion_addition(self, value):
        self._corrosion_add[value.id] = value

    @staticmethod
    def get_intersection(node1, node2, node3, node4):
        """
        :param node1: List of [x, y, z] coordinates of the first point defining the first girder.
        :param node2: List of [x, y, z] coordinates of the second point defining the first girder.
        :param node3: List of [x, y, z] coordinates of the first point defining the second girder.
        :param node4: List of [x, y, z] coordinates of the second point defining the second girder.
        :return: Coordinates of the intersection of the first and second girder, type numpy.ndarray.
        """
        vector1 = node2 - node1                 # First girder defined by node1 and node2
        vector2 = node4 - node3                 # Second girder defined by node3 and node4
        magnitude1 = np.linalg.norm(vector1)    # Length of the first girder
        magnitude2 = np.linalg.norm(vector2)    # Length of the second girder
        direction1 = vector1 / magnitude1       # Direction vector of the first girder
        direction2 = vector2 / magnitude2       # Direction vector of the second girder

        a = ([[direction1[0], -direction2[0]], [direction1[1], -direction2[1]]])
        b = ([[node3[0] - node1[0]], [node3[1] - node1[1]]])
        det = np.linalg.det(a)
        if det == 0:  # Determinant of matrix a is 0 if there is no intersection
            pass
        else:
            solve = np.hstack(np.linalg.solve(a, b))
            x_inter = solve[0] * direction1[0] + node1[0]
            y_inter = solve[0] * direction1[1] + node1[1]
            z_inter = solve[0] * direction1[2] + node1[2]
            coords = np.array([x_inter, y_inter, z_inter])  # x, y, z coordinates of the intersection point
            return coords

    @staticmethod
    def get_midpoint(node1, node2):
        """
        :param node1: List of [x, y, z] coordinates for the first point.
        :param node2: List of [x, y, z] coordinates for the second point.
        :return: Coordinates of the midpoint between node1 and node2, type numpy.ndarray.
        """
        vector = np.subtract(node2, node1)
        unit_vector = vector / np.linalg.norm(vector)
        half_dist = np.linalg.norm(vector) / 2
        midpoint = unit_vector * half_dist + node1
        return midpoint

    def generate_prim_supp_members(self):
        # Generate primary supporting members
        #       Generates initial nodes and assigns primary supporting members
        #       based on basic data input values for uniform girder spacing
        delta_xrel = 1 / (self._N_transverse - 1)       # Initial uniform relative spacing of transverse primary supporting members
        delta_yrel = 1 / (self._N_longitudinal - 1)     # Initial uniform relative spacing of longitudinal primary supporting members

        # Generate longitudinal primary supporting members
        for i_long in range(1, self._N_longitudinal + 1):
            current_member = PrimarySuppMem(i_long, BeamDirection.LONGITUDINAL, (i_long - 1) * delta_yrel, self)
            self._longitudinal_memb[i_long] = current_member

        # Generate transverse primary supporting members
        for i_trans in range(1, self._N_transverse + 1):
            current_member = PrimarySuppMem(i_trans, BeamDirection.TRANSVERSE, (i_trans - 1) * delta_xrel, self)
            self._transverse_memb[i_trans] = current_member

    def generate_segments(self, beamprop_longitudinal, beamprop_transverse, beamprop_edge):
        for beam in self._longitudinal_memb.values():
            if beam.rel_dist == 0.0 or beam.rel_dist == 1.0:
                beam.generate_segment(beamprop_edge)
            else:
                beam.generate_segment(beamprop_longitudinal)

        for beam in self._transverse_memb.values():
            if beam.rel_dist == 0.0 or beam.rel_dist == 1.0:
                beam.generate_segment(beamprop_edge)
            else:
                beam.generate_segment(beamprop_transverse)

    def get_oposite_dir_beams(self, direction):
        if direction == BeamDirection.LONGITUDINAL:  # Longitudinal primary supporting members
            return self._transverse_memb
        elif direction == BeamDirection.TRANSVERSE:
            return self._longitudinal_memb
        else:
            print('Unkonown beam direction', direction)

    def generate_plating(self, plateprop1, stifflayout1, stiff_dir, ref_edge=Ref.EDGE1):
        # Plating zone defined by longitudinal and transverse segments along its edge
        plate_id = 1
        for i_long in range(1, len(self._longitudinal_memb)):
            for i_segment in range(0, self._N_transverse - 1):
                long_seg1 = self._longitudinal_memb[i_long].segments[i_segment]
                long_seg2 = self._longitudinal_memb[i_long + 1].segments[i_segment]
                transversal1 = self._longitudinal_memb[i_long].segments[i_segment].cross_member1.id
                transversal2 = self._longitudinal_memb[i_long].segments[i_segment].cross_member2.id
                trans_seg1 = self._transverse_memb[transversal1].segments[i_long - 1]
                trans_seg2 = self._transverse_memb[transversal2].segments[i_long - 1]
                curr_plate = Plate(plate_id, plateprop1, long_seg1, trans_seg1, long_seg2, trans_seg2, stifflayout1,
                                   stiff_dir, ref_edge)
                self._plating[curr_plate.id] = curr_plate
                plate_id += 1

    @staticmethod
    def get_segment_nodes(segment):
        node1 = Segment.get_segment_node1(segment)
        node2 = Segment.get_segment_node2(segment)
        return node1, node2

    @staticmethod
    def get_segment_length(segment):
        length = Segment.segment_len(segment)
        return length

    def assign_symmetric_members(self):
        # Assign symmetric primary supporting members
        for i_long in self._longitudinal_memb.values():
            if i_long.rel_dist < 0.5:
                symmetric_id = len(self._longitudinal_memb) - i_long.id + 1     # Longitudinal axis of symmetry
                i_long.symmetric_member = self._longitudinal_memb[symmetric_id]

        for i_trans in self._transverse_memb.values():
            if i_trans.rel_dist < 0.5:
                symmetric_id = len(self._transverse_memb) - i_trans.id + 1      # Transverse axis of symmetry
                i_trans.symmetric_member = self._transverse_memb[symmetric_id]

    def assign_symmetric_plating(self):
        # Assign symmetric plating zones
        for id_plate in range(1, (self._N_longitudinal - 1) * (self._N_transverse - 1) + 1):
            plate_id_list = np.arange(1, len(self._plating) + 1)
            plate_id_matrix = np.reshape(plate_id_list, (self._N_longitudinal - 1, self._N_transverse - 1))
            ir = np.where(plate_id_matrix == id_plate)[0]       # Row index of selected plating zone in id_matrix
            ic = np.where(plate_id_matrix == id_plate)[1]       # Column index of selected plating zone in id_matrix

            matrix_longsym = np.flip(plate_id_matrix, axis=0)
            matrix_transsym = np.flip(plate_id_matrix, axis=1)
            matrix_combsym = np.flip(matrix_transsym, axis=0)

            id_plate_longsym = int(matrix_longsym[ir, ic])      # Longitudinal axis of symmetry condition
            id_plate_transsym = int(matrix_transsym[ir, ic])    # Transverse axis of symmetry condition
            id_plate_combsym = int(matrix_combsym[ir, ic])      # Combined longitudinal and transverse symmetry conditions

            self._plating[id_plate].symmetric_plate_zones = self._plating[id_plate_longsym]
            self._plating[id_plate].symmetric_plate_zones = self._plating[id_plate_transsym]
            self._plating[id_plate].symmetric_plate_zones = self._plating[id_plate_combsym]

    def getGrillageMass(self):
        # Longitudinal segments mass
        mass_long = 0.0
        for i_long in self._longitudinal_memb.keys():
            for i_segment in range(0, self._N_transverse - 1):
                curr_segment = self._longitudinal_memb[i_long].segments[i_segment]
                area = curr_segment.beam_prop.getArea * 0.0001   # Beam area, [m2]
                density = curr_segment.beam_prop.mat.ro          # Material density, [kg/m3]
                length = Segment.segment_len(curr_segment)       # Segment length, [m]
                mass = area * length * density
                mass_long += mass

        # Transverse segments mass
        mass_trans = 0.0
        for i_tran in self._transverse_memb.keys():
            for i_segment in range(0, self._N_longitudinal - 1):
                curr_segment = self._transverse_memb[i_tran].segments[i_segment]
                area = curr_segment.beam_prop.getArea * 0.0001                    # Beam area, [m2]
                density = curr_segment.beam_prop.mat.ro                           # Material density, [kg/m3]
                length = Segment.segment_len(curr_segment)                        # Segment length, [m]
                mass = area * length * density
                mass_trans += mass

        # Plating mass
        mass_plate = 0.0
        for id_plate in range(1, (self._N_longitudinal - 1) * (self._N_transverse - 1) + 1):
            curr_plate = self._plating[id_plate]
            tp = curr_plate.plate_prop.tp * 0.001               # Plate thickness, [m]
            lPlate = Plate.plate_longitudinal_dim(curr_plate)   # Longitudinal plate dimension, [m]
            bPlate = Plate.plate_transverse_dim(curr_plate)     # Transverse plate dimension, [m]
            density = curr_plate.plate_prop.plate_mat.ro        # Plate material density, [kg/m3]
            mass = tp * lPlate * bPlate * density               # Mass of selected plate, [kg]
            mass_plate += mass

        # Stiffener mass
        mass_stiff = 0.0
        for id_plate in range(1, (self._N_longitudinal - 1) * (self._N_transverse - 1) + 1):
            curr_plate = self._plating[id_plate]
            stiff_num = Plate.get_stiffener_number(curr_plate)
            stiff_leng = self.get_segment_length(Plate.get_reference_segment(curr_plate))
            density = curr_plate.stiff_layout.beam_prop.mat.ro
            area = curr_plate.stiff_layout.beam_prop.getArea * 0.0001
            mass = stiff_num * stiff_leng * area * density
            mass_stiff += mass

        return mass_long + mass_trans + mass_plate + mass_stiff

    def get_elementary_plate_dimensions(self, plate_id):
        plate = self.plating()[plate_id]
        stiff_direction = plate.stiff_dir
        stiff_spacing = plate.get_stiffener_spacing()
        dim_x = plate.plate_longitudinal_dim()
        dim_y = plate.plate_transverse_dim()

        if stiff_direction == BeamDirection.LONGITUDINAL:
            if plate.stiff_layout.beam_prop is HatBeamProperty:
                if dim_x > stiff_spacing:
                    ss = stiff_spacing - 0.5 * HatBeamProperty.getS1_Hat(plate.stiff_layout.beam_prop, self.corrosion_addition()[1])
                    ls = dim_x
                else:
                    ss = dim_x
                    ls = stiff_spacing
            else:
                if dim_x > stiff_spacing:
                    ss = stiff_spacing
                    ls = dim_x
                else:
                    ss = dim_x
                    ls = stiff_spacing
            return ss, ls   # Returns the length, in [m], of the shorter and longer side of the plate panel

        elif stiff_direction == BeamDirection.TRANSVERSE:
            if plate.stiff_layout.beam_prop is HatBeamProperty:
                if dim_y > stiff_spacing:
                    ss = stiff_spacing - 0.5 * HatBeamProperty.getS1_Hat(plate.stiff_layout.beam_prop, self.corrosion_addition()[1])
                    ls = dim_y
                else:
                    ss = dim_y
                    ls = stiff_spacing
            else:
                if dim_y > stiff_spacing:
                    ss = stiff_spacing
                    ls = dim_y
                else:
                    ss = dim_y
                    ls = stiff_spacing
            return ss, ls   # Returns the length, in [m], of the shorter and longer side of the plate panel

    @staticmethod
    def plate_common_segment(plate1: Plate, plate2: Plate):
        # Returns the common segment of two given plating zones
        segments_list1 = plate1.plate_segments
        segments_list2 = plate2.plate_segments
        common_segment = list(set(segments_list1).intersection(segments_list2))
        if not common_segment:
            return None
        else:
            return common_segment[0]

    def get_adjacent_plates(self, test_plate: Plate):
        # Returns a list of Plate objects which have a common segment with test_plate
        adjacent_plate_list = []
        for plate_id in self._plating.keys():
            plate_n = self._plating[plate_id]
            if test_plate != plate_n:
                common_segment = self.plate_common_segment(test_plate, plate_n)
                if common_segment is not None:
                    adjacent_plate_list.append(plate_n)
        return adjacent_plate_list

    def get_long_adjacent_plates(self, test_plate: Plate):
        # Returns a list of Plate objects which are longitudinally adjacent and
        # have a common transverse primary supporting member with test_plate
        adjacent_plate_list = []
        adjacent_plates = self.get_adjacent_plates(test_plate)
        for i in range(0, len(adjacent_plates)):
            common_segment = self.plate_common_segment(test_plate, adjacent_plates[i])
            if common_segment.primary_supp_mem.direction == BeamDirection.TRANSVERSE:
                adjacent_plate_list.append(adjacent_plates[i])
        return adjacent_plate_list

    def get_tran_adjacent_plates(self, test_plate: Plate):
        # Returns a list of Plate objects which are transversely adjacent and
        # have a common longitudinal primary supporting member with test_plate
        adjacent_plate_list = []
        adjacent_plates = self.get_adjacent_plates(test_plate)
        for i in range(0, len(adjacent_plates)):
            common_segment = self.plate_common_segment(test_plate, adjacent_plates[i])
            if common_segment.primary_supp_mem.direction == BeamDirection.LONGITUDINAL:
                adjacent_plate_list.append(adjacent_plates[i])
        return adjacent_plate_list

    def plating_zones_between_psm(self, test_member1: PrimarySuppMem, test_member2: PrimarySuppMem):
        # Plating zones between adjacent Primary Supporting Members
        identified_zones_list = []
        for plate_id in self._plating.keys():
            plate = self._plating[plate_id]
            if Plate.test_plate_between_psm(plate, test_member1, test_member2):
                identified_zones_list.append(plate)
        return identified_zones_list

    def segments_between_psm(self, test_member1: PrimarySuppMem, test_member2: PrimarySuppMem):
        # Returns all segments between adjacent Primary Supporting Members
        identified_segments_set = set()
        psm_direction = test_member1.direction
        plates = self.plating_zones_between_psm(test_member1, test_member2)

        if psm_direction == BeamDirection.LONGITUDINAL:
            for plate in plates:
                identified_segments_set.add(plate.trans_seg1)
                identified_segments_set.add(plate.trans_seg2)

        elif psm_direction == BeamDirection.TRANSVERSE:
            for plate in plates:
                identified_segments_set.add(plate.long_seg1)
                identified_segments_set.add(plate.long_seg2)

        return identified_segments_set

    def set_longitudinal_PSM_spacing(self, prim_supp_member_id: int, ref_member_id: int, spacing: float):
        """
        :param prim_supp_member_id: ID of the Primary Supporting Member to have relative distance changed based on spacing.
        :param ref_member_id: ID of the Primary Supporting Member in relation to which the spacing value is given.
        :param spacing: Desired Primary Supporting Member spacing in [m].
        :return: Changes the Primary Supporting Member relative distance value based on input spacing value.
                If both prim_supp_member_id and ref_member_id are the same, the Primary Supporting Member is moved from its
                original position by the value of spacing in [m]. Edge members are fixed and can not be moved, their position is
                determined by overall grillgeo dimensions.
        """
        primary_supp_member = self._longitudinal_memb[prim_supp_member_id]
        adjacent_member = self._longitudinal_memb[ref_member_id]
        new_rel_dist = primary_supp_member.rel_dist
        if prim_supp_member_id == 1 or prim_supp_member_id == len(self._longitudinal_memb):
            print("ERROR: Edge Primary Supporting Member positions are fixed! Spacing of this member can not be changed.")
        else:
            if prim_supp_member_id > ref_member_id:
                new_rel_dist = adjacent_member.rel_dist + (spacing / self._B_overall)
            elif prim_supp_member_id < ref_member_id:
                new_rel_dist = adjacent_member.rel_dist - (spacing / self._B_overall)
            else:
                new_rel_dist = primary_supp_member.rel_dist + (spacing / self._B_overall)
        primary_supp_member.rel_dist = new_rel_dist

    def set_transverse_PSM_spacing(self, prim_supp_member_id: int, ref_member_id: int, spacing: float):
        """
        :param prim_supp_member_id: ID of the Primary Supporting Member to have relative distance changed based on spacing.
        :param ref_member_id: ID of the Primary Supporting Member in relation to which the spacing value is given.
        :param spacing: Desired Primary Supporting Member spacing in [m].
        :return: Changes the Primary Supporting Member relative distance value based on input spacing value.
                If both prim_supp_member_id and ref_member_id are the same, the Primary Supporting Member is moved from its
                original position by the value of spacing in [m]. Edge members are fixed and can not be moved, their position is
                determined by overall grillgeo dimensions.
        """
        primary_supp_member = self._transverse_memb[prim_supp_member_id]
        adjacent_member = self._transverse_memb[ref_member_id]
        new_rel_dist = primary_supp_member.rel_dist
        if prim_supp_member_id == 1 or prim_supp_member_id == len(self._transverse_memb):
            print("ERROR: Edge Primary Supporting Member positions are fixed! Spacing of this member can not be changed.")
        else:
            if prim_supp_member_id > ref_member_id:
                new_rel_dist = adjacent_member.rel_dist + (spacing / self._L_overall)
            elif prim_supp_member_id < ref_member_id:
                new_rel_dist = adjacent_member.rel_dist - (spacing / self._L_overall)
            else:
                new_rel_dist = primary_supp_member.rel_dist + (spacing / self._L_overall)
        primary_supp_member.rel_dist = new_rel_dist

    def set_all_longitudinal_PSM(self, *args):
        """
        :param args: Spacing values between all longitudinal Primary Supporting Members in [m], starting from csy origin.
        :return: Sets longitudinal Primary Supporting Member relative distance values. Does not change the position of edge PSM.
        """
        n_inputs = len(args)
        spacing_n = 0

        if n_inputs == self.N_longitudinal - 2:
            for psm_id in range(2, self.N_longitudinal):
                Grillage.set_longitudinal_PSM_spacing(self, psm_id, psm_id - 1, args[spacing_n])
                spacing_n += 1

        elif n_inputs == self.N_longitudinal - 1:
            sumcheck = 0
            for i in range(0, n_inputs):
                sumcheck += args[i]

            if sumcheck != self.B_overall:
                print("ERROR: The sum of entered spacing values of", sumcheck, "m does not match overall grillgeo width of", self.B_overall,
                      "m! Adjust spacing values to match or input one less spacing value.")
            else:
                for psm_id in range(2, self.N_longitudinal):
                    Grillage.set_longitudinal_PSM_spacing(self, psm_id, psm_id - 1, args[spacing_n])
                    spacing_n += 1

        elif n_inputs >= self.N_longitudinal:
            print("ERROR: Number of input spacing values for setting all longitudinal primary supporting member positions:", n_inputs,
                  ", expected at most", self.N_longitudinal - 1, "spacing values for grillgeo with",
                  self.N_longitudinal, "longitudinal primary supporting members.")

        else:
            print("ERROR: Number of input spacing values for setting all longitudinal primary supporting member positions:", n_inputs,
                  ", expected at least", self.N_longitudinal - 2, "spacing values for grillgeo with",
                  self.N_longitudinal, "longitudinal primary supporting members.")

    def set_all_transverse_PSM(self, *args):
        """
        :param args: Spacing values between all transverse Primary Supporting Members in [m], starting from csy origin.
        :return: Sets transverse Primary Supporting Member relative distance values. Does not change the position of edge PSM.
        """
        n_inputs = len(args)
        spacing_n = 0

        if n_inputs == self.N_transverse - 2:
            for psm_id in range(2, self.N_transverse):
                Grillage.set_transverse_PSM_spacing(self, psm_id, psm_id - 1, args[spacing_n])
                spacing_n += 1

        elif n_inputs == self.N_transverse - 1:
            sumcheck = 0
            for i in range(0, n_inputs):
                sumcheck += args[i]

            if sumcheck != self.L_overall:
                print("ERROR: The sum of entered spacing values of", sumcheck, "m does not match overall grillgeo length of", self.L_overall,
                      "m! Adjust spacing values to match or input one less spacing value.")
            else:
                for psm_id in range(2, self.N_transverse):
                    Grillage.set_transverse_PSM_spacing(self, psm_id, psm_id - 1, args[spacing_n])
                    spacing_n += 1

        elif n_inputs >= self._N_transverse:
            print("ERROR: Number of input spacing values for setting all transverse primary supporting member positions:", n_inputs,
                  ", expected at most", self.N_transverse - 1, "spacing values for grillgeo with",
                  self.N_transverse, "transverse primary supporting members.")

        else:
            print("ERROR: Number of input spacing values for setting all transverse primary supporting member positions:", n_inputs,
                  ", expected at least", self.N_transverse - 2, "spacing values for grillgeo with",
                  self.N_transverse, "transverse primary supporting members.")

    def set_plating_prop_longitudinals(self, plate_id, plating_property, property_object):
        """
        :param plate_id: ID of the Plating zone to be changed, along with all other plating zones between adjacent longitudinal
                primary supporting members. Defines longitudinal primary supporting members by its longitudinal segments.
        :param plating_property: Plating zone property to be changed: Plate property (plate_prop), Stiffener layout (stiff_layout),
                or Stiffener direction (stiff_dir).
        :param property_object: PlateProperty, StiffenerLayout or BeamOrientation
                (either BeamOrientation.LONGITUDINAL or BeamOrientation.TRANSVERSE).
        :return: Changes all plating zone plating_property between longitudinal primary supporting members to property_value.
        """
        plate = self.plating()[plate_id]
        long_1 = plate.long_seg1.primary_supp_mem
        long_2 = plate.long_seg2.primary_supp_mem

        plating_zones = self.plating_zones_between_psm(long_1, long_2)

        if plating_property == "plate_prop":
            for i in range(0, len(plating_zones)):
                plate_id = plating_zones[i].id
                self.plating()[plate_id].plate_prop = property_object

        elif plating_property == "stiff_layout":
            for i in range(0, len(plating_zones)):
                plate_id = plating_zones[i].id
                self.plating()[plate_id].stiff_layout = property_object

        elif plating_property == "stiff_dir":
            for i in range(0, len(plating_zones)):
                plate_id = plating_zones[i].id
                self.plating()[plate_id].stiff_dir = property_object

        else:
            print("ERROR! Unknown plating property. Choose plate_prop, stiff_layout or stiff_dir")

    def set_plating_prop_transversals(self, plate_id, plating_property, property_object):
        """
        :param plate_id: ID of the Plating zone to be changed, along with all other plating zones between adjacent transverse
                primary supporting members. Defines transverse primary supporting members by its transverse segments.
        :param plating_property: Plating zone property to be changed: Plate property (plate_prop), Stiffener layout (stiff_layout),
                or Stiffener direction (stiff_dir).
        :param property_object: PlateProperty, StiffenerLayout or BeamOrientation
                (either BeamOrientation.LONGITUDINAL or BeamOrientation.TRANSVERSE).
        :return: Changes all plating zone plating_property between transverse primary supporting members to property_value.
        """
        plate = self.plating()[plate_id]
        tran_1 = plate.trans_seg1.primary_supp_mem
        tran_2 = plate.trans_seg2.primary_supp_mem

        plating_zones = self.plating_zones_between_psm(tran_1, tran_2)

        if plating_property == "plate_prop":
            for i in range(0, len(plating_zones)):
                plate_id = plating_zones[i].id
                self.plating()[plate_id].plate_prop = property_object

        elif plating_property == "stiff_layout":
            for i in range(0, len(plating_zones)):
                plate_id = plating_zones[i].id
                self.plating()[plate_id].stiff_layout = property_object

        elif plating_property == "stiff_dir":
            for i in range(0, len(plating_zones)):
                plate_id = plating_zones[i].id
                self.plating()[plate_id].stiff_dir = property_object

        else:
            print("ERROR! Unknown plating property. Choose plate_prop, stiff_layout or stiff_dir")

    def set_plating_prop_symmetric(self, plate_id, plating_property, property_object):
        """
        :param plate_id: ID of the Plating zone to be changed, along with all other symmetric plating zones.
        :param plating_property: Plating zone property to be changed: Plate property (plate_prop), Stiffener layout (stiff_layout),
                or Stiffener direction (stiff_dir).
        :param property_object: PlateProperty, StiffenerLayout or BeamOrientation
                (either BeamOrientation.LONGITUDINAL or BeamOrientation.TRANSVERSE).
        :return: Changes all symmetric plating zone plating_property to property_value.
        """
        plate = self.plating()[plate_id]
        symmetric_plate_list = plate.symmetric_plate_zones

        if plating_property == "plate_prop":
            plate.plate_prop = property_object
            for plate in symmetric_plate_list:
                plate.plate_prop = property_object

        elif plating_property == "stiff_layout":
            plate.stiff_layout = property_object
            for plate in symmetric_plate_list:
                plate.stiff_layout = property_object

        elif plating_property == "stiff_dir":
            plate.stiff_dir = property_object
            for plate in symmetric_plate_list:
                plate.stiff_dir = property_object

        else:
            print("ERROR! Unknown plating property. Choose plate_prop, stiff_layout or stiff_dir")

    def assign_symmetric_segments(self):
        for i_long in self._longitudinal_memb.values():
            for i_segment in i_long.segments:
                symmetric_id = len(i_long.segments) - i_segment.id + 1
                i_segment.symmetric_segment = i_long.segments[symmetric_id - 1]

        for i_tran in self._transverse_memb.values():
            for i_segment in i_tran.segments:
                symmetric_id = len(i_tran.segments) - i_segment.id + 1
                i_segment.symmetric_segment = i_tran.segments[symmetric_id - 1]

    def set_long_symm_segment_beam_property(self, prim_supp_member_id: int, segment_id: int, beam_property: BeamProperty):
        """
        :param prim_supp_member_id: ID of a Primary Supporting Member the Segment belongs to.
        :param segment_id: ID of the Segment to have beam property changed. Segments of all Primary Supporting Members are
                enumerated starting from 1, in the direction of the global x (longitudinal) or y (transverse) axis.
        :param beam_property: Beam property to be assigned to the selected segment and its symmetric segmenets.
        :return: Changes all symmetric Segment beam properties belonging to the selected (prim_supp_member_id) and symmetric
                Primary Supporting Member.
        """
        segment = self._longitudinal_memb[prim_supp_member_id].segments[segment_id - 1]
        symmetric_member = self._longitudinal_memb[prim_supp_member_id].symmetric_member

        if symmetric_member is None:
            symmetric_segment_1 = segment.symmetric_segment
            segment.beam_prop = beam_property
            symmetric_segment_1.beam_prop = beam_property

        else:
            symmetric_segment_1 = segment.symmetric_segment
            symmetric_segment_2 = symmetric_member.segments[segment_id - 1]
            symmetric_segment_3 = symmetric_segment_2.symmetric_segment

            segment.beam_prop = beam_property
            symmetric_segment_1.beam_prop = beam_property
            symmetric_segment_2.beam_prop = beam_property
            symmetric_segment_3.beam_prop = beam_property

        """
        segment = self._longitudinal_memb[prim_supp_member_id].segments[segment_id - 1]
        symmetric_member = self._longitudinal_memb[prim_supp_member_id].symmetric_member
        symmetric_segment_1 = segment.symmetric_segment
        symmetric_segment_2 = symmetric_member.segments[segment_id - 1]
        symmetric_segment_3 = symmetric_segment_2.symmetric_segment

        segment.beam_prop = beam_property
        symmetric_segment_1.beam_prop = beam_property
        symmetric_segment_2.beam_prop = beam_property
        symmetric_segment_3.beam_prop = beam_property
        """

    def set_long_member_beam_property(self, prim_supp_member_id: int, beam_property: BeamProperty):
        """
        :param prim_supp_member_id: ID of the Primary Supporting Member to have all Segment beam properties changed.
        :param beam_property: New BeamProperty for the entire Primary Supporting Member
        :return: Sets all Segment BeamProperty belonging to the selected Primary Supporting Member.
        """
        prim_supp_mem = self.longitudinal_members()[prim_supp_member_id]

        for segment in prim_supp_mem.segments:
            segment.beam_prop = beam_property

    def set_tran_member_beam_property(self, prim_supp_member_id: int, beam_property: BeamProperty):
        """
        :param prim_supp_member_id: ID of the Primary Supporting Member to have all Segment beam properties changed.
        :param beam_property: New BeamProperty for the entire Primary Supporting Member
        :return: Sets all Segment BeamProperty belonging to the selected Primary Supporting Member.
        """
        prim_supp_mem = self.transverse_members()[prim_supp_member_id]

        for segment in prim_supp_mem.segments:
            segment.beam_prop = beam_property

    def hc_variant_check(self):
        hc_check = True

        # 1.)   Plating zones between two common Primary Supporting Members may not have the same stiffener orientation
        #       and different stiffener spacing
        # Between longitudinal primary supporting members:
        for psm_id in range(1, self._N_longitudinal):
            psm_1 = self._longitudinal_memb[psm_id]
            psm_2 = self._longitudinal_memb[psm_id + 1]
            plate_list = self.plating_zones_between_psm(psm_1, psm_2)
            plate_combinations = itertools.combinations(plate_list, 2)
            # Plate combinations compares each plating zone stored inside plate_list with others, without duplicate checks
            for plate in list(plate_combinations):
                plate1 = plate[0]   # First plating zone of plate combinations
                plate2 = plate[1]   # Second plating zone of plate combinations
                spacing1 = Plate.get_stiffener_spacing(plate1)
                spacing2 = Plate.get_stiffener_spacing(plate2)
                if plate1.stiff_dir == BeamDirection.LONGITUDINAL and plate2.stiff_dir == BeamDirection.LONGITUDINAL and spacing1 != spacing2:
                    print("Stiffener spacing of longitudinal stiffeners on plating zones between adjacent longitudinal primary supporting members",
                          psm_1.id, "and", psm_2.id, "do not match! \n", "     Plating zone", plate1.id, "has stiffener spacing of",
                          spacing1, "m", ", while plating zone", plate2.id, "has spacing of", spacing2, "m")
                    hc_check = False

        # Between transverse primary supporting members:
        for psm_id in range(1, self._N_transverse):
            psm_1 = self._transverse_memb[psm_id]
            psm_2 = self._transverse_memb[psm_id + 1]
            plate_list = self.plating_zones_between_psm(psm_1, psm_2)
            plate_combinations = itertools.combinations(plate_list, 2)
            # Plate combinations compares each plating zone stored inside plate_list with others, without duplicate checks
            for plate in list(plate_combinations):
                plate1 = plate[0]   # First plating zone of plate combinations
                plate2 = plate[1]   # Second plating zone of plate combinations
                spacing1 = Plate.get_stiffener_spacing(plate1)
                spacing2 = Plate.get_stiffener_spacing(plate2)
                if plate1.stiff_dir == BeamDirection.TRANSVERSE and plate2.stiff_dir == BeamDirection.TRANSVERSE and spacing1 != spacing2:
                    print("Stiffener spacing of transverse stiffeners on plating zones between adjacent transverse primary supporting members",
                          psm_1.id, "and", psm_2.id, "do not match! \n", "     Plating zone", plate1.id, "has stiffener spacing of",
                          spacing1, "m", ", while plating zone", plate2.id, "has spacing of", spacing2, "m")
                    hc_check = False
        return hc_check


class GrillageModelData:
    def __init__(self, filename: str):
        self._filename = filename

    # Save grillgeo model to a file
    def write_file(self, grillage: Grillage):
        n_longitudinal = str(len(Grillage.longitudinal_members(grillage).keys()))
        n_transverse = str(len(Grillage.transverse_members(grillage).keys()))
        grillage_length = str(grillage.L_overall)
        grillage_width = str(grillage.B_overall)

        # Number of materials, corrosion additions, beam and plate properties, stiffener layouts
        n_material_prop = str(len(Grillage.material_props(grillage).keys()))
        n_corrosion_add = str(len(Grillage.corrosion_addition(grillage).keys()))
        n_beam_prop = str(len(Grillage.beam_props(grillage).keys()))
        n_plate_prop = str(len(Grillage.plate_props(grillage).keys()))
        n_stiff_layout = str(len(Grillage.stiffener_layouts(grillage).keys()))

        # Number of longitudinal and transverse segments
        n_long_segments = str((len(Grillage.transverse_members(grillage).keys()) - 1) *
                              len(Grillage.longitudinal_members(grillage).keys()))
        n_tran_segments = str((len(Grillage.longitudinal_members(grillage).keys()) - 1) *
                              len(Grillage.transverse_members(grillage).keys()))

        with open(self._filename, "w") as f:
            # Write input for grillgeo object - dimensions and number of primary supporting members
            f.write(grillage_length + "," + grillage_width + "," + n_longitudinal + "," + n_transverse + "\n")  # Grillage object

            # Write number of saved properties
            f.write(n_material_prop + "," + n_corrosion_add + "," + n_beam_prop + "," + n_plate_prop + "," + n_stiff_layout + ","
                    + n_longitudinal + "," + n_transverse + "," + n_long_segments + "," + n_tran_segments + "\n")

            # Write material properties
            for i in grillage.material_props():
                mat_property = grillage.material_props()[i]
                line = str(mat_property.id)
                line += ',' + str(mat_property.E)
                line += ',' + str(mat_property.v)
                line += ',' + str(mat_property.ro)
                line += ',' + str(mat_property.Reh)
                line += ',' + str(mat_property.name) + "\n"
                f.write(line)

            # Write corrosion additions
            for i in grillage.corrosion_addition():
                corr_add = grillage.corrosion_addition()[i]
                line = str(corr_add.id)
                line += ',' + str(corr_add.tc) + "\n"
                f.write(line)

            # Write beam properties
            for i in grillage.beam_props():
                beam_property = grillage.beam_props()[i]

                if isinstance(beam_property, TBeamProperty) and not\
                        (isinstance(beam_property, FBBeamProperty) or isinstance(beam_property, LBeamProperty)):
                    line = str(beam_property.id)
                    line += ",T"
                    line += ',' + str(beam_property.hw)
                    line += ',' + str(beam_property.tw)
                    line += ',' + str(beam_property.bf)
                    line += ',' + str(beam_property.tf)
                    line += ',' + str(beam_property.mat.id) + "\n"
                    f.write(line)

                elif isinstance(beam_property, FBBeamProperty):
                    line = str(beam_property.id)
                    line += ",FB"
                    line += ',' + str(beam_property.hw)
                    line += ',' + str(beam_property.tw)
                    line += ',' + str(beam_property.bf)
                    line += ',' + str(beam_property.tf)
                    line += ',' + str(beam_property.mat.id) + "\n"
                    f.write(line)

                elif isinstance(beam_property, LBeamProperty):
                    line = str(beam_property.id)
                    line += ",L"
                    line += ',' + str(beam_property.hw)
                    line += ',' + str(beam_property.tw)
                    line += ',' + str(beam_property.bf)
                    line += ',' + str(beam_property.tf)
                    line += ',' + str(beam_property.mat.id) + "\n"
                    f.write(line)

                elif isinstance(beam_property, BulbBeamProperty):
                    line = str(beam_property.id)
                    line += ",HP"
                    line += ',' + str(beam_property.hw_HP)
                    line += ',' + str(beam_property.tw_HP)
                    line += ',' + str(beam_property.mat.id) + "\n"
                    f.write(line)

                elif isinstance(beam_property, HatBeamProperty):
                    line = str(beam_property.id)
                    line += ",HAT"
                    line += ',' + str(beam_property.h)
                    line += ',' + str(beam_property.t)
                    line += ',' + str(beam_property.bf)
                    line += ',' + str(beam_property.fi)
                    line += ',' + str(beam_property.mat.id) + "\n"
                    f.write(line)

            # Write plate properties
            for i in grillage.plate_props():
                plate_property = grillage.plate_props()[i]
                line = str(plate_property.id)
                line += ',' + str(plate_property.tp)
                line += ',' + str(plate_property.plate_mat.id) + "\n"
                f.write(line)

            # Write stiffener layouts
            for i in grillage.stiffener_layouts():
                stiff_layout = grillage.stiffener_layouts()[i]
                line = str(stiff_layout.id)
                line += ',' + str(stiff_layout.beam_prop.id)
                line += ',' + str(stiff_layout.definition_type)
                line += ',' + str(stiff_layout.definition_value) + "\n"
                f.write(line)

            # Write longitudinal primary supporting members
            for i in grillage.longitudinal_members():
                long_member = grillage.longitudinal_members()[i]
                line = str(long_member.id)
                line += ',' + str(long_member.direction).split(".")[1]
                line += ',' + str(long_member.rel_dist) + "\n"
                f.write(line)

            # Write transverse primary supporting members
            for i in grillage.transverse_members():
                tran_member = grillage.transverse_members()[i]
                line = str(tran_member.id)
                line += ',' + str(tran_member.direction).split(".")[1]
                line += ',' + str(tran_member.rel_dist) + "\n"
                f.write(line)

            # Write longitudinal segments
            for i in grillage.longitudinal_members():
                for segment in range(0, int(n_transverse) - 1):
                    curr_segment = grillage.longitudinal_members()[i].segments[segment]
                    line = str(curr_segment.id)
                    line += ',' + str(curr_segment.beam_prop.id)
                    line += ',' + str(curr_segment.primary_supp_mem.id)
                    line += ',' + str(curr_segment.cross_member1.id)
                    line += ',' + str(curr_segment.cross_member2.id) + "\n"
                    f.write(line)

            # Write transverse segments
            for i in grillage.transverse_members():
                for segment in range(0, int(n_longitudinal) - 1):
                    curr_segment = grillage.transverse_members()[i].segments[segment]
                    line = str(curr_segment.id)
                    line += ',' + str(curr_segment.beam_prop.id)
                    line += ',' + str(curr_segment.primary_supp_mem.id)
                    line += ',' + str(curr_segment.cross_member1.id)
                    line += ',' + str(curr_segment.cross_member2.id) + "\n"
                    f.write(line)

            # Write plating zones
            for i in grillage.plating().keys():
                plate = grillage.plating()[i]
                line = str(plate.id)
                line += ',' + str(plate.plate_prop.id)
                line += ',' + str(plate.long_seg1.primary_supp_mem.id)
                line += ',' + str(plate.long_seg1.id)
                line += ',' + str(plate.trans_seg1.primary_supp_mem.id)
                line += ',' + str(plate.trans_seg1.id)
                line += ',' + str(plate.long_seg2.primary_supp_mem.id)
                line += ',' + str(plate.long_seg2.id)
                line += ',' + str(plate.trans_seg2.primary_supp_mem.id)
                line += ',' + str(plate.trans_seg2.id)
                line += ',' + str(plate.stiff_layout.id)
                line += ',' + str(plate.stiff_dir).split(".")[1]
                line += ',' + str(plate.ref_edge).split(".")[1] + "\n"
                f.write(line)
            f.close()

    # Generate model from a saved file
    def read_file(self):
        with open(self._filename, "r") as f:
            # Generate Grillage object - read from the first line
            lines = f.readlines()
            iline = 0
            line = lines[iline]
            split = line.split(",")
            grillage_variant = Grillage(float(split[0]), float(split[1]), int(split[2]), int(split[3]))
            n_plate_zones = (int(split[2]) - 1) * (int(split[3]) - 1)

            # Number of all saved properties - read from the second line
            iline = 1
            line = lines[iline]
            split = line.split(",")
            n_material_prop = int(split[0])
            n_corrosion_add = int(split[1])
            n_beam_prop = int(split[2])
            n_plate_prop = int(split[3])
            n_stiff_layout = int(split[4])
            n_longitudinal = int(split[5])
            n_transverse = int(split[6])
            n_long_segments = int(split[7])
            n_tran_segments = int(split[8])

            # Generate MaterialProperty - read from the third line
            line_start = 2
            line_end = line_start + n_material_prop
            for i in range(line_start, line_end):
                line = lines[i]
                split = line.split(",")
                material_prop = MaterialProperty(int(split[0]), float(split[1]), float(split[2]), float(split[3]), float(split[4]), str(split[5].strip()))
                Grillage.add_material(grillage_variant, material_prop)

            # Generate CorrosionAddition
            line_start = line_end
            line_end = line_start + n_corrosion_add
            for i in range(line_start, line_end):
                line = lines[i]
                split = line.split(",")
                corr_add = CorrosionAddition(int(split[0]), float(split[1]))
                Grillage.add_corrosion_addition(grillage_variant, corr_add)

            # Generate BeamProperty
            line_start = line_end
            line_end = line_start + n_beam_prop
            for i in range(line_start, line_end):
                line = lines[i]
                split = line.split(",")
                if split[1] == "FB":
                    material = Grillage.material_props(grillage_variant)[int(split[6])]
                    beam_prop = FBBeamProperty(int(split[0]), float(split[2]), float(split[3]), material, float(split[4]), float(split[5]),)
                    Grillage.add_beam_prop(grillage_variant, beam_prop)
                if split[1] == "L":
                    material = Grillage.material_props(grillage_variant)[int(split[6])]
                    beam_prop = LBeamProperty(int(split[0]), float(split[2]), float(split[3]), float(split[4]), float(split[5]), material)
                    Grillage.add_beam_prop(grillage_variant, beam_prop)
                if split[1] == "T":
                    material = Grillage.material_props(grillage_variant)[int(split[6])]
                    beam_prop = TBeamProperty(int(split[0]), float(split[2]), float(split[3]), float(split[4]), float(split[5]), material)
                    Grillage.add_beam_prop(grillage_variant, beam_prop)
                if split[1] == "HP":
                    material = Grillage.material_props(grillage_variant)[int(split[4])]
                    beam_prop = BulbBeamProperty(int(split[0]), float(split[2]), float(split[3]), material)
                    Grillage.add_beam_prop(grillage_variant, beam_prop)
                if split[1] == "HAT":
                    material = Grillage.material_props(grillage_variant)[int(split[6])]
                    beam_prop = HatBeamProperty(int(split[0]), float(split[2]), float(split[3]), float(split[4]), float(split[5]), material)
                    Grillage.add_beam_prop(grillage_variant, beam_prop)

            # Generate PlateProptery
            line_start = line_end
            line_end = line_start + n_plate_prop
            for i in range(line_start, line_end):
                line = lines[i]
                split = line.split(",")
                material = Grillage.material_props(grillage_variant)[int(split[2])]
                plate_prop = PlateProperty(int(split[0]), float(split[1]), material)
                Grillage.add_plate_prop(grillage_variant, plate_prop)

            # Generate StiffenerLayout
            line_start = line_end
            line_end = line_start + n_stiff_layout
            for i in range(line_start, line_end):
                line = lines[i]
                split = line.split(",")
                beam_prop = Grillage.beam_props(grillage_variant)[int(split[1])]
                stiff_layout = StiffenerLayout(int(split[0]), beam_prop, split[2], float(split[3]))
                Grillage.add_stiffener_layout(grillage_variant, stiff_layout)

            # Generate longitudinal Primary supporting members
            line_start = line_end
            line_end = line_start + n_longitudinal
            i_member = 1
            for i in range(line_start, line_end):
                line = lines[i]
                split = line.split(",")
                current_member = PrimarySuppMem(int(split[0]), BeamDirection[split[1].strip()], float(split[2]), grillage_variant)
                Grillage.longitudinal_members(grillage_variant)[i_member] = current_member
                i_member += 1

            # Generate transverse Primary supporting members
            line_start = line_end
            line_end = line_start + n_transverse
            i_member = 1
            for i in range(line_start, line_end):
                line = lines[i]
                split = line.split(",")
                current_member = PrimarySuppMem(int(split[0]), BeamDirection[split[1].strip()], float(split[2]), grillage_variant)
                Grillage.transverse_members(grillage_variant)[i_member] = current_member
                i_member += 1

            # Generate longitudinal Segments
            line_start = line_end
            line_end = line_start + n_long_segments
            for i in range(line_start, line_end):
                line = lines[i]
                split = line.split(",")
                beam_prop = Grillage.beam_props(grillage_variant)[int(split[1])]
                primary_supp_mem = Grillage.longitudinal_members(grillage_variant)[int(split[2])]
                cross_member1 = Grillage.transverse_members(grillage_variant)[int(split[3])]
                cross_member2 = Grillage.transverse_members(grillage_variant)[int(split[4])]
                current_segment = Segment(int(split[0]), beam_prop, primary_supp_mem, cross_member1, cross_member2)
                Grillage.longitudinal_members(grillage_variant)[int(split[2])].segments.append(current_segment)

            # Generate transverse Segments
            line_start = line_end
            line_end = line_start + n_tran_segments
            for i in range(line_start, line_end):
                line = lines[i]
                split = line.split(",")
                beam_prop = Grillage.beam_props(grillage_variant)[int(split[1])]
                primary_supp_mem = Grillage.transverse_members(grillage_variant)[int(split[2])]
                cross_member1 = Grillage.longitudinal_members(grillage_variant)[int(split[3])]
                cross_member2 = Grillage.longitudinal_members(grillage_variant)[int(split[4])]
                current_segment = Segment(int(split[0]), beam_prop, primary_supp_mem, cross_member1, cross_member2)
                Grillage.transverse_members(grillage_variant)[int(split[2])].segments.append(current_segment)

            # Generate Plating
            line_start = line_end
            line_end = line_start + n_plate_zones
            i_plate = 1
            for i in range(line_start, line_end):
                line = lines[i]
                split = line.split(",")
                plate_prop = Grillage.plate_props(grillage_variant)[int(split[1])]
                long_seg1 = Grillage.longitudinal_members(grillage_variant)[int(split[2])].segments[int(split[3]) - 1]
                trans_seg1 = Grillage.transverse_members(grillage_variant)[int(split[4])].segments[int(split[5]) - 1]
                long_seg2 = Grillage.longitudinal_members(grillage_variant)[int(split[6])].segments[int(split[7]) - 1]
                trans_seg2 = Grillage.transverse_members(grillage_variant)[int(split[8])].segments[int(split[9]) - 1]
                stifflayout = Grillage.stiffener_layouts(grillage_variant)[int(split[10])]
                stiff_dir = BeamDirection[split[11].strip()]
                ref_edge = Ref[split[12].strip()]
                curr_plate = Plate(int(split[0]), plate_prop, long_seg1, trans_seg1, long_seg2, trans_seg2, stifflayout, stiff_dir, ref_edge)
                Grillage.add_plating(grillage_variant, curr_plate)
                i_plate += 1
            return grillage_variant
