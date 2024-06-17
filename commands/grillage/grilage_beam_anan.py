import numpy as np

from grillage.grillage_model import Grillage,BeamProperty,MaterialProperty,PrimarySuppMem,Plate, PlateProperty,Segment
import anandir.beam_analysis as anan
from typing import List,Dict
import math

class BeamPropertAnAn(anan.BeamProperty):
    def __init__(self,prop:BeamProperty,bp,tp,cor_add):
        super(BeamPropertAnAn, self).__init__(prop.id)
        self._prop:BeamProperty = prop
        self._bp = bp
        self._tp = tp
        self._cor_add=cor_add


    def getW(self):
        Wmin = self._prop.getWmin(self._bp,self._tp,self._cor_add)
        return Wmin
    def getIy(self):
        Iy = self._prop.get_Iy_I(self._bp,self._tp,self._cor_add)
        return Iy

    def getE(self):
        return self._prop.mat.E

    @property
    def prop(self)->BeamProperty:
        return self._prop

class NodeAnAn(anan.Node):
    def __init__(self, id: int, cords:np.ndarray):
        super(NodeAnAn, self).__init__(id,cords = cords)



def generate_grillage_analysis(grillage:Grillage)->anan.GrillageAnalysis:
    grillan = anan.GrillageAnalysis()
    n_long = len(Grillage.longitudinal_members(grillage).keys())
    n_tran = len(Grillage.transverse_members(grillage).keys())
    length = grillage.L_overall
    width = grillage.B_overall

    # in current variant main beams are always longitudinal
    corr_add=grillage.corrosion_addition()[1]
    idb = 0
    idl = 0
    idn = 0
    #main beams
    long_beam_plate_data = {}
    for psm1 in grillage.longitudinal_members().values():
        for psm2 in grillage.longitudinal_members().values():
            plates:List[Plate] = grillage.plating_zones_between_psm(psm1, psm2)
            if len(plates) > 0:
                propp:PlateProperty = plates[0].plate_prop
                psm1_data = long_beam_plate_data.get(psm1.id)
                psm2_data = long_beam_plate_data.get(psm2.id)
                bp = abs(psm1.rel_dist - psm2.rel_dist) / 2.0 * width
                tp = propp.tp
                if psm1_data is None:
                    long_beam_plate_data[psm1.id] = [bp, tp]
                else:
                    [bp_, tp_] = psm1_data
                    long_beam_plate_data[psm1.id] =  [bp+bp_, (tp+tp_)/2.0]
                if psm2_data is None:
                    long_beam_plate_data[psm2.id] = [bp, tp]
                else:
                    [bp_, tp_] = psm2_data
                    long_beam_plate_data[psm1.id] =  [bp+bp_, (tp+tp_)/2.0]

    trans_beam_plate_data = {}
    for psm1 in grillage.transverse_members().values():
        for psm2 in grillage.transverse_members().values():
            plates:List[Plate] = grillage.plating_zones_between_psm(psm1, psm2)
            if len(plates) > 0:
                psm1_data = trans_beam_plate_data.get(psm1.id)
                psm2_data = trans_beam_plate_data.get(psm2.id)
                bp = abs(psm1.rel_dist - psm2.rel_dist) / 2.0 * width
                tp = propp.tp
                if psm1_data is None:
                    trans_beam_plate_data[psm1.id] = [bp, tp]
                else:
                    [bp_, tp_] = psm1_data
                    trans_beam_plate_data[psm1.id] = [bp + bp_, (tp + tp_) / 2.0]
                if psm2_data is None:
                    trans_beam_plate_data[psm2.id] = [bp, tp]
                else:
                    [bp_, tp_] = psm2_data
                    trans_beam_plate_data[psm1.id] = [bp + bp_, (tp + tp_) / 2.0]

    for psm in grillage.longitudinal_members().values():
        segment = psm.segments[0]
        propm = segment.beam_prop
        cords1,cords2 = psm.end_nodes_coords # Grillage model nodes have the same interface as Grillage Analysis nodes
        cords1 = cords1 / 1000.0  # from m to mm
        cords2 = cords2 / 1000.0  # from m to mm
        idn+=1
        n1= NodeAnAn(idn,cords1)
        grillan.add_node(n1)
        idn += 1
        n2 = NodeAnAn(idn, cords2)
        grillan.add_node(n2)
        propa = BeamPropertAnAn(propm,long_beam_plate_data[psm.id][0],long_beam_plate_data[psm.id][1],corr_add)
        grillan.add_prop(propa)
        idb+=1
        beam = anan.AnalyticBeam(idb,propa,n1,n2)
        grillan.add_beam(beam)
        idl+=1
        load = anan.BeamLoadNoLoad(idl)
        beam.add_load(load)

    #cross beams
    for psm in grillage.transverse_members().values():
        segment = psm.segments[0]
        propm = segment.beam_prop
        cords1, cords2 = psm.end_nodes_coords  # Grillage model nodes have the same interface as Grillage Analysis nodes
        cords1 = cords1 / 1000.0 # from m to mm
        cords2 = cords2 / 1000.0 # from m to mm
        idn += 1
        n1 = NodeAnAn(idn, cords1)
        grillan.add_node(n1)
        idn += 1
        n2 = NodeAnAn(idn, cords2)
        grillan.add_node(n2)
        propa = BeamPropertAnAn(propm,trans_beam_plate_data[psm.id][0],trans_beam_plate_data[psm.id][1],corr_add)
        grillan.add_prop(propa)
        idb+=1
        beam = anan.AnalyticBeam(idb,propa,n1,n2)
        grillan.add_beam(beam)
        idl += 1
        q0 = grillage.pressure*trans_beam_plate_data[psm.id][0]
        load = anan.BeamLoadConstContLoad(idl, q0)
        beam.add_load(load)
    return grillan

