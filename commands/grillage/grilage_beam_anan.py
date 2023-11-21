from grillage.grillage_model import Grillage,BeamProperty,MaterialProperty,PrimarySuppMem,Plate, PlateProperty
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

    def getE(self):
        return self._prop.mat.E

    @property
    def prop(self)->BeamProperty:
        return self._prop


def generate_grillage_analysis(grillage:Grillage)->anan.GrillageAnalysis:
    grillan = anan.GrillageAnalysis()
    n_long = len(Grillage.longitudinal_members(grillage).keys())
    n_tran = len(Grillage.transverse_members(grillage).keys())
    length = grillage.L_overall
    width = grillage.B_overall

    # in current variant main beams are always longitudinal
    corr_add=grillage.corrosion_addition()[1]
    idb= 0
    #main beams
    beam_plate_data = {}
    for psm1 in grillage.longitudinal_members().values():
        np=0
        tp=0.0
        bp=0.0
        for psm2 in grillage.longitudinal_members().values():
            plates:List[Plate] = grillage.plating_zones_between_psm(psm1, psm2)
            if len(plates) > 0:
                propp:PlateProperty = plates[0].plate_prop
                tp+= propp.tp
                np+=1
                bp+= abs(psm1.rel_dist-psm2.rel_dist)/2.0*width
                pass
        beam_plate_data[psm1.id:[bp, tp/np]]

    for psm1 in grillage.transverse_members().values():
        np=0
        tp=0.0
        bp=0.0
        for psm2 in grillage.transverse_members().values():
            plates:List[Plate] = grillage.plating_zones_between_psm(psm1, psm2)
            if len(plates) > 0:
                propp:PlateProperty = plates[0].plate_prop
                tp+= propp.tp
                np+=1
                bp+= abs(psm1.rel_dist-psm2.rel_dist)/2.0*length
                pass


        beam_plate_data[psm1.id:[bp, tp/np]]

    for psm in grillage.longitudinal_members().values():
        segment = psm.segments[0]
        propm = segment.beam_prop
        n1,n2 = psm.end_nodes # Grillage model nodes have the same interface as Grillage Analysis nodes
        grillan.add_node(n1)
        grillan.add_node(n2)

        propa = BeamPropertAnAn(propm,beam_plate_data[psm.id][0],beam_plate_data[psm.id][1],corr_add)
        grillan.add_prop(propa)
        idb+=1
        beam = anan.AnalyticBeam(idb,propa,n1,n2)
        grillan.add_beam(beam)
    #cross beams
    for psm in grillage.transverse_members().values():
        segment = psm.segments[0]
        propm = segment.beam_prop
        n1,n2 = psm.end_nodes # Grillage model nodes have the same interface as Grillage Analysis nodes
        grillan.add_node(n1)
        grillan.add_node(n2)
        propa = BeamPropertAnAn(propm,0.0,0.0,0.0)
        grillan.add_prop(propa)
        idb+=1
        beam = anan.AnalyticBeam(idb,propa,n1,n2)
        grillan.add_beam(beam)
    return grillan

