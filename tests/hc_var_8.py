"""
TESTNI MODUL

Generacija nove topologije i spremanje iste u datoteku
"""

from grillage.grillage_model import *
from timeit import default_timer as timer

start = timer()

# Završni rad - topologija 4x10x14L
hc_var_8 = Grillage(18.5, 18.2, 4, 10)

# Lista materijala
ST24 = MaterialProperty(1, 210000, 0.3, 7850, 235, "ST24")
AH32 = MaterialProperty(2, 210000, 0.3, 7850, 315, "AH32")
AH36 = MaterialProperty(3, 210000, 0.3, 7850, 355, "AH36")

hc_var_8.add_material(ST24)
hc_var_8.add_material(AH32)
hc_var_8.add_material(AH36)

# Korozijski dodatak
tc = CorrosionAddition(1, 2)
hc_var_8.add_corrosion_addition(tc)

# Beam property
beam_LC_T_900x12x360x32 = TBeamProperty(1, 900, 12, 360, 32, AH36)
beam_L1_T_900x12x360x32 = TBeamProperty(1, 900, 12, 360, 32, AH36)       # inicijalni longitudinal T beam prop
beam_L2_T_900x12x360x32 = TBeamProperty(2, 900, 12, 360, 32, AH36)
beam_L3_T_900x10x360x30 = TBeamProperty(3, 900, 10, 360, 30, AH36)
beam_L4_T_900x10x360x20 = TBeamProperty(4, 900, 10, 360, 20, AH36)
beam_LE_T_900x10x150x10 = LBeamProperty(5, 900, 10, 150, 10, AH36)       # inicijalni rubni L beam prop
beam_TC1_T_900x14x360x32 = TBeamProperty(6, 900, 14, 360, 32, AH36)
beam_TC2_T_900x10x260x30 = TBeamProperty(7, 900, 10, 260, 30, AH36)
beam_TC3_T_900x10x150x10 = TBeamProperty(8, 900, 10, 150, 10, AH36)
beam_TC4_T_900x10x150x10 = TBeamProperty(9, 900, 10, 150, 10, AH36)
beam_T1_T_900x12x360x30 = TBeamProperty(10, 900, 12, 360, 30, AH36)       # inicijalni transverse T beam prop
beam_T2_T_900x10x260x30 = TBeamProperty(11, 900, 10, 260, 30, AH36)
beam_T3_T_900x10x150x10 = TBeamProperty(12, 900, 10, 150, 10, AH36)
beam_T4_T_900x10x150x10 = TBeamProperty(13, 900, 10, 150, 10, AH36)
beam_TE_T_900x10x150x10 = LBeamProperty(14, 900, 10, 150, 10, AH36)       # inicijalni rubni L beam prop
stiff_HP_120x6 = BulbBeamProperty(15, 120, 6, AH36)                  # inicijalna ukrepa Bulb

hc_var_8.add_beam_prop(beam_LC_T_900x12x360x32)
hc_var_8.add_beam_prop(beam_L1_T_900x12x360x32)
hc_var_8.add_beam_prop(beam_L2_T_900x12x360x32)
hc_var_8.add_beam_prop(beam_L3_T_900x10x360x30)
hc_var_8.add_beam_prop(beam_L4_T_900x10x360x20)
hc_var_8.add_beam_prop(beam_LE_T_900x10x150x10)
hc_var_8.add_beam_prop(beam_TC1_T_900x14x360x32)
hc_var_8.add_beam_prop(beam_TC2_T_900x10x260x30)
hc_var_8.add_beam_prop(beam_TC3_T_900x10x150x10)
hc_var_8.add_beam_prop(beam_TC4_T_900x10x150x10)
hc_var_8.add_beam_prop(beam_T1_T_900x12x360x30)
hc_var_8.add_beam_prop(beam_T2_T_900x10x260x30)
hc_var_8.add_beam_prop(beam_T3_T_900x10x150x10)
hc_var_8.add_beam_prop(beam_T4_T_900x10x150x10)
hc_var_8.add_beam_prop(beam_TE_T_900x10x150x10)
hc_var_8.add_beam_prop(stiff_HP_120x6)

# Plate property
plateprop1 = PlateProperty(1, 11, AH36)                                 # inicijalni plate property za cijeli poklopac
plateprop2 = PlateProperty(2, 9, AH36)
hc_var_8.add_plate_prop(plateprop1)
hc_var_8.add_plate_prop(plateprop2)

# Stiffener layouts
stifflayout1 = StiffenerLayout(1, stiff_HP_120x6, DefinitionType.NUMBER, 14)
hc_var_8.add_stiffener_layout(stifflayout1)                                   # dodavanje stiffener layouta u dictionary
stiff_dir = BeamDirection.LONGITUDINAL                                        # inicijalna orijentacija ukrepa na svim zonama oplate

# Generacija topologije
hc_var_8.generate_prim_supp_members()                                   # Generacija svih jakih nosaca
hc_var_8.generate_segments(beam_L1_T_900x12x360x32, beam_T1_T_900x12x360x30, beam_LE_T_900x10x150x10)  # Generacija svih segmenata
hc_var_8.generate_plating(plateprop1, stifflayout1, stiff_dir)  # Generacija oplate

# Pridruzivanje simetricnih elemenata
hc_var_8.assign_symmetric_members()
hc_var_8.assign_symmetric_plating()
hc_var_8.assign_symmetric_segments()

# Izmjene plating property
hc_var_8.set_plating_prop_longitudinals(1, "plate_prop", plateprop2)
hc_var_8.set_plating_prop_longitudinals(19, "plate_prop", plateprop2)
hc_var_8.set_plating_prop_symmetric(10, "plate_prop", plateprop2)
hc_var_8.set_plating_prop_symmetric(11, "plate_prop", plateprop2)

# Simetricna izmjena svojstava segmenta
hc_var_8.set_long_symm_segment_beam_property(2, 1, beam_L4_T_900x10x360x20)
hc_var_8.set_long_symm_segment_beam_property(2, 2, beam_L3_T_900x10x360x30)
hc_var_8.set_long_symm_segment_beam_property(2, 3, beam_L2_T_900x12x360x32)
hc_var_8.set_long_symm_segment_beam_property(2, 4, beam_L1_T_900x12x360x32)
hc_var_8.set_long_symm_segment_beam_property(2, 5, beam_LC_T_900x12x360x32)
hc_var_8.set_tran_symm_segment_beam_property(2, 1, beam_T4_T_900x10x150x10)
hc_var_8.set_tran_symm_segment_beam_property(3, 1, beam_T3_T_900x10x150x10)
hc_var_8.set_tran_symm_segment_beam_property(4, 1, beam_T2_T_900x10x260x30)
hc_var_8.set_tran_symm_segment_beam_property(5, 1, beam_T1_T_900x12x360x30)
hc_var_8.set_tran_symm_segment_beam_property(2, 2, beam_TC4_T_900x10x150x10)
hc_var_8.set_tran_symm_segment_beam_property(3, 2, beam_TC3_T_900x10x150x10)
hc_var_8.set_tran_symm_segment_beam_property(4, 2, beam_TC2_T_900x10x260x30)
hc_var_8.set_tran_symm_segment_beam_property(5, 2, beam_TC1_T_900x14x360x32)

# Spremanje generirane topologije
GrillageModelData('../grillage savefiles/hc_var_8_savefile.gin').write_file(hc_var_8)

end = timer()
print("***************************************************************************************************************")
print("Code run time =", end - start, "s")
