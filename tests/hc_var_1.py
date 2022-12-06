"""
TESTNI MODUL

Generacija nove topologije i spremanje iste u datoteku
"""

from grillage.grillage_model import *
from timeit import default_timer as timer

start = timer()

# Eksplicitno zadana topologija hc_var_1 za provjeru:  18.54m x 18.18m,  mreza nosaca 5 x 5
hc_var_1 = Grillage(18.54, 18.18, 5, 5)

# Lista materijala
ST24 = MaterialProperty(1, 210000, 0.3, 7.85 * 10**(-9), 235, "ST24")
AH32 = MaterialProperty(2, 210000, 0.3, 7.85 * 10**(-9), 315, "AH32")
AH36 = MaterialProperty(3, 210000, 0.3, 7.85 * 10**(-9), 355, "AH36")

hc_var_1.add_material(ST24)
hc_var_1.add_material(AH32)
hc_var_1.add_material(AH36)

# Korozijski dodatak
tc = CorrosionAddition(1, 2)
hc_var_1.add_corrosion_addition(tc)

# Beam property
initial_longitudinal_beam = TBeamProperty(1, 1089, 10, 230, 16, ST24)       # inicijalni longitudinal T beam prop
initial_transverse_beam = TBeamProperty(2, 1089, 10, 545, 40, ST24)         # inicijalni transverse T beam prop
initial_edge_beam = LBeamProperty(3, 1089, 12, 150, 16, ST24)               # inicijalni rubni L beam prop
initial_stiffener = HatBeamProperty(4, 300, 6, 220, 80, AH36)             # inicijalna ukrepa Hat
# initial_stiffener = BulbBeamProperty(4, 240, 10, AH36)                    # inicijalna ukrepa Bulb
# initial_stiffener = TBeamProperty(4, 250, 8, 90, 12, ST24)                # inicijalna ukrepa T
# initial_stiffener = FBBeamProperty(4, 250, 8, ST24)                       # inicijalna ukrepa FB
# stiffener2 = BulbBeamProperty(6, 240, 10, AH36)                    # inicijalna ukrepa Bulb


center_girder = TBeamProperty(5, 1089, 10, 560, 40, AH32)
center_longitudinal = TBeamProperty(6, 1089, 10, 305, 16, AH32)
center_long_ends = TBeamProperty(7, 1089, 12, 305, 16, AH32)
long_edge_beam = LBeamProperty(8, 1089, 10, 150, 16, ST24)
FB_beam = FBBeamProperty(9, 1089, 10, ST24)

# DRUGI SET - FB JAKI NOSAČI
# initial_longitudinal_beam = FBBeamProperty(1, 1089, 10, ST24)
# initial_transverse_beam = FBBeamProperty(2, 1089, 10, ST24)
# initial_edge_beam = FBBeamProperty(3, 1089, 10, ST24)
# initial_stiffener = HatBeamProperty(4, 220, 6, 220, 80, AH36)               # inicijalna ukrepa
# center_girder = FBBeamProperty(5, 1089, 10, ST24)

hc_var_1.add_beam_prop(initial_longitudinal_beam)
hc_var_1.add_beam_prop(initial_transverse_beam)
hc_var_1.add_beam_prop(initial_edge_beam)
hc_var_1.add_beam_prop(initial_stiffener)
hc_var_1.add_beam_prop(center_girder)
hc_var_1.add_beam_prop(center_longitudinal)
hc_var_1.add_beam_prop(center_long_ends)
hc_var_1.add_beam_prop(long_edge_beam)
hc_var_1.add_beam_prop(FB_beam)
# hc_var_1.add_beam_prop(stiffener2)

# Plate property
plateprop1 = PlateProperty(1, 10, ST24)     # inicijalni plate property za cijeli poklopac
plateprop2 = PlateProperty(2, 11, AH32)
plateprop3 = PlateProperty(3, 9, ST24)
plateprop4 = PlateProperty(4, 9, ST24)
hc_var_1.add_plate_prop(plateprop1)
hc_var_1.add_plate_prop(plateprop2)     # Dodatni plate property za testiranje identifikacije jednakih svojstava
hc_var_1.add_plate_prop(plateprop3)
hc_var_1.add_plate_prop(plateprop4)

# Stiffener layouts
stifflayout1 = StiffenerLayout(1, initial_stiffener, DefinitionType.SPACING, 0.873)
stifflayout2 = StiffenerLayout(2, initial_stiffener, DefinitionType.SPACING, 0.935)
stifflayout3 = StiffenerLayout(3, initial_stiffener, DefinitionType.SPACING, 1.500)

hc_var_1.add_stiffener_layout(stifflayout1)                                   # dodavanje stiffener layouta u dictionary
hc_var_1.add_stiffener_layout(stifflayout2)
hc_var_1.add_stiffener_layout(stifflayout3)


stiff_dir = BeamDirection.TRANSVERSE                  # inicijalna orijentacija ukrepa na svim zonama oplate

# Generacija topologije
hc_var_1.generate_prim_supp_members()                 # Generacija svih jakih nosaca
hc_var_1.generate_segments(initial_longitudinal_beam, initial_transverse_beam, initial_edge_beam)  # Generacija svih segmenata
hc_var_1.generate_plating(plateprop1, stifflayout1, stiff_dir)  # Generacija oplate
hc_var_1.generate_elementary_plate_panels()                     # Generacija neukrepljenih (elementarnih) polja oplate

# hc_var_1.plating()[6].set_intercostal_stiffeners(4, FB_beam)    # Dodavanje interkostalnih ukrepa na sva neukrepljena polja zone 6
# hc_var_1.plating()[7].set_intercostal_stiffeners(4, FB_beam)
# hc_var_1.plating()[10].set_intercostal_stiffeners(4, FB_beam)
# hc_var_1.plating()[11].set_intercostal_stiffeners(4, FB_beam)

# hc_var_1.plating()[1].elementary_plate_panels[1].intercostal_stiffener_num = 1    # Dodavanje interkostalnih ukrepa na jedno elementarno polje
# hc_var_1.plating()[1].elementary_plate_panels[1].beam_prop = FB_beam

# hc_var_1.plating()[6].regenerate_elementary_plate_panel()       # Regeneracija - brisanje svih interkostala sa zone oplate 6

# Pridruzivanje simetricnih elemenata
hc_var_1.assign_symmetric_members()
hc_var_1.assign_symmetric_plating()
hc_var_1.assign_symmetric_segments()

# Izmjena polozaja jakih nosaca
Grillage.set_all_longitudinal_PSM(hc_var_1, 4.5, 4.59, 4.59)
Grillage.set_all_transverse_PSM(hc_var_1, 4.325, 4.935, 4.935)

# Izmjene plating property - postavljanje drugačijih svojstava svim poljima oplate između poprečnih nosača koji definiraju polja oplate 1 i 4
hc_var_1.set_plating_prop_transversals(1, "stiff_layout", stifflayout2)
hc_var_1.set_plating_prop_transversals(4, "stiff_layout", stifflayout2)
hc_var_1.set_plating_prop_transversals(2, "plate_prop", plateprop2)
hc_var_1.set_plating_prop_transversals(3, "plate_prop", plateprop2)

# Grillage.set_plating_prop_symmetric(hc_var_1, 1, "stiff_dir", BeamDirection.LONGITUDINAL)
# hc_var_1.plating()[6].stiff_layout = stifflayout3                             # Unos drugacijeg layouta da mesh_feasibility javi gresku
# hc_var_1.plating()[7].stiff_layout = stifflayout3                             # Unos drugacijeg layouta da mesh_feasibility javi gresku
# hc_var_1.plating()[10].stiff_layout = stifflayout3                             # Unos drugacijeg layouta da mesh_feasibility javi gresku
# hc_var_1.plating()[11].stiff_layout = stifflayout3                             # Unos drugacijeg layouta da mesh_feasibility javi gresku

# hc_var_1.plating()[1].stiff_dir = BeamDirection.LONGITUDINAL                  # izmjena da zone oplate ne budu simetricne
# hc_var_1.longitudinal_members()[1].segments[1].beam_prop = center_girder      # izmjena da segmenti ne budu simetricni
# hc_var_1.set_tran_member_beam_property(1, FB_beam)                            # Izmjena beam property za cijeli nosač

# Izmjena svojstava nosaca
hc_var_1.set_tran_member_beam_property(3, center_girder)
hc_var_1.set_long_member_beam_property(3, center_longitudinal)
hc_var_1.longitudinal_members()[3].segments[0].beam_prop = center_long_ends
hc_var_1.longitudinal_members()[3].segments[3].beam_prop = center_long_ends
# hc_var_1.longitudinal_members()[3].segments[1].beam_prop = FB_beam        # test za mesh feasibility u slucaju da je jedan segment FB

hc_var_1.set_long_member_beam_property(1, long_edge_beam)
hc_var_1.set_long_member_beam_property(5, long_edge_beam)


# Spremanje generirane topologije
GrillageModelData('../grillage savefiles/hc_var_1_savefile.gin').write_file(hc_var_1)

end = timer()
print("***************************************************************************************************************")
print("Code run time =", end - start, "s")
