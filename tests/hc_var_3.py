"""
TESTNI MODUL

Generacija nove topologije i spremanje iste u datoteku
"""

from grillage.grillage_model import *
from timeit import default_timer as timer

start = timer()

# Eksplicitno zadana topologija hc_var_3 za provjeru:  18.54m x 18.18m,  mreza nosaca 5 x 5
hc_var_3 = Grillage(18.54, 18.18, 5, 5)

# Lista materijala
ST24 = MaterialProperty(1, 210000, 0.3, 7.85 * 10**(-9), 235, "ST24")
AH32 = MaterialProperty(2, 210000, 0.3, 7.85 * 10**(-9), 315, "AH32")
AH36 = MaterialProperty(3, 210000, 0.3, 7.85 * 10**(-9), 355, "AH36")

hc_var_3.add_material(ST24)
hc_var_3.add_material(AH32)
hc_var_3.add_material(AH36)

# Korozijski dodatak
tc = CorrosionAddition(1, 2)
hc_var_3.add_corrosion_addition(tc)

# hc_var_3.add_corrosion_addition(1, 2)     # ALTERNATIVNO ZADAVANJE I SPREMANJE

# Beam property
initial_longitudinal_beam = TBeamProperty(1, 1089, 10, 230, 16, ST24)       # inicijalni longitudinal T beam prop
initial_transverse_beam = TBeamProperty(2, 1089, 10, 545, 40, ST24)         # inicijalni transverse T beam prop
initial_edge_beam = LBeamProperty(3, 1089, 10, 150, 16, ST24)               # inicijalni rubni L beam prop
initial_stiffener = HatBeamProperty(4, 220, 6, 220, 80, AH36)               # inicijalna ukrepa
center_girder = TBeamProperty(5, 1089, 10, 560, 40, ST24)

hc_var_3.add_beam_prop(initial_longitudinal_beam)
hc_var_3.add_beam_prop(initial_transverse_beam)
hc_var_3.add_beam_prop(initial_edge_beam)
hc_var_3.add_beam_prop(initial_stiffener)
hc_var_3.add_beam_prop(center_girder)

# Plate property
plateprop1 = PlateProperty(1, 10, ST24)                                 # inicijalni plate property za cijeli poklopac
hc_var_3.add_plate_prop(plateprop1)

# Stiffener layouts
stifflayout1 = StiffenerLayout(1, initial_stiffener, DefinitionType.SPACING, 0.873)
stifflayout2 = StiffenerLayout(2, initial_stiffener, DefinitionType.NUMBER, 2)
hc_var_3.add_stiffener_layout(stifflayout1)                                   # dodavanje stiffener layouta u dictionary
hc_var_3.add_stiffener_layout(stifflayout2)

stiff_dir = BeamDirection.TRANSVERSE                                        # inicijalna orijentacija ukrepa na svim zonama oplate

# Generacija topologije
hc_var_3.generate_prim_supp_members()                                   # Generacija svih jakih nosaca
hc_var_3.generate_segments(initial_longitudinal_beam, initial_transverse_beam, initial_edge_beam)  # Generacija svih segmenata
hc_var_3.generate_plating(plateprop1, stifflayout1, stiff_dir)  # Generacija oplate

# Pridruzivanje simetricnih elemenata
hc_var_3.assign_symmetric_members()
hc_var_3.assign_symmetric_plating()

# Izmjena polozaja jakih nosaca
Grillage.set_all_longitudinal_PSM(hc_var_3, 4.5, 4.59, 4.59)
Grillage.set_all_transverse_PSM(hc_var_3, 4.325, 4.935, 4.935)

# Izmjene plating property
Grillage.set_plating_prop_transversals(hc_var_3, 1, "stiff_layout", stifflayout2)
Grillage.set_plating_prop_transversals(hc_var_3, 4, "stiff_layout", stifflayout2)

Grillage.set_plating_prop_symmetric(hc_var_3, 2, "stiff_dir", BeamDirection.LONGITUDINAL)
Grillage.set_plating_prop_symmetric(hc_var_3, 2, "stiff_layout", stifflayout2)

# Izmjena svojstava nosaca
Grillage.set_tran_member_beam_property(hc_var_3, 3, center_girder)

# Spremanje generirane topologije
GrillageModelData('../grillage savefiles/hc_var_3_savefile.gin').write_file(hc_var_3)

end = timer()
print("***************************************************************************************************************")
print("Code run time =", end - start, "s")
