"""
TESTNI MODUL

Generacija nove topologije i spremanje iste u datoteku
"""

from grillage.grillage_model import *
from timeit import default_timer as timer

start = timer()

hc_var_7 = Grillage(18.54, 18.18, 5, 4)

# Lista materijala
ST24 = MaterialProperty(1, 210000, 0.3, 7.85 * 10**(-9), 235, "ST24")
AH32 = MaterialProperty(2, 210000, 0.3, 7.85 * 10**(-9), 315, "AH32")
AH36 = MaterialProperty(3, 210000, 0.3, 7.85 * 10**(-9), 355, "AH36")

hc_var_7.add_material(ST24)
hc_var_7.add_material(AH32)
hc_var_7.add_material(AH36)

# Korozijski dodatak
tc = CorrosionAddition(1, 2)
hc_var_7.add_corrosion_addition(tc)

# Beam property
initial_longitudinal_beam = TBeamProperty(1, 900, 10, 230, 16, ST24)       # inicijalni longitudinal T beam prop
initial_transverse_beam = TBeamProperty(2, 900, 10, 545, 40, ST24)         # inicijalni transverse T beam prop
initial_edge_beam = LBeamProperty(3, 900, 10, 150, 16, ST24)               # inicijalni rubni L beam prop
initial_stiffener = FBBeamProperty(4, 250, 8, ST24)                         # inicijalna ukrepa FB

hc_var_7.add_beam_prop(initial_longitudinal_beam)
hc_var_7.add_beam_prop(initial_transverse_beam)
hc_var_7.add_beam_prop(initial_edge_beam)
hc_var_7.add_beam_prop(initial_stiffener)

# Plate property
plateprop1 = PlateProperty(1, 11, ST24)                                 # inicijalni plate property za cijeli poklopac
hc_var_7.add_plate_prop(plateprop1)

# Stiffener layouts
stifflayout1 = StiffenerLayout(1, initial_stiffener, DefinitionType.SPACING, 0.800)
stifflayout2 = StiffenerLayout(2, initial_stiffener, DefinitionType.SPACING, 0.750)
stifflayout3 = StiffenerLayout(3, initial_stiffener, DefinitionType.SPACING, 0.700)
hc_var_7.add_stiffener_layout(stifflayout1)                                   # dodavanje stiffener layouta u dictionary
hc_var_7.add_stiffener_layout(stifflayout2)
hc_var_7.add_stiffener_layout(stifflayout3)

stiff_dir = BeamDirection.LONGITUDINAL                                        # inicijalna orijentacija ukrepa na svim zonama oplate

# Generacija topologije
hc_var_7.generate_prim_supp_members()                                   # Generacija svih jakih nosaca
hc_var_7.generate_segments(initial_longitudinal_beam, initial_transverse_beam, initial_edge_beam)  # Generacija svih segmenata
hc_var_7.generate_plating(plateprop1, stifflayout1, stiff_dir)  # Generacija oplate

# Pridruzivanje simetricnih elemenata
hc_var_7.assign_symmetric_members()
hc_var_7.assign_symmetric_plating()
hc_var_7.assign_symmetric_segments()

# Izmjena polozaja jakih nosaca
hc_var_7.set_all_longitudinal_PSM(4.5, 4.59, 4.59)
hc_var_7.set_all_transverse_PSM(6.335, 5.870)

# Izmjene plating property
hc_var_7.set_plating_prop_symmetric(4, "stiff_layout", stifflayout2)
hc_var_7.set_plating_prop_symmetric(5, "stiff_layout", stifflayout3)
hc_var_7.set_plating_prop_symmetric(4, "stiff_dir", BeamDirection.TRANSVERSE)
hc_var_7.set_plating_prop_symmetric(5, "stiff_dir", BeamDirection.TRANSVERSE)


# Simetricna izmjena svojstava segmenta - FB za nosac 2 segment 2 - TEST ZA V2 *****************
# hc_var_5.set_long_symm_segment_beam_property(2, 2, flatbar_property)


# Spremanje generirane topologije
GrillageModelData('../grillage savefiles/hc_var_7_savefile.gin').write_file(hc_var_7)

end = timer()
print("***************************************************************************************************************")
print("Code run time =", end - start, "s")
