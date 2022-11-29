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
initial_longitudinal_beam = TBeamProperty(1, 500, 8, 250, 15, AH36)       # inicijalni longitudinal T beam prop
initial_transverse_beam = TBeamProperty(2, 500, 8, 250, 15, AH36)         # inicijalni transverse T beam prop
initial_edge_beam = LBeamProperty(3, 500, 8, 150, 15, AH36)               # inicijalni rubni L beam prop
initial_stiffener = BulbBeamProperty(4, 200, 10, AH36)                    # inicijalna ukrepa Bulb

hc_var_8.add_beam_prop(initial_longitudinal_beam)
hc_var_8.add_beam_prop(initial_transverse_beam)
hc_var_8.add_beam_prop(initial_edge_beam)
hc_var_8.add_beam_prop(initial_stiffener)

# Plate property
plateprop1 = PlateProperty(1, 11, AH36)                                 # inicijalni plate property za cijeli poklopac
hc_var_8.add_plate_prop(plateprop1)

# Stiffener layouts
stifflayout1 = StiffenerLayout(1, initial_stiffener, DefinitionType.NUMBER, 14)
hc_var_8.add_stiffener_layout(stifflayout1)                                   # dodavanje stiffener layouta u dictionary
stiff_dir = BeamDirection.LONGITUDINAL                                        # inicijalna orijentacija ukrepa na svim zonama oplate

# Generacija topologije
hc_var_8.generate_prim_supp_members()                                   # Generacija svih jakih nosaca
hc_var_8.generate_segments(initial_longitudinal_beam, initial_transverse_beam, initial_edge_beam)  # Generacija svih segmenata
hc_var_8.generate_plating(plateprop1, stifflayout1, stiff_dir)  # Generacija oplate

# Spremanje generirane topologije
GrillageModelData('../grillage savefiles/hc_var_8_savefile.gin').write_file(hc_var_8)

end = timer()
print("***************************************************************************************************************")
print("Code run time =", end - start, "s")
