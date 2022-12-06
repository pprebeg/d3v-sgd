"""
TESTNI MODUL

Generacija nove topologije i spremanje iste u datoteku
"""

from grillage.grillage_model import *
from timeit import default_timer as timer

start = timer()

# Eksplicitno zadana topologija hc_var_6 za provjeru:  18.54 x 18.18m,  mreza nosaca 5 x 5
# Poklopac je obostrano simetričan za test identifikacije segmenata na simetralnim osima

hc_var_6 = Grillage(18.54, 18.18, 5, 5)

# Lista materijala
ST24 = MaterialProperty(1, 210000, 0.3, 7.85 * 10**(-9), 235, "ST24")
AH32 = MaterialProperty(2, 210000, 0.3, 7.85 * 10**(-9), 315, "AH32")
AH36 = MaterialProperty(3, 210000, 0.3, 7.85 * 10**(-9), 355, "AH36")

hc_var_6.add_material(ST24)
hc_var_6.add_material(AH32)
hc_var_6.add_material(AH36)

# Korozijski dodatak
tc = CorrosionAddition(1, 2)
hc_var_6.add_corrosion_addition(tc)

# Beam property
initial_longitudinal_beam = TBeamProperty(1, 1089, 10, 230, 16, ST24)       # inicijalni longitudinal T beam prop
initial_transverse_beam = TBeamProperty(2, 1089, 10, 545, 40, ST24)         # inicijalni transverse T beam prop
initial_edge_beam = LBeamProperty(3, 1089, 10, 150, 16, ST24)               # inicijalni rubni L beam prop
initial_stiffener = HatBeamProperty(4, 220, 6, 220, 80, AH36)               # inicijalna ukrepa

hc_var_6.add_beam_prop(initial_longitudinal_beam)
hc_var_6.add_beam_prop(initial_transverse_beam)
hc_var_6.add_beam_prop(initial_edge_beam)
hc_var_6.add_beam_prop(initial_stiffener)


# Plate property
plateprop1 = PlateProperty(1, 10, ST24)                                 # inicijalni plate property za cijeli poklopac
hc_var_6.add_plate_prop(plateprop1)

# Stiffener layouts
stifflayout1 = StiffenerLayout(1, initial_stiffener, DefinitionType.NUMBER, 5)
stifflayout2 = StiffenerLayout(2, initial_stiffener, DefinitionType.SPACING, 0.7)
stifflayout3 = StiffenerLayout(3, initial_stiffener, DefinitionType.SPACING, 0.64)
hc_var_6.add_stiffener_layout(stifflayout1)                                   # dodavanje stiffener layouta u dictionary
hc_var_6.add_stiffener_layout(stifflayout2)
hc_var_6.add_stiffener_layout(stifflayout3)

stiff_dir = BeamDirection.TRANSVERSE                                        # inicijalna orijentacija ukrepa na svim zonama oplate

# Generacija topologije
hc_var_6.generate_prim_supp_members()                                   # Generacija svih jakih nosaca
hc_var_6.generate_segments(initial_longitudinal_beam, initial_transverse_beam, initial_edge_beam)  # Generacija svih segmenata
hc_var_6.generate_plating(plateprop1, stifflayout1, stiff_dir)  # Generacija oplate

# Pridruzivanje simetricnih elemenata
hc_var_6.assign_symmetric_members()
hc_var_6.assign_symmetric_plating()
hc_var_6.assign_symmetric_segments()

# Izmjena polozaja jakih nosaca
Grillage.set_all_longitudinal_PSM(hc_var_6, 5, 4.09, 4.09)
Grillage.set_all_transverse_PSM(hc_var_6, 5, 4.27, 4.27)

# Izmjene plating property - postavljanje drugačijih svojstava svim poljima oplate između poprečnih nosača koji definiraju polja oplate 1 i 5
Grillage.set_plating_prop_symmetric(hc_var_6, 2, "stiff_layout", stifflayout2)
Grillage.set_plating_prop_symmetric(hc_var_6, 2, "stiff_dir", BeamDirection.LONGITUDINAL)
Grillage.set_plating_prop_symmetric(hc_var_6, 6, "stiff_layout", stifflayout3)
Grillage.set_plating_prop_symmetric(hc_var_6, 6, "stiff_dir", BeamDirection.LONGITUDINAL)


# Spremanje generirane topologije
GrillageModelData('../grillage savefiles/hc_var_6_savefile.gin').write_file(hc_var_6)

end = timer()
print("***************************************************************************************************************")
print("Code run time =", end - start, "s")
