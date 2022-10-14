"""
Modul za testiranje diskretizacije učitanog modela

"""

from femdir.grillage_mesh import *
from timeit import default_timer as timer

start = timer()

# Ucitavanje topologije iz datoteke
filename = '../grillage savefiles\\hc_var_2_savefile.txt'
hc_variant = GrillageModelData(filename).read_file()

# Generacija mreze
mesh1 = GrillageMesh(hc_variant)    # Izrada GrillageMesh objekta
mesh1.flange_aspect_ratio = 7       # Postavljanje aspektnog odnosa elemenata prirubnica jakih nosača i oplate uz struk jakih nosača
mesh1.plate_aspect_ratio = 4        # Postavljanje aspektnog odnosa elemenata oplate i strukova jakih nosača

"""
# mesh1.GenerateNodes(hc_variant)               # Generacija cvorova
# mesh1.GenerateElements(hc_variant)            # Generacija elemenata
# Zapis mreze u datoteku
# MeshData("mesh1_savefile.txt").write_file(mesh1)
"""


def Test_element_size_stiffener_spacing(grillage, plate_id):
    plate = grillage.plating()[plate_id]
    dims = GrillageMesh.element_size_stiffener_spacing(mesh1, plate)

    print("Dimenzije quad elementa oplate uz jedan element između ukrepa, uzdužno:", dims[0], "mm, poprečno:", dims[1], "mm")


def Test_element_size_flange_width(grillage, direction: BeamDirection, psm_id, segment_id, aspect_ratio):
    segment = None

    if direction == BeamDirection.LONGITUDINAL:
        segment = grillage.longitudinal_members()[psm_id].segments[segment_id]

    elif direction == BeamDirection.TRANSVERSE:
        segment = grillage.transverse_members()[psm_id].segments[segment_id]

    dim = GrillageMesh.element_size_flange_width(grillage, segment, aspect_ratio)
    print("Maksimalna dimenzija elementa prirubnice prema aspektnom odnosu:", dim, "mm")


def Test_element_size_plating_zone(grillage, plate_id):
    plate = grillage.plating()[plate_id]
    dims = GrillageMesh.element_size_plating_zone(mesh1, grillage, plate)

    dim_x = dims[0]
    dim_y = dims[1]
    print(dim_x, dim_y)


def Test_ALL_element_size_plating_zone(grillage):
    n_x = int(grillage.N_transverse - 1)    # Broj polja u uzduznom smjeru
    n_y = int(grillage.N_longitudinal - 1)  # Broj polja u poprecnom smjeru

    polje1 = np.zeros((n_y, n_x))
    # dimenzije karakteristicnih elemenata na cijeloj oplati po x osi
    plate_id = 1
    for stupac in range(0, n_y):
        for redak in range(0, n_x):
            plate = grillage.plating()[plate_id]
            polje1[stupac, redak] = GrillageMesh.element_size_plating_zone(mesh1, grillage, plate)[0]
            plate_id += 1
    print("Sve x dimenzije elemenata: \n", polje1, "\n")

    # dimenzije karakteristicnih elemenata na cijeloj oplati po y osi
    plate_id = 1
    for stupac in range(0, n_y):
        for redak in range(0, n_x):
            plate = grillage.plating()[plate_id]
            polje1[stupac, redak] = GrillageMesh.element_size_plating_zone(mesh1, grillage, plate)[1]
            plate_id += 1
    print("Sve y dimenzije elemenata: \n", polje1)


def Test_element_size_mesh(grillage):
    GrillageMesh.element_size_mesh(mesh1, grillage)
    print("Konačno odabrane dimenzije mreže po x:", mesh1.mesh_dim_x)
    print("Konačno odabrane dimenzije mreže po y:", mesh1.mesh_dim_y)


def Test_get_mesh_dim_x(grillage, plate_id):
    GrillageMesh.element_size_mesh(mesh1, grillage)
    plate = grillage.plating()[plate_id]
    print(GrillageMesh.get_mesh_dim_x(mesh1, plate))


def Test_get_mesh_dim_y(grillage, plate_id):
    GrillageMesh.element_size_mesh(mesh1, grillage)
    plate = grillage.plating()[plate_id]
    print(GrillageMesh.get_mesh_dim_y(mesh1, plate))


def Test_all_plating_zones_mesh_dimensions(grillage):
    GrillageMesh.element_size_mesh(mesh1, grillage)
    for plate in grillage.plating().values():
        dim_x = GrillageMesh.get_mesh_dim_x(mesh1, plate)
        dim_y = GrillageMesh.get_mesh_dim_y(mesh1, plate)
        print("Zona oplate ID:", plate.id, ",   dim_x =", "{:.1f}".format(dim_x), "mm", ",   dim_y =", "{:.1f}".format(dim_y), "mm")


"""
def Test_get_number_of_quads(grillgeo, plate_id):
    GrillageMesh.element_size_mesh(mesh1, grillgeo)

    plate = grillgeo.plating()[plate_id]
    GrillageMesh.get_number_of_quads(mesh1, grillgeo, plate)
"""

"""
def TestGenerateNodes():
    for node in FEMesh.nodes(mesh1).keys():
        print("Node ID:", node, " Koordinate:", mesh1.nodes()[node].coords)


def TestCheckNodeOverlap(grillgeo):
    FEMesh.GenerateNodes(mesh1, grillgeo)
    FEMesh.CheckNodeOverlap(mesh1)


def MeshTest():
    # Provjera izrađenih elemenata
    for i in FEMesh.elements(mesh1):
        element = FEMesh.elements(mesh1)[i]
        print("ID Elementa:", element.id, ", property ID:", element.property_id,
              ", ID cvorova:", element.node1.id, element.node2.id, element.node3.id, element.node4.id,
              ", koordinate cvorova:", element.node1.coords, element.node2.coords, element.node3.coords, element.node4.coords, )


def TestPlatingNodes():
    nodes_list = FEMesh.plating_nodes(mesh1)
    n_polja = len(nodes_list)
    print("Upisano polja:", n_polja)
    for polje_id in range(0, n_polja):
        print("Polje", polje_id + 1, ": \n", nodes_list[polje_id])
"""

# print(FEMesh.find_closest_divisor(4935, 660))
# Test_element_size_stiffener_spacing(hc_variant, 1)   # Dimenzije elementa oplate za odabranu zonu prema razmaku ukrepa
# Test_element_size_flange_width(hc_variant, BeamOrientation.TRANSVERSE, 1, 1)    # Dimenzije elementa prirubnice za odabrani segment
# Test_element_size_plating_zone(hc_variant, 1)         # Dimenzije elemenata na odabranoj zoni oplate prema razmaku ukrepa i ar prirubnice
# Test_ALL_element_size_plating_zone(hc_variant)
# Test_element_size_mesh(hc_variant)                    # Konacno odabrane dimenzije mreze po x i y
# Test_get_mesh_dim_x(hc_variant, 2)                    # Odabrana x dimenzija za neko polje oplate
# Test_get_mesh_dim_y(hc_variant, 2)                    # Odabrana y dimenzija za neko polje oplate
# Test_all_plating_zones_mesh_dimensions(hc_variant)    # Odabrane x i y dimenzije elemenata za sva polja oplate

# ***** WIP ******
# Test_get_number_of_quads(hc_variant, 2)
# TestGenerateNodes()                           # Ispis cvorova i njihovih koordinata
# TestCheckNodeOverlap(hc_variant)
# MeshTest()                                    # Ispis elemenata, property ID, cvorova elemenata i koordinata
# TestPlatingNodes()                            # Matricni zapis svih cvorova oplate
# np.allclose()

end = timer()

print("***************************************************************************************************************")
print("Code run time =", end - start, "s")
