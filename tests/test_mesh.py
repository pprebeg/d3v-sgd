"""
Modul za prikaz diskretizaciju ucitanog modela

"""

# from FEM_mesh import *
from grillage.FEM_mesh import *
from timeit import default_timer as timer

start = timer()

# Ucitavanje topologije iz datoteke
filename = '..\\examples\\hc_var_1.txt'
hc_variant = GrillageModelData(filename).read_file()

# Generacija mreze
mesh1 = FEMesh(hc_variant)                  # Izrada FEMesh objekta

"""
# mesh1.GenerateNodes(hc_variant)               # Generacija cvorova
# mesh1.GenerateElements(hc_variant)            # Generacija elemenata
# Zapis mreze u datoteku
# MeshData("mesh1_savefile.txt").write_file(mesh1)
"""


def Test_element_size_stiffener_spacing(grillage, plate_id):
    plate = grillage.plating()[plate_id]
    dims = FEMesh.element_size_stiffener_spacing(mesh1, plate)

    print("Dimenzije quad elementa oplate uz jedan element između ukrepa, uzdužno:", dims[0], "mm, poprečno:", dims[1], "mm")


def Test_element_size_flange_width(grillage, direction: BeamOrientation, psm_id, segment_id):
    segment = None
    if direction == BeamOrientation.LONGITUDINAL:
        segment = grillage.longitudinal_members()[psm_id].segments[segment_id]
    elif direction == BeamOrientation.TRANSVERSE:
        segment = grillage.longitudinal_members()[psm_id].segments[segment_id]

    dim = FEMesh.element_size_flange_width(grillage, segment)
    print("Maksimalna dimenzija elementa prirubnice prema aspektnom odnosu:", dim, "mm")


def Test_element_size_plating_zone(grillage, plate_id):
    plate = grillage.plating()[plate_id]
    dims = FEMesh.element_size_plating_zone(mesh1, grillage, plate)

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
            polje1[stupac, redak] = FEMesh.element_size_plating_zone(mesh1, grillage, plate)[0]
            plate_id += 1
    print("Sve x dimenzije elemenata: \n", polje1, "\n")

    # dimenzije karakteristicnih elemenata na cijeloj oplati po y osi
    plate_id = 1
    for stupac in range(0, n_y):
        for redak in range(0, n_x):
            plate = grillage.plating()[plate_id]
            polje1[stupac, redak] = FEMesh.element_size_plating_zone(mesh1, grillage, plate)[1]
            plate_id += 1
    print("Sve y dimenzije elemenata: \n", polje1)


def Test_element_size_mesh(grillage):
    FEMesh.element_size_mesh(mesh1, grillage)
    print("Konačno odabrane dimenzije mreže po x:", mesh1.mesh_dim_x)
    print("Konačno odabrane dimenzije mreže po y:", mesh1.mesh_dim_y)


def Test_get_mesh_dim_x(grillage, plate_id):
    FEMesh.element_size_mesh(mesh1, grillage)
    plate = grillage.plating()[plate_id]
    print(FEMesh.get_mesh_dim_x(mesh1, plate))


def Test_get_mesh_dim_y(grillage, plate_id):
    FEMesh.element_size_mesh(mesh1, grillage)
    plate = grillage.plating()[plate_id]
    print(FEMesh.get_mesh_dim_y(mesh1, plate))


def Test_get_number_of_quads(grillage, plate_id):
    FEMesh.element_size_mesh(mesh1, grillage)

    plate = grillage.plating()[plate_id]
    FEMesh.get_number_of_quads(mesh1, grillage, plate)


"""
def TestGenerateNodes():
    for node in FEMesh.nodes(mesh1).keys():
        print("Node ID:", node, " Koordinate:", mesh1.nodes()[node].coords)


def TestCheckNodeOverlap(grillage):
    FEMesh.GenerateNodes(mesh1, grillage)
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
# Test_element_size_stiffener_spacing(hc_variant, 2)                                # Dimenzije elementa oplate za odabranu zonu
# Test_element_size_flange_width(hc_variant, BeamOrientation.LONGITUDINAL, 1, 1)    # Dimenzije elementa prirubnice za odabrani segment
# Test_element_size_plating_zone(hc_variant, 1)
# Test_ALL_element_size_plating_zone(hc_variant)
Test_element_size_mesh(hc_variant)                                                # Konacno odabrane dimenzije mreze po x i y
# Test_get_mesh_dim_x(hc_variant, 2)                                                # Odabrana x dimenzija za neko polje oplate
# Test_get_mesh_dim_y(hc_variant, 2)                                                # Odabrana y dimenzija za neko polje oplate
# Test_get_number_of_quads(hc_variant, 2)

# TestGenerateNodes()                           # Ispis cvorova i njihovih koordinata
# TestCheckNodeOverlap(hc_variant)
# MeshTest()                                    # Ispis elemenata, property ID, cvorova elemenata i koordinata
# TestPlatingNodes()                            # Matricni zapis svih cvorova oplate
# np.allclose()

end = timer()

print("***************************************************************************************************************")
print("Code run time =", end - start, "s")
