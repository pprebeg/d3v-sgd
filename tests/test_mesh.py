"""
Modul za testiranje diskretizacije učitanog modela

"""

from femdir.grillage_mesh import *
from timeit import default_timer as timer

start = timer()

# Ucitavanje topologije iz datoteke
filename = '../grillage savefiles/hc_var_1_savefile.gin'
hc_variant = GrillageModelData(filename).read_file()

# Generacija mreze
mesh1 = GrillageMesh(hc_variant, MeshSolution.V1)    # Izrada GrillageMesh objekta
mesh1.min_num_ebs = 1               # Postavljanje minimalnog broja elemenata između ukrepa
# mesh1.min_num_eweb = 3              # Postavljanje minimalnog broja elemenata duž visine struka
mesh1.flange_aspect_ratio = 7       # Postavljanje aspektnog odnosa elemenata prirubnica jakih nosača i oplate uz struk jakih nosača
mesh1.plate_aspect_ratio = 4        # Postavljanje aspektnog odnosa elemenata oplate i strukova jakih nosača

"""
# mesh1.GenerateNodes(hc_variant)               # Generacija cvorova
# mesh1.GenerateElements(hc_variant)            # Generacija elemenata
# Zapis mreze u datoteku
# MeshData("mesh1_savefile.txt").write_file(mesh1)
"""


def Test_get_reduced_plate_dim(grillage, plate_id):
    plate = grillage.plating()[plate_id]
    reduced_dim = GrillageMesh.get_reduced_plate_dim(mesh1, plate)
    print("Reducirana dimenzija paralalna s ukrepama na zoni oplate", plate_id, "iznosi:", reduced_dim)


def Test_element_size_stiffener_spacing(grillage, plate_id):
    plate = grillage.plating()[plate_id]
    dims = GrillageMesh.element_size_stiffener_spacing(mesh1, plate)
    n_elem = mesh1.min_num_ebs
    print("Dimenzije quad elementa oplate uz", n_elem, "element između ukrepa, uzdužno:", dims[0], "mm, poprečno:", dims[1], "mm")


def Test_get_flange_el_length(grillage, direction: BeamDirection, psm_id, segment_id):
    segment = None

    if direction == BeamDirection.LONGITUDINAL:
        segment = grillage.longitudinal_members()[psm_id].segments[segment_id]

    elif direction == BeamDirection.TRANSVERSE:
        segment = grillage.transverse_members()[psm_id].segments[segment_id]

    dim = GrillageMesh.get_flange_el_length(mesh1, segment)
    print("Maksimalna dimenzija elementa prirubnice prema aspektnom odnosu:", dim, "mm")


def Test_element_size_plating_zone(grillage, plate_id):
    plate = grillage.plating()[plate_id]
    dims = GrillageMesh.element_size_plating_zone(mesh1, plate)
    print("Osnovne dimenzije elementa na zoni oplate", plate_id)
    dim_x = dims[0]
    dim_y = dims[1]
    print("dim_x =", dim_x, ", dim_y =", dim_y)


def Test_element_aspect_ratio(dim_x, dim_y):
    ar = GrillageMesh.element_aspect_ratio(dim_x, dim_y)
    print("Aspektni odnos =", ar[0], ", gdje je veća dimenzija", ar[1])


def Test_refine_plate_element(length, dim_limit, perpendicular_dim):
    dim = GrillageMesh.refine_plate_element(mesh1, length, dim_limit, perpendicular_dim)
    print("Duljinu", length, "je potrebno podijeliti na jednake dijelove, tako da dimenzija elementa ne prelazi", dim_limit,
          ". \n Odabrana je dimenzija", dim, "koja istovremeno zadovoljava aspektni odnos oplate uz drugu dimenziju elementa",
          perpendicular_dim)


def Test_ALL_element_size_plating_zone(grillage):
    n_x = int(grillage.N_transverse - 1)    # Broj polja u uzduznom smjeru
    n_y = int(grillage.N_longitudinal - 1)  # Broj polja u poprecnom smjeru

    polje1 = np.zeros((n_y, n_x))
    # dimenzije karakteristicnih elemenata na cijeloj oplati po x osi
    plate_id = 1
    for stupac in range(0, n_y):
        for redak in range(0, n_x):
            plate = grillage.plating()[plate_id]
            polje1[stupac, redak] = GrillageMesh.element_size_plating_zone(mesh1, plate)[0]
            plate_id += 1
    print("Sve x dimenzije elemenata: \n", polje1, "\n")

    # dimenzije karakteristicnih elemenata na cijeloj oplati po y osi
    plate_id = 1
    for stupac in range(0, n_y):
        for redak in range(0, n_x):
            plate = grillage.plating()[plate_id]
            polje1[stupac, redak] = GrillageMesh.element_size_plating_zone(mesh1, plate)[1]
            plate_id += 1
    print("Sve y dimenzije elemenata: \n", polje1)


def Test_element_size_mesh(grillage):
    GrillageMesh.element_size_mesh(mesh1, grillage)
    print("Konačno odabrane dimenzije mreže po x:", mesh1.mesh_dim_x)
    print("Konačno odabrane dimenzije mreže po y:", mesh1.mesh_dim_y)


def Test_get_base_dim_x(grillage, plate_id):
    GrillageMesh.element_size_mesh(mesh1, grillage)
    plate = grillage.plating()[plate_id]
    print(GrillageMesh.get_base_dim_x(mesh1, plate))


def Test_get_base_dim_y(grillage, plate_id):
    GrillageMesh.element_size_mesh(mesh1, grillage)
    plate = grillage.plating()[plate_id]
    print(GrillageMesh.get_base_dim_y(mesh1, plate))


def Test_element_size_transition_x(grillage, plate_id):
    GrillageMesh.element_size_mesh(mesh1, grillage)
    plate = grillage.plating()[plate_id]
    print("Zona oplate", plate_id, ", dimenzije prijelaznih elemenata po x:", GrillageMesh.element_size_transition_x(mesh1, plate))


def Test_element_size_transition_y(grillage, plate_id):
    GrillageMesh.element_size_mesh(mesh1, grillage)
    plate = grillage.plating()[plate_id]
    print("Zona oplate", plate_id, ", dimenzije prijelaznih elemenata po y:", GrillageMesh.element_size_transition_y(mesh1, plate))


def Test_transition_element_dimensions(grillage, plate_id):
    GrillageMesh.element_size_mesh(mesh1, grillage)
    plate = grillage.plating()[plate_id]
    print("Zona oplate", plate_id, ", dimenzije prijelaznih elemenata po x:", GrillageMesh.element_size_transition_x(mesh1, plate))
    print("Zona oplate", plate_id, ", dimenzije prijelaznih elemenata po y:", GrillageMesh.element_size_transition_y(mesh1, plate))


def Test_get_flange_el_width(psm_id, segment_id):
    segment = hc_variant.longitudinal_members()[psm_id].segments[segment_id - 1]

    print("Jaki uzduzni nosac", psm_id, ", ID segmenta:", segment.id,
          " ,BeamProperty ID:", segment.beam_prop.id,
          ", tip profila: ", segment.beam_prop.beam_type,
          ", tf =", segment.beam_prop.tf,
          ", dimenzija elementa prirubnice po širini:", mesh1.get_flange_el_width(segment), "mm")


def Test_all_plating_zones_mesh_dimensions(grillage):
    GrillageMesh.element_size_mesh(mesh1, grillage)
    for plate in grillage.plating().values():
        dim_x = GrillageMesh.get_base_dim_x(mesh1, plate)
        dim_y = GrillageMesh.get_base_dim_y(mesh1, plate)
        print("Zona oplate ID:", plate.id, ",   dim_x =", "{:.2f}".format(dim_x), "mm", ",   dim_y =", "{:.2f}".format(dim_y), "mm")


def Test_identify_unique_property(grillage):
    GrillageMesh.identify_unique_property(mesh1, grillage)

    for prop in mesh1.unique_properties.values():
        plate_prop = len(prop.plate_prop)
        beam_prop = len(prop.beam_prop)
        print("Unique property ID:", prop.id, ", tp =", prop.tp, "mm, material ID", prop.mat.id,
              ", upisano istih na modelu, plate:", plate_prop, ", beam:", beam_prop)


def Test_transition_element_mesh(grillage):
    GrillageMesh.element_size_mesh(mesh1, grillage)
    GrillageMesh.transition_element_mesh(mesh1, grillage)


def Test_get_tr_dim_x(grillage, plate_id):
    plate = grillage.plating()[plate_id]
    GrillageMesh.element_size_mesh(mesh1, grillage)
    GrillageMesh.transition_element_mesh(mesh1, grillage)
    dims = GrillageMesh.get_tr_dim_x(mesh1, plate)
    print("Dimenzije x prijelaznih elemenata na zoni oplate", plate_id, ":", dims)


def Test_get_tr_dim_y(grillage, plate_id):
    plate = grillage.plating()[plate_id]
    GrillageMesh.element_size_mesh(mesh1, grillage)
    GrillageMesh.transition_element_mesh(mesh1, grillage)
    dims = GrillageMesh.get_tr_dim_y(mesh1, plate)
    print("Dimenzije y prijelaznih elemenata na zoni oplate", plate_id, ":", dims)


def Test_get_tr_element_num(grillage, plate_id):
    plate = grillage.plating()[plate_id]
    GrillageMesh.element_size_mesh(mesh1, grillage)
    GrillageMesh.transition_element_mesh(mesh1, grillage)
    n_elemenata = GrillageMesh.get_tr_element_num(mesh1, plate)
    print("Broj prijelaznih elemenata na zoni oplate", plate_id, "po x:", n_elemenata[0], ", po y:", n_elemenata[1])


def Test_get_element_number(grillage, plate_id):
    plate = grillage.plating()[plate_id]
    GrillageMesh.element_size_mesh(mesh1, grillage)
    GrillageMesh.transition_element_mesh(mesh1, grillage)
    n_elemenata = GrillageMesh.get_element_number(mesh1, plate)
    print("Broj elemenata na zoni oplate", plate_id, "po x:", n_elemenata[0], ", po y:", n_elemenata[1])


def Test_all_element_numbers(grillage):
    GrillageMesh.element_size_mesh(mesh1, grillage)
    GrillageMesh.transition_element_mesh(mesh1, grillage)
    suma = 0
    for plate in grillage.plating().values():
        n_prijelaznih = mesh1.get_tr_element_num(plate)
        n_svih = mesh1.get_element_number(plate)
        n_prirubnica = mesh1.get_flange_element_num(plate)
        ukupno = n_svih[0] * n_svih[1]
        suma += ukupno
        print("Zona oplate", plate.id, ", broj svih elemenata po x i y:", n_svih[0], n_svih[1], ", broj prijelaznih",
              n_prijelaznih[0], n_prijelaznih[1], ", broj elemenata prirubnica:", n_prirubnica[0], n_prirubnica[1])
    print("Ukupan broj quad elemenata oplate:", suma)


def Test_get_mesh_dim(grillage, plate_id):
    plate = grillage.plating()[plate_id]
    GrillageMesh.element_size_mesh(mesh1, grillage)
    GrillageMesh.transition_element_mesh(mesh1, grillage)
    print("Dimenzije x:", GrillageMesh.get_mesh_dim_x(mesh1, plate))
    print("Dimenzije y:", GrillageMesh.get_mesh_dim_y(mesh1, plate))


def Test_generate_plating_zone_elements(grillage, plate_id, start_node_id, start_element_id):
    plate = grillage.plating()[plate_id]
    GrillageMesh.element_size_mesh(mesh1, grillage)
    GrillageMesh.transition_element_mesh(mesh1, grillage)
    GrillageMesh.identify_unique_property(mesh1, grillage)

    GrillageMesh.generate_plating_zone_elements(mesh1, plate, start_node_id, start_element_id)


def Test_get_web_el_height(psm_id, segment_id):
    segment = hc_variant.longitudinal_members()[psm_id].segments[segment_id - 1]
    dim = GrillageMesh.get_web_el_height(mesh1, segment)
    print(dim)


def Test_generate_segment_web_elements(grillage, direction: BeamDirection, psm_id, segment_id, start_node_id, start_element_id):
    segment = None

    if direction == BeamDirection.LONGITUDINAL:
        segment = grillage.longitudinal_members()[psm_id].segments[segment_id - 1]

    elif direction == BeamDirection.TRANSVERSE:
        segment = grillage.transverse_members()[psm_id].segments[segment_id - 1]

    GrillageMesh.element_size_mesh(mesh1, grillage)
    GrillageMesh.transition_element_mesh(mesh1, grillage)
    GrillageMesh.generate_segment_web_elements(mesh1, segment, start_node_id, start_element_id)


def AOS_Logic_test(AxisOfSymm: AOS):
    mesh_dim_y = np.zeros(7)
    mesh_dim_x = np.zeros(10)

    row_limit = len(mesh_dim_y)
    column_limit = len(mesh_dim_x)

    if AxisOfSymm == AOS.LONGITUDINAL:          # Longitudinal axis of symmetry splits the plating zone
        if np.mod(len(mesh_dim_y), 2) == 0:     # Even number of elements along y axis
            row_limit = len(mesh_dim_y) / 2
        else:                                   # Odd number of elemenets along y axis
            row_limit = np.ceil(len(mesh_dim_y) / 2)
            # Divide coordinate y of last row of nodes by 2

    elif AxisOfSymm == AOS.TRANSVERSE:          # Transverse axis of symmetry splits the plating zone
        if np.mod(len(mesh_dim_y), 2) == 0:     # Even number of elements along x axis
            column_limit = len(mesh_dim_x) / 2
        else:                                   # Odd number of elemenets along x axis
            column_limit = np.ceil(len(mesh_dim_x) / 2)
            # Divide coordinate x of last column of nodes by 2

    elif AxisOfSymm == AOS.BOTH:  # Both longitudinal and transverse axis of symmetry splits the plating zone
        if np.mod(len(mesh_dim_y), 2) == 0:     # Even number of elements along y axis
            row_limit = len(mesh_dim_y) / 2
        else:                                   # Odd number of elemenets along y axis
            row_limit = np.ceil(len(mesh_dim_y) / 2)
            # Divide coordinate y of last row of nodes by 2

        if np.mod(len(mesh_dim_y), 2) == 0:     # Even number of elements along x axis
            column_limit = len(mesh_dim_x) / 2
        else:                                   # Odd number of elemenets along x axis
            column_limit = np.ceil(len(mesh_dim_x) / 2)
            # Divide coordinate x of last column of nodes by 2

    print(row_limit, column_limit)


def Test_generate_mesh(grillage, aos: AOS):
    GrillageMesh.element_size_mesh(mesh1, grillage)
    GrillageMesh.transition_element_mesh(mesh1, grillage)
    GrillageMesh.generate_mesh(mesh1, grillage, aos)


# AOS_Logic_test(AOS.TRANSVERSE)
# Test_get_reduced_plate_dim(hc_variant, 16)
# print(GrillageMesh.find_closest_divisor(4935, 660))
# Test_element_size_stiffener_spacing(hc_variant, 1)    # Dimenzije elementa oplate za odabranu zonu prema razmaku ukrepa
# Test_get_flange_el_length(hc_variant, BeamDisrection.LONGITUDINAL, 1, 1)    # Dimenzije elementa prirubnice za odabrani segment
# Test_element_size_plating_zone(hc_variant, 1)           # Dimenzije elemenata na odabranoj zoni oplate prema razmaku ukrepa i ar prirubnice
# Test_element_aspect_ratio(700, 100)
# Test_refine_plate_element(935, 798, 847.6)
# Test_refine_plate_element(4500, 200, 10000)           # loop break test
# Test_ALL_element_size_plating_zone(hc_variant)        # dim_x i dim_y za sve zone oplate, prikaz u matrici
# Test_element_size_mesh(hc_variant)                    # Konacno odabrane dimenzije mreze po x i y
# Test_get_base_dim_x(hc_variant, 2)                    # Odabrana osnovna x dimenzija za neko polje oplate
# Test_get_base_dim_y(hc_variant, 2)                    # Odabrana osnovna y dimenzija za neko polje oplate
# Test_element_size_transition_x(hc_variant, 1)         # x dimenzije prijelaznih elemenata
# Test_element_size_transition_y(hc_variant, 1)         # y dimenzije prijelaznih elemenata
# Test_transition_element_dimensions(hc_variant, 5)     # Obje dimenzije prijelaznih elemenata za neko polje oplate
# Test_get_flange_el_width(2, 1)
# Test_all_plating_zones_mesh_dimensions(hc_variant)    # Odabrane x i y dimenzije elemenata za sva polja oplate
# Test_identify_unique_property(hc_variant)
# Test_transition_element_mesh(hc_variant)
# Test_get_tr_dim_x(hc_variant, 2)
# Test_get_tr_dim_y(hc_variant, 1)
# Test_get_tr_element_num(hc_variant, 2)
# Test_get_element_number(hc_variant, 2)
# Test_all_element_numbers(hc_variant)
# Test_get_mesh_dim(hc_variant, 1)                        # Dimenzije x i y svih elemenata duž odabrane zone oplate
# Test_get_web_el_height(1, 1)
# Test_generate_plating_zone_elements(hc_variant, 1, 1, 1)
# Test_generate_segment_web_elements(hc_variant, BeamDirection.LONGITUDINAL, 1, 1, 1, 1)
# Test_generate_mesh(hc_variant, AOS.BOTH)

end = timer()

print("***************************************************************************************************************")
print("Code run time =", end - start, "s")
