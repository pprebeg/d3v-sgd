"""
Modul za testiranje diskretizacije učitanog modela

"""

from femdir.grillage_mesh import *
from timeit import default_timer as timer

start = timer()

# Ucitavanje topologije iz datoteke
filename = '../grillage savefiles/hc_var_1_savefile.gin'
hc_variant = GrillageModelData(filename).read_file()


# Generacija mreže
mesh1 = MeshV1(hc_variant, AOS.NONE)   # Izrada MeshSize objekta

# Kontrola mreže
mesh1.min_num_ebs = 1               # Postavljanje minimalnog broja elemenata između ukrepa
# mesh1.min_num_eweb = 3              # Postavljanje minimalnog broja elemenata duž visine struka
mesh1.flange_aspect_ratio = 7       # Postavljanje aspektnog odnosa elemenata prirubnica jakih nosača i oplate uz struk jakih nosača
mesh1.plate_aspect_ratio = 4        # Postavljanje aspektnog odnosa elemenata oplate i strukova jakih nosača
mesh1.des_plate_aspect_ratio = 3    # Postavljanje poželjnog aspektnog odnosa elemenata oplate

# Metode za izračun dimenzija pojedine mreže:
mesh1.element_base_size_mesh()      # Izračun osnovnih dimenzija mreže
mesh1.transition_element_mesh()     # Izračun dimenzija prijelaznih elemenata

"""
# mesh1.GenerateNodes(hc_variant)               # Generacija cvorova
# mesh1.GenerateElements(hc_variant)            # Generacija elemenata
# Zapis mreze u datoteku
# MeshData("mesh1_savefile.txt").write_file(mesh1)
"""


def Test_get_reduced_plate_dim(plate_id):
    plate = hc_variant.plating()[plate_id]
    reduced_dim = mesh1.get_reduced_plate_dim(plate)
    print("Reducirana dimenzija paralalna s ukrepama na zoni oplate", plate_id, "iznosi:", reduced_dim)


def Test_find_closest_divisor(length, spacing):
    div = mesh1.find_closest_divisor(length, spacing)
    if div is None:
        print("Nije pronađen djelitelj!")
    else:
        print("Pronađen je djelitelj! Broj", length, "se može podijeliti na", div, "dijelova, tako da svaki iznosi", length/div,
              ", što bi trebalo biti blizu zadanog broja", spacing)


def Test_element_size_para_to_stiffeners(plate_id):
    plate = hc_variant.plating()[plate_id]
    okomita_dimenzija = mesh1.element_size_perp_to_stiffeners(plate)
    paralelna_dimenzija = mesh1.element_size_para_to_stiffeners(plate)
    n_elem = mesh1.min_num_ebs
    ar = MeshSize.element_aspect_ratio(okomita_dimenzija, paralelna_dimenzija)
    print("Dimenzije quad elementa oplate uz", n_elem, "element između ukrepa, okomito:", okomita_dimenzija,
          "mm, paralelno:", paralelna_dimenzija, "mm", ", aspektni odnos:", ar)


def Test_get_flange_el_length(direction: BeamDirection, psm_id, segment_id):
    segment = None

    if direction == BeamDirection.LONGITUDINAL:
        segment = hc_variant.longitudinal_members()[psm_id].segments[segment_id]

    elif direction == BeamDirection.TRANSVERSE:
        segment = hc_variant.transverse_members()[psm_id].segments[segment_id]

    dim = mesh1.get_flange_el_length(segment)
    print("Maksimalna dimenzija elementa prirubnice prema aspektnom odnosu:", dim, "mm")


def Test_element_size_plating_zone(plate_id):
    plate = hc_variant.plating()[plate_id]
    dims = mesh1.element_size_plating_zone(plate)
    print("Osnovne dimenzije elementa na zoni oplate", plate_id)
    dim_x = dims[0]
    dim_y = dims[1]
    print("dim_x =", dim_x, ", dim_y =", dim_y)


def Test_element_aspect_ratio(dim_x, dim_y):
    ar = mesh1.element_aspect_ratio(dim_x, dim_y)
    print("Aspektni odnos =", ar)


def Test_refine_plate_element(length, dim_limit):
    dim = mesh1.refine_plate_element(length, dim_limit)
    print("Duljinu", length, "je potrebno podijeliti na jednake dijelove, tako da dimenzija elementa ne prelazi", dim_limit,
          ". \n Odabrana je dimenzija", dim)


def Test_ALL_element_size_plating_zone():
    n_x = int(hc_variant.N_transverse - 1)    # Broj polja u uzduznom smjeru
    n_y = int(hc_variant.N_longitudinal - 1)  # Broj polja u poprecnom smjeru

    polje1 = np.zeros((n_y, n_x))
    # dimenzije karakteristicnih elemenata na cijeloj oplati po x osi
    plate_id = 1
    for stupac in range(0, n_y):
        for redak in range(0, n_x):
            plate = hc_variant.plating()[plate_id]
            polje1[stupac, redak] = mesh1.element_size_plating_zone(plate)[0]
            plate_id += 1
    print("Sve x dimenzije elemenata: \n", polje1, "\n")

    # dimenzije karakteristicnih elemenata na cijeloj oplati po y osi
    plate_id = 1
    for stupac in range(0, n_y):
        for redak in range(0, n_x):
            plate = hc_variant.plating()[plate_id]
            polje1[stupac, redak] = mesh1.element_size_plating_zone(plate)[1]
            plate_id += 1
    print("Sve y dimenzije elemenata: \n", polje1)


def Test_element_size_mesh():
    mesh1.element_base_size_mesh()
    print("Konačno odabrane dimenzije mreže po x:", mesh1.mesh_dim_x)
    print("Konačno odabrane dimenzije mreže po y:", mesh1.mesh_dim_y)


def Test_get_flange_el_width(psm_id, segment_id):
    segment = hc_variant.longitudinal_members()[psm_id].segments[segment_id - 1]

    print("Jaki uzduzni nosac", psm_id, ", ID segmenta:", segment.id,
          " ,BeamProperty ID:", segment.beam_prop.id,
          ", tip profila: ", segment.beam_prop.beam_type,
          ", tf =", segment.beam_prop.tf,
          ", dimenzija elementa prirubnice po širini:", mesh1.get_flange_el_width(segment), "mm")


def Test_all_plating_zones_mesh_dimensions():
    for plate in hc_variant.plating().values():
        dim_x = mesh1.get_base_dim_x(plate)
        dim_y = mesh1.get_base_dim_y(plate)
        print("Zona oplate ID:", plate.id, ",   dim_x =", "{:.2f}".format(dim_x), "mm", ",   dim_y =", "{:.2f}".format(dim_y), "mm")


def Test_identify_unique_property():
    mesh1.identify_unique_property()

    for prop in mesh1.unique_properties.values():
        plate_prop = len(prop.plate_prop)
        beam_prop = len(prop.beam_prop)
        print("Unique property ID:", prop.id, ", tp =", prop.tp, "mm, material ID", prop.mat.id,
              ", upisano istih na modelu, plate:", plate_prop, ", beam:", beam_prop)


def Test_get_tr_dim_x(plate_id):
    plate = hc_variant.plating()[plate_id]
    dims = mesh1.get_tr_dim_x(plate)
    print("Dimenzije x prijelaznih elemenata na zoni oplate", plate_id, ":", dims)


def Test_get_tr_dim_y(plate_id):
    plate = hc_variant.plating()[plate_id]
    dims = mesh1.get_tr_dim_y(plate)
    print("Dimenzije y prijelaznih elemenata na zoni oplate", plate_id, ":", dims)


def Test_get_tr_element_num(plate_id):
    plate = hc_variant.plating()[plate_id]
    n_elemenata = mesh1.get_tr_element_num(plate)
    poprecni_segment = n_elemenata[0]
    uzduzni_segment = n_elemenata[1]
    print("Broj prijelaznih elemenata na zoni oplate", plate_id, "uz elemente prirubnice uzdužnog segmenta duž osi y:", uzduzni_segment,
          ", uz elemente prirubnice poprečnog segmenta duž osi x:", poprecni_segment)


def Test_get_element_number(plate_id):
    plate = hc_variant.plating()[plate_id]
    n_elemenata = mesh1.get_element_number(plate)
    print("Broj elemenata na zoni oplate", plate_id, "po x:", n_elemenata[0], ", po y:", n_elemenata[1])


def Test_all_element_numbers():
    suma = 0
    for plate in hc_variant.plating().values():
        n_prijelaznih = mesh1.get_tr_element_num(plate)
        n_svih = mesh1.get_element_number(plate)
        n_prirubnica = mesh1.get_flange_element_num(plate)
        ukupno = n_svih[0] * n_svih[1]
        suma += ukupno
        print("Zona oplate", plate.id, ", broj svih elemenata po x i y:", n_svih[0], n_svih[1], ", broj prijelaznih",
              n_prijelaznih[0], n_prijelaznih[1], ", broj elemenata prirubnica:", n_prirubnica[0], n_prirubnica[1])
    print("Ukupan broj quad elemenata oplate:", suma)


def Test_get_mesh_dim(plate_id):
    plate = hc_variant.plating()[plate_id]
    print("Dimenzije svih konačnih elemenata redom za zonu oplate", plate_id)
    print("Dimenzije x:", mesh1.get_mesh_dim_x(plate))
    print("Dimenzije y:", mesh1.get_mesh_dim_y(plate))


def Test_generate_plating_zone_elements(plate_id, start_node_id, start_element_id):
    plate = hc_variant.plating()[plate_id]
    mesh1.element_base_size_mesh()
    mesh1.transition_element_mesh()
    mesh1.identify_unique_property()

    mesh1.generate_plating_zone_elements(plate, start_node_id, start_element_id)


def Test_get_web_el_height(psm_id, segment_id):
    segment = hc_variant.longitudinal_members()[psm_id].segments[segment_id - 1]
    dim = mesh1.get_web_el_height(segment)
    print(dim)


def Test_generate_segment_web_elements(direction: BeamDirection, psm_id, segment_id, start_node_id, start_element_id):
    segment = None

    if direction == BeamDirection.LONGITUDINAL:
        segment = hc_variant.longitudinal_members()[psm_id].segments[segment_id - 1]

    elif direction == BeamDirection.TRANSVERSE:
        segment = hc_variant.transverse_members()[psm_id].segments[segment_id - 1]

    mesh1.generate_segment_web_elements(segment, start_node_id, start_element_id)


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


def Test_generate_mesh():
    mesh1.generate_mesh()


def Test_get_min_flange_el_length():
    psm1 = PrimarySuppMem(1, BeamDirection.LONGITUDINAL, 0.1, hc_variant)
    ST24 = MaterialProperty(1, 210000, 0.3, 7850, 235, "ST24")
    beam_prop1 = FBBeamProperty(1, 1089, 10, ST24)
    # beam_prop2 = FBBeamProperty(2, 1089, 10, ST24)
    beam_prop2 = TBeamProperty(2, 1089, 10, 545, 40, ST24)
    segment1 = Segment(1, beam_prop1, psm1, psm1, psm1)
    segment2 = Segment(2, beam_prop2, psm1, psm1, psm1)

    min_dim = mesh1.get_min_flange_el_length(segment1, segment2)
    print("Minimalna vrijednost je", min_dim)


def Test_get_min_flange_el_length_between_psm(member1_id, member2_id):
    member1 = hc_variant.longitudinal_members()[member1_id]
    member2 = hc_variant.longitudinal_members()[member2_id]
    min_dim = mesh1.get_min_flange_el_length_between_psm(member1, member2)
    print("Najmanja vrijednost maksimalne duljine prirubnice svih segmenata između jakih nosača:", min_dim)


def Test_find_largest_divisor(length, max_val):
    divisor = mesh1.find_largest_divisor(length, max_val)
    print("Najveći djelitelj broja", length, ", koji rezultira dimenzijom manjom ili jednakom", max_val, "je", divisor,
          ", što daje vrijednost", length / divisor)


def Test_element_size_transition(plate_id, segment_id):
    plate = hc_variant.plating()[plate_id]
    transition_dims = mesh1.element_size_transition(plate, segment_id)
    print("Dimenzije prijelaznog elemenata na zoni oplate", plate_id, "uz segment", segment_id, ":", transition_dims)


def Test_SegmentMesh(psm_id, segment_id):
    segment = hc_variant.longitudinal_members()[psm_id].segments[segment_id - 1]
    # segmentmesh1 = SegmentMesh(segment, 1, 1)


def Test_grillage_plate_extent():
    mesh1.grillage_plate_extent()
    print("Odabrana je", mesh1.axis_of_symm, "os simetrije.")
    print("ID zona oplate za izradu pune mreže:")
    for plate in mesh1.full_plate_zones.values():
        print(plate.id)
    print("ID zona oplate za izradu polovične mreže, presječenih uzdužnom osi simetrije:")
    for plate in mesh1.long_half_plate_zones.values():
        print(plate.id)
    print("ID zona oplate za izradu polovične mreže, presječenih poprečnom osi simetrije:")
    for plate in mesh1.tran_half_plate_zones.values():
        print(plate.id)
    print("ID zona oplate za izradu četvrtinske mreže:")
    for plate in mesh1.quarter_plate_zone.values():
        print(plate.id)


# AOS_Logic_test(AOS.TRANSVERSE)
# Test_get_reduced_plate_dim(16)
# Test_find_closest_divisor(4133, 935)
# Test_find_closest_divisor(4238, 935)
# Test_element_size_para_to_stiffeners(5)    # Dimenzije elementa oplate za odabranu zonu prema razmaku ukrepa
# Test_get_flange_el_length(BeamDirection.LONGITUDINAL, 1, 1)    # Dimenzije elementa prirubnice za odabrani segment
# Test_element_size_plating_zone(1)         # Dimenzije elemenata na odabranoj zoni oplate prema razmaku ukrepa i ar prirubnice
# Test_element_aspect_ratio(700, 100)
# Test_refine_plate_element(4133, 830)
# Test_ALL_element_size_plating_zone()        # dim_x i dim_y za sve zone oplate, prikaz u matrici
# Test_element_size_mesh()                    # Konacno odabrane dimenzije mreze po x i y
# Test_get_flange_el_width(1, 1)
# Test_all_plating_zones_mesh_dimensions()    # Odabrane x i y dimenzije elemenata za sva polja oplate
# Test_identify_unique_property()
# Test_get_tr_dim_x(1)
# Test_get_tr_dim_y(1)
# Test_get_tr_element_num(1)
# Test_get_element_number(2)
# Test_all_element_numbers()
# Test_get_mesh_dim(1)                        # Dimenzije x i y svih elemenata duž odabrane zone oplate
# Test_get_web_el_height(1, 1)
# Test_generate_plating_zone_elements(1, 1, 1)
# Test_generate_segment_web_elements(BeamDirection.LONGITUDINAL, 1, 1, 1, 1)
# Test_generate_mesh()
# Test_get_min_flange_el_length()
# Test_get_min_flange_el_length_between_psm(1, 2)
# Test_find_largest_divisor(4500, 1000)
# Test_element_size_transition(1, 1)
# Test_SegmentMesh(2, 1)
# Test_grillage_plate_extent()

end = timer()

print("***************************************************************************************************************")
print("Code run time =", end - start, "s")
