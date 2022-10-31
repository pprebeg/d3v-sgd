"""
Modul za testiranje diskretizacije učitanog modela

"""

from femdir.grillage_mesh import *
from timeit import default_timer as timer

start = timer()

# Ucitavanje topologije iz datoteke
filename = '../grillage savefiles/hc_var_1_savefile.gin'
hc_variant = GrillageModelData(filename).read_file()

# Izrada MeshSize objekta
mesh1 = MeshV1(hc_variant, AOS.LONGITUDINAL)    # Unos geometrije i odabir globalne osi simetrije konstrukcije

# Kontrola mreže
mesh1.min_num_ebs = 1                   # Postavljanje minimalnog broja elemenata između ukrepa
mesh1.min_num_eweb = 3                  # Postavljanje minimalnog broja elemenata duž visine struka
mesh1.num_eaf = 1                       # Postavljanje broja elemenata u smjeru širine prirubnice
mesh1.flange_aspect_ratio = 7           # Postavljanje aspektnog odnosa elemenata prirubnica jakih nosača i oplate uz struk jakih nosača
mesh1.plate_aspect_ratio = 4            # Postavljanje aspektnog odnosa elemenata oplate i strukova jakih nosača
mesh1.des_plate_aspect_ratio = 3        # Postavljanje poželjnog aspektnog odnosa elemenata oplate

mesh1.calculate_mesh_dimensions()       # Izračun svih dimenzija za odabranu mrežu


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
    plate_dim = mesh1.get_reduced_plate_dim(plate)
    paralelna_dimenzija = mesh1.element_size_para_to_stiffeners(plate, plate_dim)
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
    plate_dim = mesh1.get_reduced_plate_dim(plate)
    dims = mesh1.element_size_plating_zone(plate, plate_dim)
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
            plate_dim = mesh1.get_reduced_plate_dim(plate)
            polje1[stupac, redak] = mesh1.element_size_plating_zone(plate, plate_dim)[0]
            plate_id += 1
    print("Sve x dimenzije elemenata: \n", polje1, "\n")

    # dimenzije karakteristicnih elemenata na cijeloj oplati po y osi
    plate_id = 1
    for stupac in range(0, n_y):
        for redak in range(0, n_x):
            plate = hc_variant.plating()[plate_id]
            plate_dim = mesh1.get_reduced_plate_dim(plate)
            polje1[stupac, redak] = mesh1.element_size_plating_zone(plate, plate_dim)[1]
            plate_id += 1
    print("Sve y dimenzije elemenata: \n", polje1)


def Test_element_size_mesh():
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
    for plate in mesh1.plating_zones.values():
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
    # plate = mesh1.plating_zones[plate_id]
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
    for plate in mesh1.plating_zones.values():
        tr_el_longitudinal, tr_el_transverse = mesh1.get_tr_element_num(plate)
        n_el_x, n_el_y = mesh1.get_element_number(plate)
        flange_element_number_x, flange_element_number_y = mesh1.get_flange_element_num(plate)
        ukupno = n_el_x * n_el_y
        suma += ukupno
        print("Zona oplate", plate.id, ", broj svih elemenata po x i y:", n_el_x, n_el_y, ", broj prijelaznih",
              tr_el_longitudinal, tr_el_transverse, ", broj elemenata prirubnica:", flange_element_number_x, flange_element_number_y)
    print("Ukupan broj quad elemenata oplate:", suma)


def Test_get_mesh_dim(plate_id):
    plate = hc_variant.plating()[plate_id]
    print("Dimenzije svih konačnih elemenata redom za zonu oplate", plate_id)
    print("Dimenzije x:", mesh1.get_mesh_dim_x(plate))
    print("Dimenzije y:", mesh1.get_mesh_dim_y(plate))


def Test_get_all_mesh_dim():
    for plate in mesh1.plating_zones.values():
        print("Dimenzije svih konačnih elemenata redom za zonu oplate", plate.id)
        print("Dimenzije x:", mesh1.get_mesh_dim_x(plate))
        print("Dimenzije y:", mesh1.get_mesh_dim_y(plate))
        print("\n")


def Test_get_web_el_height(psm_id, segment_id):
    segment = hc_variant.longitudinal_members()[psm_id].segments[segment_id - 1]
    dim = mesh1.get_web_el_height(segment)
    print(dim)


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
    transition_dims = mesh1.transition_element_size_plating_zone(plate, segment_id)
    print("Dimenzije prijelaznog elemenata na zoni oplate", plate_id, "uz segment", segment_id, ":", transition_dims)


def Test_SegmentMesh(psm_id, segment_id):
    segment = hc_variant.longitudinal_members()[psm_id].segments[segment_id - 1]
    print(segment.id)
    # segmentmesh1 = SegmentMesh(segment, 1, 1)


def Test_psm_extent():
    longs = mesh1.longitudinal_psm_extent()
    trans = mesh1.transverse_psm_extent()
    print("Prema unesenoj osi simetrije", mesh1.axis_of_symm, ", od ukupno",
          len(hc_variant.longitudinal_members()), "jakih uzdužnih nosača na modelu, radi se mreža za njih", len(longs))
    print("Prema unesenoj osi simetrije", mesh1.axis_of_symm, ", od ukupno",
          len(hc_variant.transverse_members()), "jakih poprečnih nosača na modelu, radi se mreža za njih", len(trans))


def Test_segment_extent():
    mesh1.grillage_segment_extent()


def Test_grillage_plate_extent():
    mesh1.grillage_plate_extent()
    print("Odabrana je", mesh1.axis_of_symm, "os simetrije.")

    print("Ukupno identificiranih zona za izradu mreže:", len(mesh1.plating_zones))
    print("ID zona oplate za izradu mreže:")
    for key in mesh1.plating_zones:
        print("Ključ:", key, ", ID zone oplate", mesh1.plating_zones[key].id)

    print("Od tih upisanih zona, dijele se na različite tipove:")

    print("Zone oplate za izradu pune mreže:")
    for key in mesh1.full_plate_zones:
        print("Ključ:", key, ", ID zone oplate", mesh1.plating_zones[key].id)

    print("Zone oplate za izradu polovične mreže, presječenih uzdužnom osi simetrije:")
    for key in mesh1.long_half_plate_zones:
        print("Ključ:", key, ", ID zone oplate", mesh1.plating_zones[key].id)

    print("ID zona oplate za izradu polovične mreže, presječenih poprečnom osi simetrije:")
    for key in mesh1.tran_half_plate_zones:
        print("Ključ:", key, ", ID zone oplate", mesh1.plating_zones[key].id)

    print("Zone oplate za izradu četvrtinske mreže:")
    for key in mesh1.quarter_plate_zone:
        print("Ključ:", key, ", ID zone oplate", mesh1.plating_zones[key].id)


def Test_PlatingZoneMesh(plate_id, split_along=AOS.NONE):
    plate = hc_variant.plating()[plate_id]

    mesh1.grillage_plate_extent()        # Izračun koje zone oplate se meshiraju
    PlatingZoneMesh(mesh1, plate, 1, 1, split_along).generate_mesh()     # izrada mreže jedne zone oplate


def Test_plating_zones_ref_array():
    mesh1.grillage_plate_extent()
    arr = mesh1.plating_zones_ref_array

    n_redaka, n_stupaca = np.shape(arr)
    print(arr)
    print("Postoji", n_redaka, "redaka i", n_stupaca, "stupaca zona oplate koje se meshiraju")


def Test_get_plate_dim(plate_id):
    mesh1.grillage_plate_extent()
    plate = hc_variant.plating()[plate_id]
    # full_dim = plate.plate_dim_parallel_to_stiffeners() * 1000
    full_dim = mesh1.get_reduced_plate_dim(plate)
    print("Puna dimenzija:", full_dim)
    print(mesh1.get_plate_dim(plate, full_dim))


def Test_calc_element_base_size():
    print(mesh1.calc_element_base_size())


def Test_calculate_mesh_dimensions():
    x_dimenzije = mesh1.plate_edge_node_x
    y_dimenzije = mesh1.plate_edge_node_y

    print("Razmaci između čvorova (dimenzije elemenata) duž x osi:")
    for i in x_dimenzije.keys():
        print("  Zona oplate", i, ":", x_dimenzije[i])

    print("Razmaci između čvorova (dimenzije elemenata) duž y osi:")
    for i in y_dimenzije.keys():
        print("  Zona oplate", i, ":", y_dimenzije[i])


def Test_Segment_element_generation(direction: BeamDirection, psm_id, segment_id):
    segment = None
    if direction == BeamDirection.LONGITUDINAL:
        segment = hc_variant.longitudinal_members()[psm_id].segments[segment_id - 1]
    elif direction == BeamDirection.TRANSVERSE:
        segment = hc_variant.transverse_members()[psm_id].segments[segment_id - 1]

    start_node_id = 1
    start_element_id = 1
    seg_mesh = SegmentV1(mesh1, segment, start_node_id, start_element_id)
    print(seg_mesh.generate_mesh())


def Test_aos_stiffener(plate_id):
    plate = hc_variant.plating()[plate_id]
    test_btw = mesh1.aos_between_stiffeners(plate)
    test_on = mesh1.aos_on_stiffener(plate)
    print("Test osi simetrije između ukrepa:", test_btw, ", test osi simetrije na ukrepi:", test_on)


def Test_aos_on_segment(direction: BeamDirection, psm_id, segment_id):
    segment = None
    if direction == BeamDirection.LONGITUDINAL:
        segment = hc_variant.longitudinal_members()[psm_id].segments[segment_id - 1]
    elif direction == BeamDirection.TRANSVERSE:
        segment = hc_variant.transverse_members()[psm_id].segments[segment_id - 1]
    test = mesh1.aos_on_segment(segment)
    print("Test prolazi li os simetrije uz segment:", test)


def Test_get_segment_web_element_property(direction: BeamDirection, psm_id, segment_id):
    segment = None
    if direction == BeamDirection.LONGITUDINAL:
        segment = hc_variant.longitudinal_members()[psm_id].segments[segment_id - 1]
        print("Svojstva segmenta ID", segment_id, "jakog uzdužnog nosača", psm_id, ":")

    elif direction == BeamDirection.TRANSVERSE:
        segment = hc_variant.transverse_members()[psm_id].segments[segment_id - 1]
        print("Svojstva segmenta ID", segment_id, "jakog poprečnog nosača", psm_id, ":")

    seg_mesh = SegmentV1(mesh1, segment, 1, 1)
    seg_mesh.get_web_element_property()


def Test_Model_check():
    check1 = ModelCheck(hc_variant)
    hc_variant.assign_symmetric_members()

    long_psm_symm = check1.longitudinal_psm_symmetry()
    tran_psm_symm = check1.transverse_psm_symmetry()
    central_long = check1.central_longitudinal()
    central_tran = check1.central_transversal()
    long_plate = check1.longitudinal_plate_symmetry()
    tran_plate = check1.transverse_plate_symmetry()
    long_symm = check1.longitudinal_symmetry_tests()
    tran_symm = check1.transverse_symmetry_tests()
    long_segment = check1.longitudinal_segment_symmetry()
    tran_segment = check1.transverse_segment_symmetry()
    aos = check1.assign_symmetry()

    print("Uzdužna simetrija položaja nosača:", long_psm_symm, ", poprečna simetrija položaja nosača:", tran_psm_symm)
    print("Centralni uzdužni nosači:",  central_long, ", centralni poprečni nosači", central_tran)
    print("Uzdužna simetrija zona oplate:", long_plate, ", poprečna simetrija zona oplate:", tran_plate)
    print("Uzdužna simetrija segmenata:", long_segment, ", poprečna simetrija segmenata:", tran_segment)

    print("Konačna uzdužna simetrija:", long_symm, ", konačna poprečna simetrija:", tran_symm)
    print("Konačno odabrana simetrija Axis of Symmetry =", aos)


def Test_mesh_feasibility():
    check1 = ModelCheck(hc_variant)
    hc_variant.assign_symmetric_members()
    hc_var_check = check1.mesh_feasibility()
    print(hc_var_check)


# Test_get_reduced_plate_dim(16)
# Test_find_closest_divisor(4133, 935)
# Test_find_closest_divisor(4238, 935)
# Test_element_size_para_to_stiffeners(1)    # Dimenzije elementa oplate za odabranu zonu prema razmaku ukrepa
# Test_get_flange_el_length(BeamDirection.LONGITUDINAL, 1, 1)    # Dimenzije elementa prirubnice za odabrani segment
# Test_element_size_plating_zone(1)         # Dimenzije elemenata na odabranoj zoni oplate prema razmaku ukrepa i ar prirubnice

# Test_element_aspect_ratio(700, 100)
# Test_refine_plate_element(4133, 830)
# Test_ALL_element_size_plating_zone()        # dim_x i dim_y za sve zone oplate, prikaz u matrici
# Test_element_size_mesh()                    # Konacno odabrane dimenzije mreze po x i y
# Test_all_plating_zones_mesh_dimensions()    # Odabrane x i y dimenzije elemenata za sva polja oplate

# Test_get_flange_el_width(1, 1)
# Test_identify_unique_property()
# Test_get_tr_element_num(1)
# Test_get_element_number(2)
# Test_all_element_numbers()
# Test_get_mesh_dim(2)                        # Dimenzije x i y svih elemenata duž odabrane zone oplate
# Test_get_all_mesh_dim()                     # Dimenzije x i y svih elemenata duž svih generiranih zona oplate

# Test_get_web_el_height(1, 1)
# Test_get_min_flange_el_length()
# Test_get_min_flange_el_length_between_psm(1, 2)
# Test_find_largest_divisor(4500, 1000)
# Test_element_size_transition(1, 2)
# Test_get_tr_dim_x(5)
# Test_get_tr_dim_y(4)

# Test_SegmentMesh(2, 1)
# Test_psm_extent()
# Test_segment_extent()
# Test_grillage_plate_extent()
# Test_plating_zones_ref_array()
# Test_get_plate_dim(2)
# Test_calc_element_base_size()
# Test_calculate_mesh_dimensions()
# Test_PlatingZoneMesh(1, AOS.NONE)                                         # Izrada mreže jedne zone oplate
# PlateMesh(mesh1).generate_mesh()                                          # Izrada mreže cijele oplate
# Test_Segment_element_generation(BeamDirection.LONGITUDINAL, 2, 1)
# Test_aos_stiffener(4)
# Test_aos_on_segment(BeamDirection.TRANSVERSE, 3, 1)
# Test_get_segment_web_element_property(BeamDirection.LONGITUDINAL, 3, 4)
# Test_Model_check()
# print(ModelCheck(hc_variant).assign_symmetry())   # Konačno odabrana os simetrije s kojom se ide u izradu mreže
# Test_mesh_feasibility()


end = timer()

print("***************************************************************************************************************")
print("Code run time =", end - start, "s")
