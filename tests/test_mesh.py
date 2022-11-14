"""
Modul za testiranje diskretizacije učitanog modela
"""
from grillage.grillage_mesher import *
from timeit import default_timer as timer
import femdir.geofementity as gfe
import os

# Get the current working directory
cwd = os.getcwd()
dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path)
# Print the current working directory
cwd = os.getcwd()
#print("Current working directory: {0}".format(cwd))
start = timer()
hc_var = 5
filename = str("../grillage savefiles/hc_var_") + str(hc_var) + str("_savefile.gin")
hc_variant = GrillageModelData(filename).read_file()        # Učitavanje topologije iz datoteke
print("Generating FE mesh on grillage variant", hc_var)





# Potrebno računati dimenzije mreže za sve testove osim generate_grillage_mesh_v1()
# test_mesh.calculate_mesh_dimensions()       # Izračun svih dimenzija za odabranu mrežu



def test_mesh():
    # extents = MeshExtent(hc_variant, AOS.NONE)                # Opseg izrade mreže uz ručni odabir simetrije
    extents = MeshExtent(hc_variant)  # Opseg izrade mreže uz automatsko prepoznavanje simetrije
    test_mesh = ElementSizeV1(extents)  # Izračun dimenzija mreže za V1
    # test_mesh = ElementSizeV2(extents)       # Izračun dimenzija mreže za V2
    # Kontrola mreže
    test_mesh.min_num_ebs = 1  # Postavljanje minimalnog broja elemenata između ukrepa
    test_mesh.min_num_eweb = 3  # Postavljanje minimalnog broja elemenata duž visine struka
    test_mesh.num_eaf = 1  # Postavljanje broja elemenata u smjeru širine prirubnice
    test_mesh.flange_aspect_ratio = 8  # Postavljanje aspektnog odnosa elemenata prirubnica jakih nosača i oplate uz struk jakih nosača
    test_mesh.plate_aspect_ratio = 4  # Postavljanje aspektnog odnosa elemenata oplate i strukova jakih nosača
    test_mesh.des_plate_aspect_ratio = 3  # Postavljanje poželjnog aspektnog odnosa elemenata oplate
    gm = GrillageMesh(test_mesh)
    grillage_test_mesh = gm.generate_grillage_mesh_v1("test_mesh")
    grillage_test_mesh.merge_coincident_nodes()
    # grillage_test_mesh.full_model_node_overlap_check()
    return grillage_test_mesh

def Test_get_reduced_plate_dim(plate_id):
    plate = hc_variant.plating()[plate_id]
    reduced_dim = test_mesh.get_reduced_plate_dim(plate)
    print("Reducirana dimenzija paralalna s ukrepama na zoni oplate", plate_id, "iznosi:", reduced_dim)


def Test_find_closest_divisor(length, spacing):
    div = test_mesh.find_closest_divisor(length, spacing)
    if div is None:
        print("Nije pronađen djelitelj!")
    else:
        print("Pronađen je djelitelj! Broj", length, "se može podijeliti na", div, "dijelova, tako da svaki iznosi", length/div,
              ", što bi trebalo biti blizu zadanog broja", spacing)


def Test_element_size_para_to_stiffeners(plate_id):
    plate = hc_variant.plating()[plate_id]
    okomita_dimenzija = test_mesh.element_size_perp_to_stiffeners(plate)
    plate_dim = test_mesh.get_reduced_plate_dim(plate)
    paralelna_dimenzija = test_mesh.element_size_para_to_stiffeners(plate, plate_dim)
    n_elem = test_mesh.min_num_ebs
    ar = MeshSize.element_aspect_ratio(okomita_dimenzija, paralelna_dimenzija)
    print("Dimenzije quad elementa oplate uz", n_elem, "element između ukrepa, okomito:", okomita_dimenzija,
          "mm, paralelno:", paralelna_dimenzija, "mm", ", aspektni odnos:", ar)


def Test_get_flange_el_length(direction: BeamDirection, psm_id, segment_id):
    segment = None

    if direction == BeamDirection.LONGITUDINAL:
        segment = hc_variant.longitudinal_members()[psm_id].segments[segment_id]

    elif direction == BeamDirection.TRANSVERSE:
        segment = hc_variant.transverse_members()[psm_id].segments[segment_id]

    dim = test_mesh.get_flange_el_length(segment)
    print("Maksimalna dimenzija elementa prirubnice prema aspektnom odnosu:", dim, "mm")


def Test_element_size_plating_zone(plate_id):
    plate = hc_variant.plating()[plate_id]
    plate_dim = test_mesh.get_reduced_plate_dim(plate)
    dims = test_mesh.element_size_plating_zone(plate, plate_dim)
    print("Osnovne dimenzije elementa na zoni oplate", plate_id)
    dim_x = dims[0]
    dim_y = dims[1]
    print("dim_x =", dim_x, ", dim_y =", dim_y)


def Test_element_aspect_ratio(dim_x, dim_y):
    ar = test_mesh.element_aspect_ratio(dim_x, dim_y)
    print("Aspektni odnos =", ar)


def Test_refine_plate_element(length, dim_limit):
    dim = test_mesh.refine_plate_element(length, dim_limit)
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
            plate_dim = test_mesh.get_reduced_plate_dim(plate)
            polje1[stupac, redak] = test_mesh.element_size_plating_zone(plate, plate_dim)[0]
            plate_id += 1
    print("Sve x dimenzije elemenata: \n", polje1, "\n")

    # dimenzije karakteristicnih elemenata na cijeloj oplati po y osi
    plate_id = 1
    for stupac in range(0, n_y):
        for redak in range(0, n_x):
            plate = hc_variant.plating()[plate_id]
            plate_dim = test_mesh.get_reduced_plate_dim(plate)
            polje1[stupac, redak] = test_mesh.element_size_plating_zone(plate, plate_dim)[1]
            plate_id += 1
    print("Sve y dimenzije elemenata: \n", polje1)


def Test_element_size_mesh():
    print("Konačno odabrane dimenzije mreže po x:", test_mesh.mesh_dim_x)
    print("Konačno odabrane dimenzije mreže po y:", test_mesh.mesh_dim_y)


def Test_get_flange_el_width(psm_id, segment_id):
    segment = hc_variant.longitudinal_members()[psm_id].segments[segment_id - 1]

    print("Jaki uzduzni nosac", psm_id, ", ID segmenta:", segment.id,
          " ,BeamProperty ID:", segment.beam_prop.id,
          ", tip profila: ", segment.beam_prop.beam_type,
          ", tf =", segment.beam_prop.tf,
          ", dimenzija elementa prirubnice po širini:", test_mesh.get_flange_el_width(segment), "mm")


def Test_all_plating_zones_mesh_dimensions():
    extents.grillage_mesh_extent()
    test_mesh.calc_element_base_size_mesh()
    for plate in extents.all_plating_zones.values():
        dim_x = test_mesh.get_base_dim_x(plate)
        dim_y = test_mesh.get_base_dim_y(plate)
        print("Zona oplate ID:", plate.id, ",   dim_x =", "{:.2f}".format(dim_x), "mm", ",   dim_y =", "{:.2f}".format(dim_y), "mm")


def Test_get_tr_dim_x(plate_id):
    # plate = mesh1.plating_zones[plate_id]
    plate = hc_variant.plating()[plate_id]
    dims = test_mesh.get_tr_dim_x(plate)
    print("Dimenzije x prijelaznih elemenata na zoni oplate", plate_id, ":", dims)


def Test_get_tr_dim_y(plate_id):
    plate = hc_variant.plating()[plate_id]
    dims = test_mesh.get_tr_dim_y(plate)
    print("Dimenzije y prijelaznih elemenata na zoni oplate", plate_id, ":", dims)


def Test_get_tr_element_num(plate_id):
    plate = hc_variant.plating()[plate_id]
    poprecni_segment = test_mesh.get_long_tr_element_num(plate, plate.trans_seg1)
    uzduzni_segment = test_mesh.get_long_tr_element_num(plate, plate.trans_seg2)
    print("Broj prijelaznih elemenata na zoni oplate", plate_id, "uz elemente prirubnice uzdužnog segmenta duž osi y:", uzduzni_segment,
          ", uz elemente prirubnice poprečnog segmenta duž osi x:", poprecni_segment)


def Test_get_element_number(plate_id):
    plate = hc_variant.plating()[plate_id]
    n_elemenata = test_mesh.get_base_element_number(plate)
    print("Broj elemenata osnovnih dimenzija:", plate_id, "po x:", n_elemenata[0], ", po y:", n_elemenata[1])


def Test_get_mesh_dim(plate_id):
    print("Dimenzije svih konačnih elemenata redom za zonu oplate", plate_id)
    print("Dimenzije x:", test_mesh.plate_edge_node_x[plate_id])
    print("Dimenzije y:", test_mesh.plate_edge_node_y[plate_id])


def Test_get_all_mesh_dim():
    for plate in extents.all_plating_zones.values():
        print("Dimenzije svih konačnih elemenata redom za zonu oplate", plate.id)
        print("Dimenzije x:", test_mesh.plate_edge_node_spacing_x(plate))
        print("Dimenzije y:", test_mesh.plate_edge_node_spacing_y(plate))
        print("\n")


def Test_get_web_el_height(psm_id, segment_id):
    segment = hc_variant.longitudinal_members()[psm_id].segments[segment_id - 1]
    dim = test_mesh.get_web_el_height(segment)
    print(dim)


def Test_get_min_flange_el_length():
    psm1 = PrimarySuppMem(1, BeamDirection.LONGITUDINAL, 0.1, hc_variant)
    ST24 = MaterialProperty(1, 210000, 0.3, 7850, 235, "ST24")
    beam_prop1 = FBBeamProperty(1, 1089, 10, ST24)
    # beam_prop2 = FBBeamProperty(2, 1089, 10, ST24)
    beam_prop2 = TBeamProperty(2, 1089, 10, 545, 40, ST24)
    segment1 = Segment(1, beam_prop1, psm1, psm1, psm1)
    segment2 = Segment(2, beam_prop2, psm1, psm1, psm1)

    min_dim = test_mesh.get_min_fl_el_len(segment1, segment2)
    print("Minimalna vrijednost je", min_dim)


def Test_get_min_flange_el_length_between_psm(member1_id, member2_id):
    member1 = hc_variant.longitudinal_members()[member1_id]
    member2 = hc_variant.longitudinal_members()[member2_id]
    min_dim = test_mesh.get_min_fl_el_len_between_psm(member1, member2)
    print("Najmanja vrijednost maksimalne duljine prirubnice svih segmenata između jakih nosača:", min_dim)


def Test_find_largest_divisor(length, max_val):
    divisor = test_mesh.find_largest_divisor(length, max_val)
    print("Najveći djelitelj broja", length, ", koji rezultira dimenzijom manjom ili jednakom", max_val, "je", divisor,
          ", što daje vrijednost", length / divisor)


def Test_transition_element_size_plating_zone(plate_id, segment_id):
    plate = hc_variant.plating()[plate_id]
    transition_dims = test_mesh.tr_element_size_plating_zone(plate, segment_id)
    print("Dimenzije prijelaznog elemenata na zoni oplate", plate_id, "uz segment", segment_id, ":", transition_dims)


def Test_psm_extent():
    longs = extents.longitudinal_psm_extent()
    trans = extents.transverse_psm_extent()
    print("Prema unesenoj osi simetrije", extents.axis_of_symm, ", od ukupno",
          len(hc_variant.longitudinal_members()), "jakih uzdužnih nosača na modelu, radi se mreža za njih", len(longs))
    print("Prema unesenoj osi simetrije", extents.axis_of_symm, ", od ukupno",
          len(hc_variant.transverse_members()), "jakih poprečnih nosača na modelu, radi se mreža za njih", len(trans))


def Test_segment_extent():
    extents.grillage_segment_extent()
    print("Ukupno identificiranih segmenata za punu izradu mreže:", len(extents.full_segments))
    for item in extents.full_segments.items():
        key, segment = item
        print("  ", key, "Jaki", segment.primary_supp_mem.direction.name,
              segment.beam_prop.beam_type.name, "nosač broj", segment.primary_supp_mem.id,
              ", segmenta broj", segment.id)
    print("Ukupno identificiranih segmenata za polovičnu izradu mreže:", len(extents.half_segments))
    for item in extents.half_segments.items():
        key, segment = item
        print("  ", key, "Jaki", segment.primary_supp_mem.direction.name,
              segment.beam_prop.beam_type.name, "nosač broj", segment.primary_supp_mem.id,
              ", segmenta broj", segment.id)


def Test_all_segment_extent():
    print("Ukupno identificiranih segmenata za izradu mreže:", len(extents.all_segments))
    print("ID segmenata za izradu mreže:")
    for item in extents.all_segments.items():
        key, segment = item
        print(key, "Jaki", segment.primary_supp_mem.direction.name,
              segment.beam_prop.beam_type.name, "nosač broj", segment.primary_supp_mem.id,
              ", segmenta broj", segment.id)


def Test_grillage_plate_extent():
    print("Odabrana je", extents.axis_of_symm, "os simetrije.")

    print("Ukupno identificiranih zona za izradu mreže:", len(extents.all_plating_zones))
    print("ID zona oplate za izradu mreže:")
    for key in extents.all_plating_zones:
        print("Ključ:", key, ", ID zone oplate", extents.all_plating_zones[key].id)

    print("Od tih upisanih zona, dijele se na različite tipove:")

    print("Zone oplate za izradu pune mreže:")
    for key in extents.full_plate_zones:
        print("Ključ:", key, ", ID zone oplate", extents.all_plating_zones[key].id)

    print("Zone oplate za izradu polovične mreže, presječenih uzdužnom osi simetrije:")
    for key in extents.long_half_plate_zones:
        print("Ključ:", key, ", ID zone oplate", extents.all_plating_zones[key].id)

    print("ID zona oplate za izradu polovične mreže, presječenih poprečnom osi simetrije:")
    for key in extents.tran_half_plate_zones:
        print("Ključ:", key, ", ID zone oplate", extents.all_plating_zones[key].id)

    print("Zone oplate za izradu četvrtinske mreže:")
    for key in extents.quarter_plate_zone:
        print("Ključ:", key, ", ID zone oplate", extents.all_plating_zones[key].id)


def Test_PlatingZoneMesh(plate_id, split_along=AOS.NONE):
    plate = hc_variant.plating()[plate_id]

    extents.grillage_plate_extent()        # Izračun koje zone oplate se meshiraju
    PlateMesh(test_mesh, plate, 1, 1, 1, split_along).generate_mesh(grillage_fem)     # izrada mreže jedne zone oplate


def Test_plating_zones_ref_array():
    extents.grillage_plate_extent()
    arr = extents.plating_zones_ref_array

    n_redaka, n_stupaca = np.shape(arr)
    print(arr)
    print("Postoji", n_redaka, "redaka i", n_stupaca, "stupaca zona oplate koje se meshiraju")


def Test_get_plate_dim(plate_id):
    extents.grillage_plate_extent()
    plate = hc_variant.plating()[plate_id]
    # full_dim = plate.plate_dim_parallel_to_stiffeners() * 1000
    full_dim = test_mesh.get_reduced_plate_dim(plate)
    print("Puna dimenzija:", full_dim)
    print(extents.get_plate_dim(plate, full_dim))


def Test_calc_element_base_size():
    print("Osnovne dimenzije mreže dim_x i dim_y za sve stupce i retke zona oplate koji se meshiraju:")
    print(test_mesh.calc_element_base_size())


def Test_calculate_mesh_dimensions():
    x_dimenzije = test_mesh.plate_edge_node_x
    y_dimenzije = test_mesh.plate_edge_node_y

    print("Razmaci između čvorova (dimenzije elemenata) duž x osi:")
    for i in x_dimenzije.keys():
        print("  Zona oplate", i, ":", x_dimenzije[i])

    print("Razmaci između čvorova (dimenzije elemenata) duž y osi:")
    for i in y_dimenzije.keys():
        print("  Zona oplate", i, ":", y_dimenzije[i])


def Test_saved_plate_edge_node_dims():
    for plate_id in test_mesh.plate_edge_node_x.keys():
        print("Zona oplate:", plate_id, ", upisane x dimenzije:", test_mesh.plate_edge_node_x[plate_id])
    for plate_id in test_mesh.plate_edge_node_y.keys():
        print("Zona oplate:", plate_id, ", upisane y dimenzije:", test_mesh.plate_edge_node_y[plate_id])


def Test_Segment_element_generation(direction: BeamDirection, psm_id, segment_id):
    segment = None
    if direction == BeamDirection.LONGITUDINAL:
        segment = hc_variant.longitudinal_members()[psm_id].segments[segment_id - 1]
    elif direction == BeamDirection.TRANSVERSE:
        segment = hc_variant.transverse_members()[psm_id].segments[segment_id - 1]

    start_node_id = 1
    start_element_id = 1
    seg_mesh = SegmentMeshV1(test_mesh, segment, start_node_id, start_element_id)
    extents.generate_FEM_material(grillage_fem)
    extents.generate_FEM_plate_property(grillage_fem)
    seg_mesh.generate_mesh(grillage_fem)


def Test_edge_segment_node_generation(direction: BeamDirection, psm_id, segment_id):
    segment = None
    if direction == BeamDirection.LONGITUDINAL:
        segment = hc_variant.longitudinal_members()[psm_id].segments[segment_id - 1]
    elif direction == BeamDirection.TRANSVERSE:
        segment = hc_variant.transverse_members()[psm_id].segments[segment_id - 1]

    start_node_id = 1
    start_element_id = 1
    seg_mesh = SegmentMeshV1(test_mesh, segment, start_node_id, start_element_id)
    seg_mesh.get_plate_edge_nodes()
    seg_mesh.generate_web_nodes(grillage_fem)
    last_node = seg_mesh.generate_flange_nodes(grillage_fem, FlangeDirection.INWARD, start_node_id)
    print("ID koji se prenosi na idući segment: za čvor", last_node)


def Test_get_segment_element_property(direction: BeamDirection, psm_id, segment_id):
    if direction == BeamDirection.LONGITUDINAL:
        segment = hc_variant.longitudinal_members()[psm_id].segments[segment_id - 1]
    else:
        segment = hc_variant.transverse_members()[psm_id].segments[segment_id - 1]
    extents.generate_FEM_material(grillage_fem)
    extents.generate_FEM_plate_property(grillage_fem)

    start_node_id = 1
    start_element_id = 1
    seg_mesh = SegmentMeshV1(test_mesh, segment, start_node_id, start_element_id)
    web_fem_prop = seg_mesh.get_web_element_property_id(grillage_fem)
    flange_fem_prop = seg_mesh.get_flange_element_property_id(grillage_fem)

    print("  ", "Jaki", segment.primary_supp_mem.direction.name,
          segment.beam_prop.beam_type.name, "nosač broj", segment.primary_supp_mem.id,
          ", segment broj", segment.id)
    print("     ID upisanog svojstva struka u rječniku GeoFEM properties:", web_fem_prop.id)
    print("     Debljina materijala:", web_fem_prop.tp, "mm, materijal:", web_fem_prop.material.name)
    print("     ID upisanog svojstva prirubnice u rječniku GeoFEM properties:", flange_fem_prop.id)
    print("     Debljina materijala:", flange_fem_prop.tp, "mm, materijal:", flange_fem_prop.material.name)


def Test_ALL_segment_element_property():
    extents.generate_FEM_material(grillage_fem)
    extents.generate_FEM_plate_property(grillage_fem)
    for segment in extents.all_segments.values():
        start_node_id = 1
        start_element_id = 1

        seg_mesh = SegmentMeshV1(test_mesh, segment, start_node_id, start_element_id)
        web_fem_prop = seg_mesh.get_web_element_property_id(grillage_fem)
        flange_fem_prop = seg_mesh.get_flange_element_property_id(grillage_fem)

        print("  ", "Jaki", segment.primary_supp_mem.direction.name,
              segment.beam_prop.beam_type.name, "nosač broj", segment.primary_supp_mem.id,
              ", segment broj", segment.id)
        print("     ID BeamProperty u modelu:", segment.beam_prop.id)
        print("     ID upisanog svojstva struka u rječniku GeoFEM properties:", web_fem_prop.id)
        print("         Debljina materijala:", web_fem_prop.tp, "mm, materijal:", web_fem_prop.material.name)
        print("     ID upisanog svojstva prirubnice u rječniku GeoFEM properties:", flange_fem_prop.id)
        print("         Debljina materijala:", flange_fem_prop.tp, "mm, materijal:", flange_fem_prop.material.name, "\n")


def Test_aos_stiffener(plate_id):
    plate = hc_variant.plating()[plate_id]
    test_btw = extents.aos_between_stiffeners(plate)
    test_on = extents.aos_on_stiffener(plate)
    print("Test osi simetrije između ukrepa:", test_btw, ", test osi simetrije na ukrepi:", test_on)


def Test_aos_on_segment(direction: BeamDirection, psm_id, segment_id):
    segment = None
    if direction == BeamDirection.LONGITUDINAL:
        segment = hc_variant.longitudinal_members()[psm_id].segments[segment_id - 1]
    elif direction == BeamDirection.TRANSVERSE:
        segment = hc_variant.transverse_members()[psm_id].segments[segment_id - 1]
    test = extents.aos_on_segment(segment)
    print("Test prolazi li os simetrije uz segment:", test)


def Test_Model_check():
    check1 = ModelCheck(hc_variant)
    hc_variant.assign_symmetric_members()

    long_psm_symm = check1.longitudinal_psm_symmetry()
    tran_psm_symm = check1.transverse_psm_symmetry()
    central_long = check1.central_longitudinal()
    central_tran = check1.central_transversal()
    long_plate = check1.longitudinal_plate_symmetry()
    tran_plate = check1.transverse_plate_symmetry()
    long_symm = check1.long_symmetry_tests()
    tran_symm = check1.tran_symmetry_tests()
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


def Test_identify_split_element_zones():
    print("Uzdužna os simetrije prolazi između ukrepa na zonama:", extents.long_e_split_zone)
    print("Poprečna os simetrije prolazi između ukrepa na zonama:", extents.tran_e_split_zone)


def Test_get_split_elements_number(plate_id):
    plate = hc_variant.plating()[plate_id]
    ldim = test_mesh.get_long_split_element_num(plate)
    tdim = test_mesh.get_tran_split_element_num(plate)
    print("Uzdužna os simeterije prolazi između ukrepa, siječe broj elemenata na pola i stavlja element dimenzije", ldim)
    print("Poprečna os simeterije prolazi između ukrepa, siječe broj elemenata na pola i ostavlja element dimenzije", tdim)


def Test_generate_inward_flange_nodes(direction: BeamDirection, psm_id, segment_id, flange_dir: FlangeDirection):
    segment = None
    if direction == BeamDirection.LONGITUDINAL:
        segment = hc_variant.longitudinal_members()[psm_id].segments[segment_id - 1]
    elif direction == BeamDirection.TRANSVERSE:
        segment = hc_variant.transverse_members()[psm_id].segments[segment_id - 1]

    start_node_id = 1
    start_element_id = 1
    seg_mesh = SegmentMeshV1(test_mesh, segment, start_node_id, start_element_id)
    seg_mesh.get_plate_edge_nodes()
    seg_mesh.generate_web_nodes(grillage_fem)
    seg_mesh.generate_flange_nodes(grillage_fem, flange_dir, start_node_id)


def Test_PlatingZoneMesh_beam_elements(plate_id, split_along=AOS.NONE):
    plate = hc_variant.plating()[plate_id]

    extents.grillage_plate_extent()        # Izračun koje zone oplate se meshiraju
    PlateMesh(test_mesh, plate, 1, 1, 1, split_along).generate_beam_elements(grillage_fem)


def Test_flange_ref_array(direction: BeamDirection, psm_id, segment_id):
    segment = None
    if direction == BeamDirection.LONGITUDINAL:
        segment = hc_variant.longitudinal_members()[psm_id].segments[segment_id - 1]
    elif direction == BeamDirection.TRANSVERSE:
        segment = hc_variant.transverse_members()[psm_id].segments[segment_id - 1]

    start_node_id = 1
    start_element_id = 1
    seg_mesh = SegmentMeshV1(test_mesh, segment, start_node_id, start_element_id)
    seg_mesh.get_plate_edge_nodes()
    flange_start_node = 56
    ref_array = seg_mesh.ref_flange_node_ID_array(flange_start_node)
    print(ref_array)


def Test_get_opposite_flange_width(direction: BeamDirection, psm_id, segment_id):
    segment = None
    if direction == BeamDirection.LONGITUDINAL:
        segment = hc_variant.longitudinal_members()[psm_id].segments[segment_id - 1]
    elif direction == BeamDirection.TRANSVERSE:
        segment = hc_variant.transverse_members()[psm_id].segments[segment_id - 1]
    bf1, bf2 = test_mesh.get_opposite_flange_width(segment)
    print(bf1, bf2)


def Test_generate_element_row(direction: BeamDirection, psm_id, segment_id):
    segment = None
    if direction == BeamDirection.LONGITUDINAL:
        segment = hc_variant.longitudinal_members()[psm_id].segments[segment_id - 1]
    elif direction == BeamDirection.TRANSVERSE:
        segment = hc_variant.transverse_members()[psm_id].segments[segment_id - 1]

    start_node_id = 1
    start_element_id = 101
    SegmentMeshV2(test_mesh, segment, start_node_id, start_element_id).generate_element_row(1, start_element_id)


def Test_T_Profile_BeamProperty(hw, tw, bf, tf):
    test_beam = gfe.T_Profile_BeamProperty()
    test_beam.hw = hw
    test_beam.tw = tw
    test_beam.bf = bf
    test_beam.tf = tf
    print(" Visina struka,                          hw = ", test_beam.hw, "mm")
    print(" Debljina struka,                        tw = ", test_beam.tw, "mm")
    print(" Sirina prirubnice,                      bf = ", test_beam.bf, "mm")
    print(" Debljina prirubnice,                    tf = ", test_beam.tf, "mm")
    print(" Povrsina T profila,                      A = ", test_beam.area, "mm2")
    print(" Polozaj neutralne linije,             Z_na = ", test_beam.z_na, "mm")
    print(" Moment inercije (net),                  Iy = ", test_beam.Iy, "mm4")


def Test_Hat_Profile_BeamProperty(h, t, bf, fi):
    test_beam = gfe.Hat_Profile_BeamProperty()
    test_beam.h = h
    test_beam.t = t
    test_beam.bf = bf
    test_beam.fi = fi

    print(" Visina profila,                                     h = ", test_beam.h, "mm")
    print(" Debljina struka i prirubnice,                       t = ", test_beam.t, "mm")
    print(" Sirina prirubnice,                                 bf = ", test_beam.bf, "mm")
    print(" Kut nagiba struka profila,                         fi = ", test_beam.fi, "°")
    print(" Povrsina Hat profila,                               A = ", test_beam.area, "mm2")
    print(" Težište Hat profila od simetrale prirubnice,     z_NA = ", test_beam.z_na, "mm")
    print(" Moment inercije (net),                             Iy = ", test_beam.Iy, "mm4")


def Test_Bulb_Profile_BeamProperty(hw_ekv, tw_ekv, bf_ekv, tf_ekv):
    test_beam = gfe.Bulb_Profile_BeamProperty()
    test_beam.hw_ekv = hw_ekv
    test_beam.tw_ekv = tw_ekv
    test_beam.bf_ekv = bf_ekv
    test_beam.tf_ekv = tf_ekv
    print(" Visina struka,                          hw = ", test_beam.hw_ekv, "mm")
    print(" Debljina struka,                        tw = ", test_beam.tw_ekv, "mm")
    print(" Sirina prirubnice,                      bf = ", test_beam.bf_ekv, "mm")
    print(" Debljina prirubnice,                    tf = ", test_beam.tf_ekv, "mm")
    print(" Povrsina T profila,                      A = ", test_beam.area, "mm2")
    print(" Polozaj neutralne linije,             Z_na = ", test_beam.z_na, "mm")
    print(" Moment inercije (net),                  Iy = ", test_beam.Iy, "mm4")


def Test_FB_Profile_BeamProperty(hw, tw):
    test_beam = gfe.FB_Profile_BeamProperty()
    test_beam.hw = hw
    test_beam.tw = tw
    print(" Visina struka,                          hw = ", test_beam.hw, "mm")
    print(" Debljina struka,                        tw = ", test_beam.tw, "mm")
    print(" Povrsina FB profila,                     A = ", test_beam.area, "mm2")
    print(" Polozaj neutralne linije,             Z_na = ", test_beam.z_na, "mm")
    print(" Moment inercije (net),                  Iy = ", test_beam.Iy, "mm4")


def Test_GeoFEM_materials():
    extents.generate_FEM_material(grillage_fem)
    num = grillage_fem.num_materials
    print("Lista upisanih materijala u GeoFEM materials:")
    for material_id in range(1, num + 1):
        mat = grillage_fem.getMaterial(material_id)
        print(mat.id, mat.name, mat.E, mat.rho, mat.ReH, mat.ni)


def Test_GeoFEM_plate_property():
    extents.generate_FEM_material(grillage_fem)
    extents.generate_FEM_plate_property(grillage_fem)
    num = grillage_fem.num_properties
    print("Lista upisanih plate property u GeoFEM properties:")
    for prop_id in range(1, num + 1):
        prop = grillage_fem.getProperty(prop_id)
        print(prop.id, " tp =", prop.tp, prop.material.name)


def Test_GeoFEM_T_L_beam_property():
    extents.generate_FEM_material(grillage_fem)
    extents.generate_FEM_beam_property(grillage_fem)
    num = grillage_fem.num_properties
    print("Ukupno upisanih:", num)
    print("Lista upisanih beam property u GeoFEM properties:")
    for prop_id in range(1, num + 1):
        prop = grillage_fem.getProperty(prop_id)
        print(prop.id, prop.material.name, prop.name, ", hw =", prop.hw, ", tw =", prop.tw,
              ", bf =", prop.bf, ", tf =", prop.tf, ", Iy =", prop.Iy, ", A =", prop.area)


def Test_GeoFEM_Bulb_beam_property():
    extents.generate_FEM_material(grillage_fem)
    extents.generate_FEM_beam_property(grillage_fem)
    num = grillage_fem.num_properties
    print("Ukupno upisanih property:", num)

    for prop_id in range(1, num + 1):
        prop = grillage_fem.getProperty(prop_id)
        num_vals = prop.get_num_vals()
        print(prop.id, prop.name, prop.material.name)
        print("   ", "Iy =", prop.Iy, "mm4")
        print("   ",  "A =", prop.area, "mm2")
        for descriptor in range(0, num_vals):
            print("   ", prop.get_desc_names()[descriptor], "=", prop.descriptors[descriptor])
