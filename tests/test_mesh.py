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


hc_var = 5
filename = str("../grillage savefiles/hc_var_") + str(hc_var) + str("_savefile.gin")
hc_variant = GrillageModelData(filename).read_file()
print("Testing FE mesh for grillage variant", hc_var)

# extents = MeshExtent(hc_variant, AOS.NONE)    # Calculate mesh extents with Axis of Symmetry override
extents = MeshExtent(hc_variant)                # Calculate mesh extents with automatic Axis of Symmetry discovery

# test_mesh_size = ElementSizeV1(extents)      # Calculate mesh dimensions for mesh variant V1
test_mesh_size = ElementSizeV2(extents)     # Calculate mesh dimensions for mesh variant V2

# Mesh Control
test_mesh_size.min_num_ebs = 1               # Minimum number of elements between stiffeners; default = 1
test_mesh_size.min_num_eweb = 3              # Minimum number of elements along psm web height; default = 3
test_mesh_size.num_eaf = 1                   # Number of elements across the psm flange; default = 1
test_mesh_size.flange_aspect_ratio = 8       # Max flange aspect ratio; default = 8
test_mesh_size.plate_aspect_ratio = 4        # Max plate aspect ratio; default = 4
test_mesh_size.des_plate_aspect_ratio = 3    # Desired plate aspect ratio; default = 3
gm_test = GrillageMesh(test_mesh_size)

# grillage_test_mesh = gm_test.generate_grillage_mesh_v1("test_mesh")


def generate_test_mesh_v1():
    start = timer()
    mesh_v1 = ElementSizeV1(extents)  # Calculate mesh dimensions for mesh variant V2

    # Mesh Control
    mesh_v1.min_num_ebs = 1  # Minimum number of elements between stiffeners; default = 1
    mesh_v1.min_num_eweb = 3  # Minimum number of elements along psm web height; default = 3
    mesh_v1.num_eaf = 1  # Number of elements across the psm flange; default = 1
    mesh_v1.flange_aspect_ratio = 8  # Max flange aspect ratio; default = 8
    mesh_v1.plate_aspect_ratio = 4  # Max plate aspect ratio; default = 4
    mesh_v1.des_plate_aspect_ratio = 3  # Desired plate aspect ratio; default = 3
    test_mesh_v1 = GrillageMesh(mesh_v1)

    grillage_test_mesh = test_mesh_v1.generate_grillage_mesh_v1("test_mesh_v1")
    grillage_test_mesh.merge_coincident_nodes()
    grillage_test_mesh.merge_coincident_elements()
    end = timer()
    print("Mesh generation time:", end - start, "s")
    # grillage_test_mesh.full_model_node_overlap_check()
    return grillage_test_mesh


def generate_test_mesh_v2():
    start = timer()
    mesh_v2 = ElementSizeV2(extents)  # Calculate mesh dimensions for mesh variant V2

    # Mesh Control
    mesh_v2.min_num_ebs = 1  # Minimum number of elements between stiffeners; default = 1
    mesh_v2.min_num_eweb = 3  # Minimum number of elements along psm web height; default = 3
    mesh_v2.num_eaf = 1  # Number of elements across the psm flange; default = 1
    mesh_v2.flange_aspect_ratio = 8  # Max flange aspect ratio; default = 8
    mesh_v2.plate_aspect_ratio = 4  # Max plate aspect ratio; default = 4
    mesh_v2.des_plate_aspect_ratio = 3  # Desired plate aspect ratio; default = 3
    test_mesh_v2 = GrillageMesh(mesh_v2)

    grillage_test_mesh = test_mesh_v2.generate_grillage_mesh_v2("test_mesh_v2")

    grillage_test_mesh.merge_coincident_nodes()
    grillage_test_mesh.merge_coincident_elements()
    end = timer()
    print("Mesh generation time:", end - start, "s")
    # grillage_test_mesh.full_model_node_overlap_check()
    return grillage_test_mesh


def Test_calculate_mesh_dimensions():
    test_mesh_size.calculate_mesh_dimensions()
    for plate in test_mesh_size.mesh_extent.all_plating_zones.values():
        dim_x = test_mesh_size.plate_edge_node_spacing_x(plate)
        dim_y = test_mesh_size.plate_edge_node_spacing_y(plate)
        print("Plating zone ID", plate.id, ":")
        print(" Distance between plating edge nodes (dim_x):", dim_x)
        print(" Distance between plating edge nodes (dim_y):", dim_y)


def Test_GeoFEM_T_L_beam_property(grillage_test_mesh):
    num = grillage_test_mesh.num_properties
    print("Ukupno upisanih property:", num)

    for prop_id in range(1, num + 1):
        prop = grillage_test_mesh.getProperty(prop_id)
        if not isinstance(prop, gfe.PlateProperty):
            num_vals = prop.get_num_vals()
            print(prop.id, prop.name, prop.material.name)
            print("   ", "Iy =", prop.Iy, "mm4")
            print("   ",  "A =", prop.area, "mm2")
            print("   ",  "z_na =", prop.z_na, "mm2")

            for descriptor in range(0, num_vals):
                print("   ", prop.get_desc_names()[descriptor], "=", prop.descriptors[descriptor])


def Test_current():
    start = timer()
    grillage_test_mesh = gm_test.generate_grillage_mesh_v1("test_mesh")
    grillage_test_mesh.merge_coincident_nodes()
    grillage_test_mesh.identify_plating_nodes()
    # grillage_test_mesh.check_node_overlap_np()
    # grillage_test_mesh.check_node_overlap()
    # grillage_test_mesh.full_model_node_overlap_check()

    # for node in grillage_test_mesh.nodes.values():
    #     print(node.p)

    end = timer()
    print("Mesh generation time:", end - start, "s")


def Test_MeshVariant_V2_basesize():
    test_mesh_size.mesh_extent.grillage_mesh_extent()
    test_mesh_size.calc_element_base_size_mesh()

    for plate in extents.all_plating_zones.values():
        dim_x = test_mesh_size.get_base_dim_x(plate)
        dim_y = test_mesh_size.get_base_dim_y(plate)
        print("Zona oplate ID:", plate.id, ",   dim_x =", "{:.2f}".format(dim_x), "mm", ",   dim_y =", "{:.2f}".format(dim_y), "mm")


def Test_MeshVariant_V2_transition():
    test_mesh_size.mesh_extent.grillage_mesh_extent()
    test_mesh_size.calc_element_base_size_mesh()

    """
    # TRANSITION DIM X,X; Y,Y
    for plate in extents.all_plating_zones.values():
        dimx1, dimx2 = test_mesh_size.transition_dim_x(plate)
        dimy1, dimy2 = test_mesh_size.transition_dim_y(plate)
        print("Dimenzije prijelaznog elemenata na zoni oplate", plate.id, "X:", dimx1, dimx2, ", Y:", dimy1, dimy2)
    """
    for plate in extents.all_plating_zones.values():
        tr_dim_x = test_mesh_size.get_tr_dim_x(plate)
        tr_dim_y = test_mesh_size.get_tr_dim_y(plate)
        print("Dimenzije globalno usklađenog prijelaznog elemenata na zoni oplate", plate.id, "\n",
              "X:", tr_dim_x, "\n", "Y:", tr_dim_y)

    # TRANSITION FLANGE ELEMENT DIMENSIONS
    for segment in test_mesh_size.mesh_extent.all_segments.values():
        fl_tr_1, fl_tr_2 = test_mesh_size.flange_transition_dim(segment)
        num_end1, num_end2 = test_mesh_size.opposite_flange_element_num(segment)
        psm_id = segment.primary_supp_mem.id
        direct = segment.primary_supp_mem.direction.name
        psm_type = segment.beam_prop.beam_type.name
        print("Jaki", direct, psm_type, "nosač broj", psm_id, ", segmenta broj", segment.id,
              " , dimenzije prijelaznih elemenata prirubnice:", fl_tr_1, fl_tr_2,
              ", broj elemenata prirubnice:", num_end1, num_end2)


    # OPPOSITE FLANGE WIDTH
    print("\n", "OPPOSITE FLANGE WIDTH")
    for segment in test_mesh_size.mesh_extent.all_segments.values():
        bf_max1, bf_max2 = test_mesh_size.opposite_flange_width(segment)
        psm_id = segment.primary_supp_mem.id
        direct = segment.primary_supp_mem.direction.name
        psm_type = segment.beam_prop.beam_type.name
        is_central = hc_variant.central_segment(segment)
        print("Jaki", direct, psm_type, "nosač broj", psm_id, ", segmenta broj", segment.id,
              " , dimenzije:", bf_max1, bf_max2, ", na sredini:", is_central)

def Test_MeshVariant_V2_element_number():
    test_mesh_size.mesh_extent.grillage_mesh_extent()
    test_mesh_size.calc_element_base_size_mesh()
    """
    print("**** BROJ ELEMENATA NA OPLATI ****")
    for plate in extents.all_plating_zones.values():
        n_x, n_y = test_mesh_size.get_base_element_number(plate)

        print("Zona oplate:", plate.id, "Broj elemenata osnovnih dimenzija po x:", n_x, ", po y:", n_y)

        ldim = test_mesh_size.get_long_split_element_num(plate)
        tdim = test_mesh_size.get_tran_split_element_num(plate)
        print(" Uzdužna os simeterije prolazi između ukrepa, siječe broj elemenata na pola", ldim)
        print(" Poprečna os simeterije prolazi između ukrepa, siječe broj elemenata na pola", tdim)
    """

    print("\n", "**** BROJ ELEMENATA NA PRIRUBNICI ****")
    for segment in test_mesh_size.mesh_extent.all_segments.values():
        psm_id = segment.primary_supp_mem.id
        direct = segment.primary_supp_mem.direction.name
        psm_type = segment.beam_prop.beam_type.name
        base_dim_num = test_mesh_size.flange_base_element_num(segment)
        flange_tr_num = test_mesh_size.get_flange_transition_num(segment)
        print("Jaki", direct, psm_type, "nosač broj", psm_id, ", segmenta broj", segment.id,
              " , broj elemenata osnovne mreže prirubnice:", base_dim_num,
              ", broj prijelaznih elemenata:", flange_tr_num)


def Test_MeshVariant_V2_flange_edge_nodes():
    test_mesh_size.mesh_extent.grillage_mesh_extent()
    test_mesh_size.calc_element_base_size_mesh()
    for segment in test_mesh_size.mesh_extent.all_segments.values():
        psm_id = segment.primary_supp_mem.id
        direct = segment.primary_supp_mem.direction.name
        psm_type = segment.beam_prop.beam_type.name
        dims = test_mesh_size.flange_edge_node_spacing(segment)
        print("Jaki", direct, psm_type, "nosač broj", psm_id, ", segmenta broj", segment.id,
              " , dimenzije mreže prirubnice:", dims)

def TestFullModelOverlap():
    start = timer()
    grillage_test_mesh = gm_test.generate_grillage_mesh_v1("test_mesh")
    grillage_test_mesh.full_model_node_overlap_check()
    end = timer()
    print("Full model overlap check time:", end - start, "s")


def Test_get_reduced_plate_dim(plate_id):
    plate = hc_variant.plating()[plate_id]
    reduced_dim = test_mesh_size.get_reduced_plate_dim(plate)
    print("Reducirana dimenzija paralalna s ukrepama na zoni oplate", plate_id, "iznosi:", reduced_dim)


def Test_find_closest_divisor(length, spacing):
    div = test_mesh_size.find_closest_divisor(length, spacing)
    if div is None:
        print("Nije pronađen djelitelj!")
    else:
        print("Pronađen je djelitelj! Broj", length, "se može podijeliti na", div, "dijelova, tako da svaki iznosi", length/div,
              ", što bi trebalo biti blizu zadanog broja", spacing)


def Test_element_size_para_to_stiffeners(plate_id):
    plate = hc_variant.plating()[plate_id]
    okomita_dimenzija = test_mesh_size.element_size_perp_to_stiffeners(plate)
    plate_dim = test_mesh_size.get_reduced_plate_dim(plate)
    paralelna_dimenzija = test_mesh_size.element_size_para_to_stiffeners(plate, plate_dim)
    n_elem = test_mesh_size.min_num_ebs
    ar = MeshSize.element_aspect_ratio(okomita_dimenzija, paralelna_dimenzija)
    print("Dimenzije quad elementa oplate uz", n_elem, "element između ukrepa, okomito:", okomita_dimenzija,
          "mm, paralelno:", paralelna_dimenzija, "mm", ", aspektni odnos:", ar)


def Test_get_flange_el_length(direction: BeamDirection, psm_id, segment_id):
    segment = None

    if direction == BeamDirection.LONGITUDINAL:
        segment = hc_variant.longitudinal_members()[psm_id].segments[segment_id]

    elif direction == BeamDirection.TRANSVERSE:
        segment = hc_variant.transverse_members()[psm_id].segments[segment_id]

    dim = test_mesh_size.get_flange_el_length(segment)
    print("Maksimalna dimenzija elementa prirubnice prema aspektnom odnosu:", dim, "mm")


def Test_element_size_plating_zone(plate_id):
    plate = hc_variant.plating()[plate_id]
    plate_dim = test_mesh_size.get_reduced_plate_dim(plate)
    dims = test_mesh_size.element_size_plating_zone(plate, plate_dim)
    print("Osnovne dimenzije elementa na zoni oplate", plate_id)
    dim_x = dims[0]
    dim_y = dims[1]
    print("dim_x =", dim_x, ", dim_y =", dim_y)


def Test_element_aspect_ratio(dim_x, dim_y):
    ar = test_mesh_size.element_aspect_ratio(dim_x, dim_y)
    print("Aspektni odnos =", ar)


def Test_refine_plate_element(length, dim_limit):
    dim = test_mesh_size.refine_plate_element(length, dim_limit)
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
            plate_dim = test_mesh_size.get_reduced_plate_dim(plate)
            polje1[stupac, redak] = test_mesh_size.element_size_plating_zone(plate, plate_dim)[0]
            plate_id += 1
    print("Sve x dimenzije elemenata: \n", polje1, "\n")

    # dimenzije karakteristicnih elemenata na cijeloj oplati po y osi
    plate_id = 1
    for stupac in range(0, n_y):
        for redak in range(0, n_x):
            plate = hc_variant.plating()[plate_id]
            plate_dim = test_mesh_size.get_reduced_plate_dim(plate)
            polje1[stupac, redak] = test_mesh_size.element_size_plating_zone(plate, plate_dim)[1]
            plate_id += 1
    print("Sve y dimenzije elemenata: \n", polje1)


def Test_element_size_mesh():
    print("Konačno odabrane dimenzije mreže po x:", test_mesh_size.mesh_dim_x)
    print("Konačno odabrane dimenzije mreže po y:", test_mesh_size.mesh_dim_y)


def Test_get_flange_el_width(psm_id, segment_id):
    segment = hc_variant.longitudinal_members()[psm_id].segments[segment_id - 1]

    print("Jaki uzduzni nosac", psm_id, ", ID segmenta:", segment.id,
          " ,BeamProperty ID:", segment.beam_prop.id,
          ", tip profila: ", segment.beam_prop.beam_type,
          ", tf =", segment.beam_prop.tf,
          ", dimenzija elementa prirubnice po širini:", test_mesh_size.get_flange_el_width(segment), "mm")


def Test_all_plating_zones_mesh_dimensions():
    extents.grillage_mesh_extent()
    test_mesh_size.calc_element_base_size_mesh()
    for plate in extents.all_plating_zones.values():
        dim_x = test_mesh_size.get_base_dim_x(plate)
        dim_y = test_mesh_size.get_base_dim_y(plate)
        print("Zona oplate ID:", plate.id, ",   dim_x =", "{:.2f}".format(dim_x), "mm", ",   dim_y =", "{:.2f}".format(dim_y), "mm")


def Test_get_tr_dim_x(plate_id):
    # plate = mesh1.plating_zones[plate_id]
    plate = hc_variant.plating()[plate_id]
    dims = test_mesh_size.get_tr_dim_x(plate)
    print("Dimenzije x prijelaznih elemenata na zoni oplate", plate_id, ":", dims)


def Test_get_tr_dim_y(plate_id):
    plate = hc_variant.plating()[plate_id]
    dims = test_mesh_size.get_tr_dim_y(plate)
    print("Dimenzije y prijelaznih elemenata na zoni oplate", plate_id, ":", dims)


def Test_get_tr_element_num(plate_id):
    plate = hc_variant.plating()[plate_id]
    poprecni_segment = test_mesh_size.get_long_tr_element_num(plate, plate.trans_seg1)
    uzduzni_segment = test_mesh_size.get_long_tr_element_num(plate, plate.trans_seg2)
    print("Broj prijelaznih elemenata na zoni oplate", plate_id, "uz elemente prirubnice uzdužnog segmenta duž osi y:", uzduzni_segment,
          ", uz elemente prirubnice poprečnog segmenta duž osi x:", poprecni_segment)


def Test_get_element_number(plate_id):
    plate = hc_variant.plating()[plate_id]
    n_elemenata = test_mesh_size.get_base_element_number(plate)
    print("Broj elemenata osnovnih dimenzija:", plate_id, "po x:", n_elemenata[0], ", po y:", n_elemenata[1])


def Test_get_all_mesh_dim():
    for plate in extents.all_plating_zones.values():
        print("Dimenzije svih konačnih elemenata redom za zonu oplate", plate.id)
        print("Dimenzije x:", test_mesh_size.plate_edge_node_spacing_x(plate))
        print("Dimenzije y:", test_mesh_size.plate_edge_node_spacing_y(plate))
        print("\n")


def Test_get_web_el_height(psm_id, segment_id):
    segment = hc_variant.longitudinal_members()[psm_id].segments[segment_id - 1]
    dim = test_mesh_size.get_web_el_height(segment)
    print(dim)


def Test_get_min_flange_el_length():
    psm1 = PrimarySuppMem(1, BeamDirection.LONGITUDINAL, 0.1, hc_variant)
    ST24 = MaterialProperty(1, 210000, 0.3, 7850, 235, "ST24")
    beam_prop1 = FBBeamProperty(1, 1089, 10, ST24)
    # beam_prop2 = FBBeamProperty(2, 1089, 10, ST24)
    beam_prop2 = TBeamProperty(2, 1089, 10, 545, 40, ST24)
    segment1 = Segment(1, beam_prop1, psm1, psm1, psm1)
    segment2 = Segment(2, beam_prop2, psm1, psm1, psm1)

    min_dim = test_mesh_size.get_min_fl_el_len(segment1, segment2)
    print("Minimalna vrijednost je", min_dim)


def Test_get_min_flange_el_length_between_psm(member1_id, member2_id):
    member1 = hc_variant.longitudinal_members()[member1_id]
    member2 = hc_variant.longitudinal_members()[member2_id]
    min_dim = test_mesh_size.get_min_fl_el_len_between_psm(member1, member2)
    print("Najmanja vrijednost maksimalne duljine prirubnice svih segmenata između jakih nosača:", min_dim)


def Test_find_largest_divisor(length, max_val):
    divisor = test_mesh_size.find_largest_divisor(length, max_val)
    print("Najveći djelitelj broja", length, ", koji rezultira dimenzijom manjom ili jednakom", max_val, "je", divisor,
          ", što daje vrijednost", length / divisor)


def Test_transition_element_size_plating_zone(plate_id, segment_id):
    plate = hc_variant.plating()[plate_id]
    transition_dims = test_mesh_size.tr_element_size_plating_zone(plate, segment_id)
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
    PlateMesh(test_mesh_size, plate, 1, 1, 1, split_along).generate_mesh(grillage_test_mesh)     # izrada mreže jedne zone oplate


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
    full_dim = test_mesh_size.get_reduced_plate_dim(plate)
    print("Puna dimenzija:", full_dim)
    print(extents.get_plate_dim(plate, full_dim))


def Test_calc_element_base_size():
    print("Osnovne dimenzije mreže dim_x i dim_y za sve stupce i retke zona oplate koji se meshiraju:")
    print(test_mesh_size.calc_element_base_size())


def Test_Segment_element_generation(direction: BeamDirection, psm_id, segment_id):
    segment = None
    if direction == BeamDirection.LONGITUDINAL:
        segment = hc_variant.longitudinal_members()[psm_id].segments[segment_id - 1]
    elif direction == BeamDirection.TRANSVERSE:
        segment = hc_variant.transverse_members()[psm_id].segments[segment_id - 1]

    start_node_id = 1
    start_element_id = 1
    seg_mesh = SegmentMeshV1(test_mesh_size, segment, start_node_id, start_element_id)
    extents.generate_FEM_material(grillage_test_mesh)
    extents.generate_FEM_plate_property(grillage_test_mesh)
    seg_mesh.generate_mesh(grillage_test_mesh)


def Test_edge_segment_node_generation(direction: BeamDirection, psm_id, segment_id):
    segment = None
    if direction == BeamDirection.LONGITUDINAL:
        segment = hc_variant.longitudinal_members()[psm_id].segments[segment_id - 1]
    elif direction == BeamDirection.TRANSVERSE:
        segment = hc_variant.transverse_members()[psm_id].segments[segment_id - 1]

    start_node_id = 1
    start_element_id = 1
    seg_mesh = SegmentMeshV1(test_mesh_size, segment, start_node_id, start_element_id)
    seg_mesh.get_plate_edge_nodes()
    seg_mesh.generate_web_nodes(grillage_test_mesh)
    last_node = seg_mesh.generate_flange_nodes(grillage_test_mesh, FlangeDirection.INWARD, start_node_id)
    print("ID koji se prenosi na idući segment: za čvor", last_node)


def Test_get_segment_element_property(direction: BeamDirection, psm_id, segment_id):
    if direction == BeamDirection.LONGITUDINAL:
        segment = hc_variant.longitudinal_members()[psm_id].segments[segment_id - 1]
    else:
        segment = hc_variant.transverse_members()[psm_id].segments[segment_id - 1]
    extents.generate_FEM_material(grillage_test_mesh)
    extents.generate_FEM_plate_property(grillage_test_mesh)

    start_node_id = 1
    start_element_id = 1
    seg_mesh = SegmentMeshV1(test_mesh_size, segment, start_node_id, start_element_id)
    web_fem_prop = seg_mesh.get_web_element_property_id(grillage_test_mesh)
    flange_fem_prop = seg_mesh.get_flange_element_property_id(grillage_test_mesh)

    print("  ", "Jaki", segment.primary_supp_mem.direction.name,
          segment.beam_prop.beam_type.name, "nosač broj", segment.primary_supp_mem.id,
          ", segment broj", segment.id)
    print("     ID upisanog svojstva struka u rječniku GeoFEM properties:", web_fem_prop.id)
    print("     Debljina materijala:", web_fem_prop.tp, "mm, materijal:", web_fem_prop.material.name)
    print("     ID upisanog svojstva prirubnice u rječniku GeoFEM properties:", flange_fem_prop.id)
    print("     Debljina materijala:", flange_fem_prop.tp, "mm, materijal:", flange_fem_prop.material.name)


def Test_ALL_segment_element_property():
    extents.generate_FEM_material(grillage_test_mesh)
    extents.generate_FEM_plate_property(grillage_test_mesh)
    for segment in extents.all_segments.values():
        start_node_id = 1
        start_element_id = 1

        seg_mesh = SegmentMeshV1(test_mesh_size, segment, start_node_id, start_element_id)
        web_fem_prop = seg_mesh.get_web_element_property_id(grillage_test_mesh)
        flange_fem_prop = seg_mesh.get_flange_element_property_id(grillage_test_mesh)

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
    ldim = test_mesh_size.get_long_split_element_num(plate)
    tdim = test_mesh_size.get_tran_split_element_num(plate)
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
    seg_mesh = SegmentMeshV1(test_mesh_size, segment, start_node_id, start_element_id)
    seg_mesh.get_plate_edge_nodes()
    seg_mesh.generate_web_nodes(grillage_test_mesh)
    seg_mesh.generate_flange_nodes(grillage_test_mesh, flange_dir, start_node_id)


def Test_PlatingZoneMesh_beam_elements(plate_id, split_along=AOS.NONE):
    plate = hc_variant.plating()[plate_id]

    extents.grillage_plate_extent()        # Izračun koje zone oplate se meshiraju
    PlateMesh(test_mesh_size, plate, 1, 1, 1, split_along).generate_beam_elements(grillage_test_mesh)


def Test_flange_ref_array(direction: BeamDirection, psm_id, segment_id):
    segment = None
    if direction == BeamDirection.LONGITUDINAL:
        segment = hc_variant.longitudinal_members()[psm_id].segments[segment_id - 1]
    elif direction == BeamDirection.TRANSVERSE:
        segment = hc_variant.transverse_members()[psm_id].segments[segment_id - 1]

    start_node_id = 1
    start_element_id = 1
    seg_mesh = SegmentMeshV1(test_mesh_size, segment, start_node_id, start_element_id)
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
    bf1, bf2 = test_mesh_size.get_opposite_flange_width(segment)
    print(bf1, bf2)


def Test_generate_element_row(direction: BeamDirection, psm_id, segment_id):
    segment = None
    if direction == BeamDirection.LONGITUDINAL:
        segment = hc_variant.longitudinal_members()[psm_id].segments[segment_id - 1]
    elif direction == BeamDirection.TRANSVERSE:
        segment = hc_variant.transverse_members()[psm_id].segments[segment_id - 1]

    start_node_id = 1
    start_element_id = 101
    SegmentMeshV2(test_mesh_size, segment, start_node_id, start_element_id).tr_web_element_row(start_element_id)


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
    extents.generate_FEM_material(grillage_test_mesh)
    num = grillage_test_mesh.num_materials
    print("Lista upisanih materijala u GeoFEM materials:")
    for material_id in range(1, num + 1):
        mat = grillage_test_mesh.getMaterial(material_id)
        print(mat.id, mat.name, mat.E, mat.rho, mat.ReH, mat.ni)


def Test_GeoFEM_plate_property():
    extents.generate_FEM_material(grillage_test_mesh)
    extents.generate_FEM_plate_property(grillage_test_mesh)
    num = grillage_test_mesh.num_properties
    print("Lista upisanih plate property u GeoFEM properties:")
    for prop_id in range(1, num + 1):
        prop = grillage_test_mesh.getProperty(prop_id)
        print(prop.id, " tp =", prop.tp, prop.material.name)


def Test_GeoFEM_Bulb_beam_property():
    extents.generate_FEM_material(grillage_test_mesh)
    extents.generate_FEM_beam_property(grillage_test_mesh)
    num = grillage_test_mesh.num_properties
    print("Ukupno upisanih property:", num)

    for prop_id in range(1, num + 1):
        prop = grillage_test_mesh.getProperty(prop_id)
        num_vals = prop.get_num_vals()
        print(prop.id, prop.name, prop.material.name)
        print("   ", "Iy =", prop.Iy, "mm4")
        print("   ",  "A =", prop.area, "mm2")
        for descriptor in range(0, num_vals):
            print("   ", prop.get_desc_names()[descriptor], "=", prop.descriptors[descriptor])

