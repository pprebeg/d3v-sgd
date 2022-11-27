"""
Modul za testiranje diskretizacije učitanog modela
"""
from grillage.grillage_mesher import *
from timeit import default_timer as timer
# import femdir.geofementity as gfe
import os

# Get the current working directory
cwd = os.getcwd()
dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path)


def generate_test_mesh_v1():
    start = timer()

    hc_var = 1
    filename = str("../grillage savefiles/hc_var_") + str(hc_var) + str("_savefile.gin")
    hc_variant = GrillageModelData(filename).read_file()
    print("Testing FE mesh variant V1 for grillage variant", hc_var)

    gril_var = hc_variant
    aos_override = None
    ebs = 1     # Number of elements between stiffeners
    eweb = 3    # Number of elements along the height of PSM web
    eaf = 1     # Number of elements across primary supporting member flange
    far = 8     # Maximum PSM flange aspect ratio
    par = 4     # Maximum plate and PSM web aspect ratio
    dpar = 3    # Desired plating aspect ratio, less than the maximum

    grillage_mesh = MeshVariantV1(gril_var, aos_override)
    grill_fem = grillage_mesh.generate_grillage_mesh("test mesh", ebs, eweb, eaf, far, par, dpar)

    end = timer()
    print("Mesh generation time:", end - start, "s")
    # grillage_test_mesh.full_model_node_overlap_check()
    return grill_fem


def generate_test_mesh_v2():
    start = timer()

    hc_var = 5
    filename = str("../grillage savefiles/hc_var_") + str(hc_var) + str("_savefile.gin")
    hc_variant = GrillageModelData(filename).read_file()
    print("Testing FE mesh variant V2 for grillage variant", hc_var)

    gril_var = hc_variant
    aos_override = None
    ebs = 1     # Number of elements between stiffeners
    eweb = 3    # Number of elements along the height of PSM web
    eaf = 1     # Number of elements across primary supporting member flange
    far = 8     # Maximum PSM flange aspect ratio
    par = 4     # Maximum plate and PSM web aspect ratio
    dpar = 3    # Desired plating aspect ratio, less than the maximum

    grillage_mesh = MeshVariantV2(gril_var, aos_override)
    grill_fem = grillage_mesh.generate_grillage_mesh("test mesh", ebs, eweb, eaf, far, par, dpar)

    end = timer()
    print("Mesh generation time:", end - start, "s")
    # grillage_test_mesh.full_model_node_overlap_check()
    return grill_fem


def Test_edge_node_spacing():
    start = timer()

    hc_var = 5
    filename = str("../grillage savefiles/hc_var_") + str(hc_var) + str("_savefile.gin")
    hc_variant = GrillageModelData(filename).read_file()
    print("")
    print("Testing mesh dimensions for grillage variant", hc_var)
    aos_override = None
    ebs = 1     # Number of elements between stiffeners
    eweb = 3    # Number of elements along the height of PSM web
    eaf = 1     # Number of elements across primary supporting member flange
    far = 2.2     # Maximum PSM flange aspect ratio
    par = 1     # Maximum plate and PSM web aspect ratio
    dpar = 1    # Desired plating aspect ratio, less than the maximum
    mesh_extent = MeshExtent(hc_variant, aos_override)
    test_mesh_size = ElementSizeV2(mesh_extent, ebs, eweb, eaf, far, par, dpar)
    test_mesh_size.calculate_mesh_dimensions()

    for plate in test_mesh_size.mesh_extent.all_plating_zones.values():
        dim_x = test_mesh_size.plate_edge_node_spacing_x(plate)
        dim_y = test_mesh_size.plate_edge_node_spacing_y(plate)
        n_elem_x = len(dim_x)
        n_elem_y = len(dim_y)

        print("Plating zone ID", plate.id, ":")
        print("     Number of elements", n_elem_x, ", edge node spacing (dim_x):", dim_x)
        print("     Number of elements", n_elem_y, ", edge node spacing (dim_y):", dim_y, "\n")

    for segment in test_mesh_size.mesh_extent.all_segments.values():
        flange_spacing = test_mesh_size.flange_edge_node_spacing(segment)
        n_elem = len(flange_spacing)
        psm_id = segment.primary_supp_mem.id
        direct = segment.primary_supp_mem.direction.name
        psm_type = segment.beam_prop.beam_type.name
        fl_tr_dim = test_mesh_size.flange_transition_dim(segment)
        fl_tr_num = test_mesh_size.get_flange_transition_num(segment)

        print(direct, psm_type, "PSM, ID:", psm_id, ", segment ID:", segment.id,
              ", Number of elements:", n_elem, ", tr dim:", fl_tr_dim, ", tr num:", fl_tr_num, "total:", sum(flange_spacing.values()),
              " , flange edge node spacing:")
        print("   ", flange_spacing, "\n")
    end = timer()
    print("Test time:", end - start, "s")
