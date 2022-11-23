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

    hc_var = 4
    filename = str("../grillage savefiles/hc_var_") + str(hc_var) + str("_savefile.gin")
    hc_variant = GrillageModelData(filename).read_file()
    print("Testing FE mesh variant V1 for grillage variant", hc_var)

    mesh_var = MeshVariant.V1
    gril_var = hc_variant
    aos_override = None
    grillage_mesh = GrillageMesh(mesh_var, gril_var, aos_override)
    grill_fem = grillage_mesh.generate_grillage_mesh("test mesh", ebs=1, eweb=3, eaf=1, far=8, par=4, dpar=3)

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

    mesh_var = MeshVariant.V2
    gril_var = hc_variant
    aos_override = None
    grillage_mesh = GrillageMesh(mesh_var, gril_var, aos_override)
    grill_fem = grillage_mesh.generate_grillage_mesh("test mesh", ebs=1, eweb=3, eaf=1, far=8, par=4, dpar=3)

    end = timer()
    print("Mesh generation time:", end - start, "s")
    # grillage_test_mesh.full_model_node_overlap_check()
    return grill_fem
