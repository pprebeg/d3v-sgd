"""
Modul za testiranje diskretizacije učitanog modela - Pokretanje testova
"""
from tests.test_mesh import *

# ********************** testovi metoda u ModelCheck **********************
# Test_Model_check()
# Test_mesh_feasibility()
# print(ModelCheck(hc_variant).assign_symmetry())   # Konačno odabrana os simetrije s kojom se ide u izradu mreže


# ********************** testovi metoda u MeshExtent **********************
# Test_plating_zones_ref_array()
# Test_get_plate_dim(2)
# Test_psm_extent()
# Test_segment_extent()
# Test_grillage_plate_extent()
# Test_aos_stiffener(2)
# Test_aos_on_segment(BeamDirection.LONGITUDINAL, 3, 1)
# Test_identify_split_element_zones()


# ********************** testovi metoda u MeshSize **********************
# Test_identify_unique_property()
# Test_find_closest_divisor(4133, 935)
# Test_find_closest_divisor(4238, 935)
# Test_find_largest_divisor(4500, 1000)
# Test_element_aspect_ratio(700, 100)
# Test_refine_plate_element(4938.3, 1100)
# Test_element_size_para_to_stiffeners(1)    # Dimenzije elementa oplate za odabranu zonu prema razmaku ukrepa
# Test_get_web_el_height(1, 1)
# Test_get_flange_el_length(BeamDirection.LONGITUDINAL, 1, 1)    # Dimenzije elementa prirubnice za odabrani segment
# Test_get_flange_el_width(1, 1)
# Test_get_min_flange_el_length()
# Test_get_min_flange_el_length_between_psm(1, 2)
# Test_element_size_plating_zone(1)         # Dimenzije elemenata na odabranoj zoni oplate prema razmaku ukrepa i ar prirubnice
# Test_ALL_element_size_plating_zone()        # dim_x i dim_y za sve zone oplate, prikaz u matrici
# Test_all_plating_zones_mesh_dimensions()    # Konačno usklađene i odabrane x i y dimenzije elemenata za sva polja oplate
# Test_saved_plate_edge_node_dims()
# Test_get_split_elements_number(5)


# ********************** testovi metoda u MeshV1 **********************
# Test_element_size_mesh()                    # Konacno odabrane osnovne dimenzije mreze dim_x i dim_y
# Test_get_reduced_plate_dim(1)
# Test_transition_element_size_plating_zone(1, 1)
# Test_get_tr_dim_x(1)
# Test_get_tr_dim_y(4)
# Test_get_tr_element_num(1)
# Test_get_element_number(2)

# Test_get_mesh_dim(5)                        # Dimenzije x i y svih elemenata duž odabrane zone oplate
# Test_get_all_mesh_dim()                     # Dimenzije x i y svih elemenata duž svih generiranih zona oplate
# Test_calc_element_base_size()


# ********************** testovi metoda u PlatingZoneMesh **********************
# Test_PlatingZoneMesh(1, AOS.NONE)                                         # Izrada mreže jedne zone oplate
# Test_PlatingZoneMesh_element_property(1)
# Test_PlatingZoneMesh_beam_elements(5, AOS.BOTH)


# ********************** testovi metoda u SegmentMesh **********************
# Test_edge_segment_node_generation(BeamDirection.LONGITUDINAL, 1, 1)
# Test_get_web_element_property(BeamDirection.LONGITUDINAL, 3, 2)
# Test_get_flange_element_property(BeamDirection.TRANSVERSE, 2, 1)
# Test_Segment_element_generation(BeamDirection.LONGITUDINAL, 2, 1)                 # Izrada mreže na jednom segmentu
# Test_generate_inward_flange_nodes(BeamDirection.LONGITUDINAsL, 1, 1, FlangeDirection.INWARD)
# Test_flange_ref_array(BeamDirection.LONGITUDINAL, 2, 1)


# ********************** testovi metoda u SegmentV2 **********************
# Test_generate_element_row(BeamDirection.LONGITUDINAL, 1, 1)


# ********************** testovi metoda u GrillageMesh **********************
# GrillageMesh(mesh1).generate_plate_mesh()         # Izrada mreže svih zona oplate
# GrillageMesh(mesh1).generate_psm_mesh_V1()           # Izrada mreže svih segmenata
# GrillageMesh(mesh1).generate_mesh()
# GrillageMesh(mesh1).check_node_overlap()


# Test nove ideje za algoritam pretraživanja preklapanja čvorova
"""
testset = [6245, 525, 1006]
testset2 = [624, 52, 106]

arraytest = [[testset, 1], [testset, 3, 4]]

print("Čvorovi na drugom setu koordinata:", arraytest[1][1:])
print(arraytest)
arraytest.append([testset2, 76, 65, 54])
print("Nakon dodavanja novog seta koordinata i čvorova:", arraytest)
arraytest[2].append(888)
print(arraytest)
print("samo koordinate testset2:", arraytest[0][0])
print("FOR TEST ************************")
print("\n")
test_coords = [624, 52, 106]
test_id = 999

exist = False
for coords in arraytest:
    print("koordinate:", coords[0], "čvorovi:", coords[1:], "cijeli redak:", coords)
    if np.allclose(test_coords, coords[0]):
        print("provjera:", test_coords, coords[0])
        coords.append(test_id)
        print("Test koordinate postoje u polju!")
        exist = True
        break

if exist is False:
    print("Test koordinate NE postoje u polju")
    arraytest.append([test_coords, test_id])
print("Polje nakon petlje:", arraytest)
"""

end = timer()

print("***************************************************************************************************************")
print("Code run time =", end - start, "s")
