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
# Test_aos_stiffener(4)
# Test_aos_on_segment(BeamDirection.TRANSVERSE, 3, 1)


# ********************** testovi metoda u MeshSize **********************
# Test_identify_unique_property()
# Test_find_closest_divisor(4133, 935)
# Test_find_closest_divisor(4238, 935)
# Test_find_largest_divisor(4500, 1000)
# Test_element_aspect_ratio(700, 100)
# Test_refine_plate_element(4133, 830)
# Test_element_size_para_to_stiffeners(1)    # Dimenzije elementa oplate za odabranu zonu prema razmaku ukrepa
# Test_get_web_el_height(1, 1)
# Test_get_flange_el_length(BeamDirection.LONGITUDINAL, 1, 1)    # Dimenzije elementa prirubnice za odabrani segment
# Test_get_flange_el_width(1, 1)
# Test_get_min_flange_el_length()
# Test_get_min_flange_el_length_between_psm(1, 2)
# Test_element_size_plating_zone(1)         # Dimenzije elemenata na odabranoj zoni oplate prema razmaku ukrepa i ar prirubnice
# Test_ALL_element_size_plating_zone()        # dim_x i dim_y za sve zone oplate, prikaz u matrici
# Test_all_plating_zones_mesh_dimensions()    # Konačno usklađene i odabrane x i y dimenzije elemenata za sva polja oplate


# ********************** testovi metoda u MeshV1 **********************
# Test_element_size_mesh()                    # Konacno odabrane osnovne dimenzije mreze dim_x i dim_y
# Test_get_reduced_plate_dim(1)
# Test_element_size_transition(1, 2)
# Test_get_tr_dim_x(5)
# Test_get_tr_dim_y(4)
# Test_get_tr_element_num(1)
# Test_get_element_number(2)
# Test_all_element_numbers()
# Test_get_mesh_dim(2)                        # Dimenzije x i y svih elemenata duž odabrane zone oplate
# Test_get_all_mesh_dim()                     # Dimenzije x i y svih elemenata duž svih generiranih zona oplate
# Test_calc_element_base_size()


# ********************** testovi metoda u PlatingZoneMesh **********************
# Test_PlatingZoneMesh(1, AOS.NONE)                                         # Izrada mreže jedne zone oplate


# ********************** testovi metoda u SegmentMesh **********************
# Test_get_segment_web_element_property(BeamDirection.LONGITUDINAL, 1, 1)
# Test_Segment_element_generation(BeamDirection.LONGITUDINAL, 1, 1)
# Test_edge_segment_node_generation(BeamDirection.LONGITUDINAL, 1, 1)


# ********************** testovi metoda u GrillageMesh **********************
# GrillageMesh(mesh1).generate_plate_mesh()         # Izrada mreže svih zona oplate
# GrillageMesh(mesh1).generate_psm_mesh()           # Izrada mreže svih segmenata - Baca grešku, neki zajeb s prenašanjem element id
# GrillageMesh(mesh1).generate_mesh()


end = timer()

print("***************************************************************************************************************")
print("Code run time =", end - start, "s")
