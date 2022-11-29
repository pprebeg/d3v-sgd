import itertools
from femdir.geofem import *
from timeit import default_timer as timer


class GeoGrillageFEM (GeoFEM):
    def __init__(self, name=''):
        """
        :param name:
        node_overlaps - Dictionary of nodes with expected overlaps.

        Conversion dictionaries for Grillage model properties [key] into GeoFEM
        properties [value]:

        plate_property_IDs - Model Plate properties into GeoFEM PlateProperty
        stiff_beam_prop_IDs - Model stiffener BeamProperty into GeoFEM Beam property
        half_stiff_beam_prop_IDs - Model stiffener BeamProperty with half original
            stiffeness for beams on Axis Of Symmetry into GeoFEM Beam property.
        web_prop_IDs - Model BeamProperty into GeoFEM PlateProperty
        flange_prop_IDs - Model BeamProperty into GeoFEM PlateProperty
        half_web_property_IDs - Model BeamProperty into GeoFEM PlateProperty,
            for Segments on AOS with half web thickness.
        """
        super().__init__(name)
        self.plate_elements = {}
        self.initial_node_overlaps = {}
        self.flange_element_overlaps = {}

        self.plate_property_IDs = {}
        self.stiff_beam_prop_IDs = {}
        self.half_stiff_beam_prop_IDs = {}
        self.web_property_IDs = {}
        self.flange_property_IDs = {}
        self.half_web_property_IDs = {}

        self.id_node_count = 1
        self.id_element_count = 1
        self.id_prop_count = 1

    def add_node(self, node_coords):
        """
        Add generated node to FEM model.
        :param node_coords:
        """
        node = Node(self.id_node_count, node_coords)
        self.addNode(node)
        self.id_node_count += 1
        return node

    def add_node_to_node_overlaps(self, node):
        """
        :param node:
        :return: Add edge nodes to FEM model overlaps dictionary.
        """
        self.initial_node_overlaps[node.id] = node

    def add_element_to_element_overlaps(self, element):
        """
        :param element:
        :return: Add flange elements to FEM model overlaps dictionary.
        """
        self.flange_element_overlaps[element.id] = element

    def add_to_plate_elements(self, element):
        """
        :param element:
        :return: Add plate elements to FEM model plate elements dictionary.
        """
        self.plate_elements[element.id] = element

    def add_element(self, elem: Element, idProp, nodeIds):
        """
        Add generated element to FEM model.
        :param elem:
        :param idProp: Element property ID.
        :param nodeIds: Node ID list.
        """
        elem.init(self.id_element_count)
        elem.property = self.getProperty(idProp)
        for idNod in nodeIds:
            node = self.getNode(idNod)
            elem.addNode(node)
        self.addElement(elem)
        self.id_element_count += 1
        return elem

    def add_quad_element(self, idProp, nodeIds):
        """
        Add generated quad element to FEM model.
        :param idProp:
        :param nodeIds:
        """
        elem = QuadElement()
        self.add_element(elem, idProp, nodeIds)
        return elem

    def add_tria_element(self, idProp, nodeIds):
        """
        Add generated triangle element to FEM model.
        :param idProp:
        :param nodeIds:
        """
        elem = TriaElement()
        self.add_element(elem, idProp, nodeIds)
        return elem

    def add_beam_element(self, idProp, nodeIds, vect_orient):
        """
        Add generated beam elemenet to FEM model.
        :param idProp:
        :param nodeIds:
        :param vect_orient:
        """
        elem = BeamElementShipStructure()
        dir_vector = BeamOrientationVector(vect_orient)
        elem.set_beam_orientation(dir_vector)
        self.add_element(elem, idProp, nodeIds)
        return elem

    def add_material(self, material_property):
        """
        :param material_property: Grillage model MaterialProperty object.
        :return: Add materials from grillage model.
            E - modulus of elasticity, [N/mm2]
            v - Poisson's ratio
            ro - material density, [kg/m3]
            Reh - yield strength, [N/mm2]
        """
        gfe_material = Material()
        gfe_material.init(material_property.id, material_property.name)
        gfe_material.E = material_property.E
        gfe_material.ni = material_property.v
        gfe_material.rho = material_property.ro
        gfe_material.ReH = material_property.Reh
        self.addMaterial(gfe_material)

    def add_property(self, prop):
        prop.id = self.id_prop_count
        self.addProperty(prop)
        self.id_prop_count += 1
        return prop

    def add_plate_property(self, prop_id, tp, mat):
        prop = PlateProperty()
        prop.init(id, 'Plate_property_' + str(prop_id))
        prop.tp = tp
        prop.material = mat
        self.add_property(prop)

    def add_T_beam_property(self, name, hw, tw, bf, tf, mat):
        prop = T_Profile_BeamProperty()
        prop.init(id, name)
        prop.hw = hw
        prop.tw = tw
        prop.bf = bf
        prop.tf = tf
        prop.material = mat
        self.add_property(prop)

    def add_half_T_beam_property(self, name, hw, tw, bf, tf, mat):
        prop = Half_T_Profile_BeamProperty()
        prop.init(id, name)
        prop.hw = hw
        prop.tw = tw
        prop.bf = bf
        prop.tf = tf
        prop.material = mat
        self.add_property(prop)

    def add_L_beam_property(self, name, hw, tw, bf, tf, mat):
        prop = L_Profile_BeamProperty()
        prop.init(id, name)
        prop.hw = hw
        prop.tw = tw
        prop.bf = bf
        prop.tf = tf
        prop.material = mat
        self.add_property(prop)

    def add_half_L_beam_property(self, name, hw, tw, bf, tf, mat):
        prop = Half_L_Profile_BeamProperty()
        prop.init(id, name)
        prop.hw = hw
        prop.tw = tw
        prop.bf = bf
        prop.tf = tf
        prop.material = mat
        self.add_property(prop)

    def add_FB_beam_property(self, name, hw, tw, mat):
        prop = FB_Profile_BeamProperty()
        prop.init(id, name)
        prop.hw = hw
        prop.tw = tw
        prop.material = mat
        self.add_property(prop)

    def add_half_FB_beam_property(self, name, hw, tw, mat):
        prop = Half_FB_Profile_BeamProperty()
        prop.init(id, name)
        prop.hw = hw
        prop.tw = tw
        prop.material = mat
        self.add_property(prop)

    def add_Hat_beam_property(self, name, h, t, bf, fi, mat):
        prop = Hat_Profile_BeamProperty()
        prop.init(id, name)
        prop.h = h
        prop.t = t
        prop.bf = bf
        prop.fi = fi
        prop.material = mat
        self.add_property(prop)

    def add_half_Hat_beam_property(self, name, h, t, bf, fi, mat):
        prop = Half_Hat_Profile_BeamProperty()
        prop.init(id, name)
        prop.h = h
        prop.t = t
        prop.bf = bf
        prop.fi = fi
        prop.material = mat
        self.add_property(prop)

    def add_Bulb_beam_property(self, name, hw_ekv, tw_ekv, bf_ekv, tf_ekv, mat):
        prop = Bulb_Profile_BeamProperty()
        prop.init(id, name)
        prop.hw_ekv = hw_ekv
        prop.tw_ekv = tw_ekv
        prop.bf_ekv = bf_ekv
        prop.tf_ekv = tf_ekv
        prop.material = mat
        self.add_property(prop)

    def add_half_Bulb_beam_property(self, name, hw_ekv, tw_ekv, bf_ekv, tf_ekv, mat):
        prop = Half_Bulb_Profile_BeamProperty()
        prop.init(id, name)
        prop.hw_ekv = hw_ekv
        prop.tw_ekv = tw_ekv
        prop.bf_ekv = bf_ekv
        prop.tf_ekv = tf_ekv
        prop.material = mat
        self.add_property(prop)

    @staticmethod
    def check_node_overlap(nodes_dict):
        """
        :return: Identifies coincident nodes in nodes_dict and returns
            overlap_list, where the first element are unique coordinates of
            overlapped nodes and the rest are overlapped geofementity Nodes.
        """
        overlap_list = []
        for node in nodes_dict.values():
            duplicate_coords = False
            for row in overlap_list:
                if np.allclose(node.p, row[0]):
                    row.append(node)
                    duplicate_coords = True
                    break
            if duplicate_coords is False:
                overlap_list.append([node.p, node])
        return overlap_list

    def check_node_overlap_np(self, nodes_dict):
        """
        :return: Identifies coincident nodes in nodes_dict and returns
            overlap_list, where the first element are unique coordinates of
            overlapped nodes and the rest are overlapped geofementity Node
            objects. Optimized version with NumPy functions.
        """
        print("Starting coincident node check...")
        start = timer()

        coords = [node.p for node in nodes_dict.values()]
        id_list = [node.id for node in nodes_dict.values()]

        coords_1 = np.expand_dims(coords, 0)
        coords_2 = np.expand_dims(coords, 1)
        # Relative node merge tolerance = 0.1mm
        boolean_array = np.isclose(coords_1, coords_2, rtol=1e-3).all(-1)
        boolean_array = np.tril(boolean_array)
        np.fill_diagonal(boolean_array, False)
        coincident_pairs = np.where(boolean_array)
        x_index, y_index = coincident_pairs
        n_overlaps = len(x_index)

        overlap_list = []
        for index in range(0, n_overlaps):
            node_1_id = id_list[int(x_index[index])]
            node_2_id = id_list[int(y_index[index])]
            node1 = self.nodes[node_1_id]
            node2 = self.nodes[node_2_id]

            duplicate_coords = False
            for row in overlap_list:
                if np.allclose(node1.p, row[0], rtol=1e-3):
                    if node1 not in row[1:]:
                        row.append(node1)
                    if node2 not in row[1:]:
                        row.append(node2)
                    duplicate_coords = True
                    break

            if duplicate_coords is False:
                overlap_list.append([node1.p, node1, node2])

        end = timer()
        print("Coincident node identification complete, found", len(overlap_list),
              "unique coordinates in", end - start, "s")

        return overlap_list

    def full_model_node_overlap_check(self):
        print("Starting full model coincident node check...")
        nodes_dict = self.nodes.values()
        coords = [node.p for node in nodes_dict]

        coords_1 = np.expand_dims(coords, 0)
        coords_2 = np.expand_dims(coords, 1)
        # Relative node merge tolerance = 0.1mm
        boolean_array = np.isclose(coords_1, coords_2, rtol=1e-3).all(-1)
        boolean_array = np.tril(boolean_array)
        np.fill_diagonal(boolean_array, False)
        x_index, y_index = np.where(boolean_array)
        n_overlaps = len(x_index)

        if n_overlaps > 1:
            print("Full model coincident node identification complete. "
                  "Node overlaps detected. Total overlaps:", n_overlaps)
        else:
            print("Full model coincident node identification complete."
                  " No node overlaps found.")

    @staticmethod
    def sorted_coincident_nodes(overlap_list):
        """
        Method sorts coincident nodes into remaining nodes and to delete nodes.

        :return: merge_nodes dictionary of nodes to be deleted (key) and
            replaced with (value), delete_list of all nodes to be deleted.
        """
        overlap_counter = 0
        merge_nodes = {}
        delete_list = []
        for row in overlap_list:
            for node_to_delete in row[2:]:
                merge_nodes[node_to_delete] = row[1]   # First node remains
                delete_list.append(node_to_delete)
            # print("Coordinates:", row[0], "nodes:",
            #       [node.id for node in row[1:]])
            overlap_counter += len(row[1:])

        print("Total coincident nodes:", overlap_counter)

        return merge_nodes, delete_list

    def merge_coincident_nodes(self):
        nodes_dict = self.initial_node_overlaps
        overlap_list = self.check_node_overlap_np(nodes_dict)
        merge_nodes, delete_list = self.sorted_coincident_nodes(overlap_list)

        start = timer()

        # Change overlapped nodes for all elements
        for element in self.elements.values():
            local_node_id = 0
            for node in element.nodes:
                if node in delete_list:
                    element.nodes[local_node_id] = merge_nodes[node]
                local_node_id += 1

        # Delete overlapped nodes
        for node in delete_list:
            del self.nodes[node.id]

        end = timer()
        print("Node merge complete, deleted", len(delete_list), "nodes in", end - start, "s")
        print("Total number of nodes:", self.num_nodes)
        print("Total number of elements before element merge:", self.num_elements)

    @staticmethod
    def check_element_overlap(element_dict):
        overlap_list = []
        element_combos = itertools.combinations(element_dict.values(), 2)
        for elements in element_combos:
            element_1, element_2 = elements
            nodes1 = [node.id for node in element_1.nodes]
            nodes2 = [node.id for node in element_2.nodes]

            if set(nodes1) == set(nodes2):
                overlap_list.append(elements)

        return overlap_list

    def merge_coincident_elements(self):
        element_dict = self.flange_element_overlaps
        overlapped_elements = self.check_element_overlap(element_dict)
        for elements in overlapped_elements:
            element_1, element_2 = elements

            tp_e1 = self.properties[element_1.prop_id].tp
            tp_e2 = self.properties[element_2.prop_id].tp

            if tp_e2 > tp_e1:
                del self.elements[element_1.id]
            else:
                del self.elements[element_2.id]

        print("Total number of elements after merge:", self.num_elements)

    def add_node_group(self, group_id: int, nodes: Dict):
        group = NodeGroup()
        group.init(group_id, "Node_group_" + str(group_id))
        for node in nodes.values():
            group.add_item(node)
        self.addGroup(group_id, group)

    def add_element_group(self, group_id: int, elements: Dict):
        group = ElementGroup()
        group.init(group_id, "Element_group_" + str(group_id))
        for element in elements.values():
            group.add_item(element)
        self.addGroup(group_id, group)

    def add_boundary_condition(self, bc_id: int, lc_id: int, values: List[float],
                               dof: List[int], set_id: int):
        """
        :param bc_id: Boundary condition ID
        :param lc_id: Load case ID
        :param values: Value for the DoF: 1 or 0; list of same length as dof
        :param dof: Degrees of Freedom: [wx=1, wy=2, wz=3, rx=4, ry=5, rz=6]
        :param set_id: GeoFEM nodal group ID
        """
        group = self.getGroup(set_id)
        bc = GroupNodalBC(bc_id, group, values, dof)
        self.addBoundaryCondition(bc)
        self.addBoundaryConditionToLoadcase(lc_id, bc)

    def add_pressure_load(self, load_id: int, lc_id: int, pressure: float,
                          set_id: int):
        """
        :param load_id: Pressure load ID
        :param lc_id: Load case ID
        :param pressure: External pressure
        :param set_id: GeoFEM element group ID
        :return:
        """
        group = self.getGroup(set_id)
        pressure_load = GroupPressureLoad(load_id, group, pressure, flip=True)
        self.addLoad(pressure_load)
        self.addLoadToLoadcase(lc_id, pressure_load)
