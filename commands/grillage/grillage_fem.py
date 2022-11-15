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

    def check_node_overlap(self):
        """
        :return: Identifies coincident nodes in initial_node_overlaps and
            returns coincident_nodes array, where first element contains ID of
            the node that will remain at the coordinates. The rest are IDs of
            overlapped nodes that will be deleted.
        """
        print("Starting coincident node check...")
        overlap_array = []
        # nodes_dict = self.nodes.values()                  # Značajno sporije!
        nodes_dict = self.initial_node_overlaps.values()

        for node in nodes_dict:
            duplicate_coords = False
            for row in overlap_array:
                if np.allclose(node.p, row[0]):
                    row.append(node)
                    duplicate_coords = True
                    break
            if duplicate_coords is False:
                overlap_array.append([node.p, node])

        overlap_counter = 0
        coincident_nodes = []
        for row in overlap_array:
            if len(row[1:]) > 1:
                # print("Coordinates:", row[0], "nodes:", [node.id for node in row[1:]])
                overlap_counter += len(row[1:])
                remaining_node = row[1]         # First node remains
                to_delete_nodes = row[2:]       # Duplicate nodes
                coincident_nodes.append([remaining_node, to_delete_nodes])

        print("Coincident node identification complete. "
              "Total coincident nodes:", overlap_counter,
              ", at", len(overlap_array), "unique coordinates.")
        return coincident_nodes

    def check_node_overlap_np(self):
        """
        :return: Identifies coincident nodes in initial_node_overlaps and
            returns coincident_nodes array, where first element contains ID of
            the node that will remain at the coordinates. The rest are IDs of
            overlapped nodes that will be deleted. Uses NumPy functions to
            identify node overlaps instead of for loops.
        """
        print("Starting coincident node check...")
        # nodes_dict = self.nodes.values()                  # Nije puno sporije!
        nodes_dict = self.initial_node_overlaps.values()
        coords = [node.p for node in nodes_dict]
        id_list = [node.id for node in nodes_dict]

        coords_1 = np.expand_dims(coords, 0)
        coords_2 = np.expand_dims(coords, 1)
        boolean_array = np.isclose(coords_1, coords_2).all(-1)
        boolean_array = np.tril(boolean_array)
        np.fill_diagonal(boolean_array, False)
        coincident_pairs = np.where(boolean_array)
        x_index, y_index = coincident_pairs
        n_overlaps = len(x_index)

        start = timer()

        overlap_array = []
        for index in range(0, n_overlaps):
            node_1_id = id_list[int(x_index[index])]
            node_2_id = id_list[int(y_index[index])]
            node1 = self.nodes[node_1_id]
            node2 = self.nodes[node_2_id]

            duplicate_coords = False
            for row in overlap_array:
                if np.allclose(node1.p, row[0]):
                    if node1 not in row[1:]:
                        row.append(node1)
                    if node2 not in row[1:]:
                        row.append(node2)
                    duplicate_coords = True
                    break

            if duplicate_coords is False:
                overlap_array.append([node1.p, node1, node2])

        end = timer()
        print("Sort coincident nodes time:", end - start, "s")

        overlap_counter = 0
        coincident_nodes = []
        for row in overlap_array:
            if len(row[1:]) > 1:
                # print("Coordinates:", row[0], "nodes:", [node.id for node in row[1:]])
                overlap_counter += len(row[1:])
                remaining_node = row[1]  # First node remains
                to_delete_nodes = row[2:]  # Duplicate nodes
                coincident_nodes.append([remaining_node, to_delete_nodes])

        print("Coincident node identification complete. "
              "Total coincident nodes:", overlap_counter,
              ", at", len(overlap_array), "unique coordinates.")
        return coincident_nodes

    def full_model_node_overlap_check(self):
        print("Starting full model coincident node check...")
        nodes_dict = self.nodes.values()

        x_coords = [node.p[0] for node in nodes_dict]
        y_coords = [node.p[1] for node in nodes_dict]
        z_coords = [node.p[2] for node in nodes_dict]

        x_boolean_array = np.isclose(x_coords, np.vstack(x_coords))
        y_boolean_array = np.isclose(y_coords, np.vstack(y_coords))
        z_boolean_array = np.isclose(z_coords, np.vstack(z_coords))

        np.fill_diagonal(x_boolean_array, False)
        np.fill_diagonal(y_boolean_array, False)
        np.fill_diagonal(z_boolean_array, False)

        x_tri_array = np.tril(x_boolean_array)
        y_tri_array = np.tril(y_boolean_array)
        z_tri_array = np.tril(z_boolean_array)

        xy_boolean_array = np.logical_and(x_tri_array, y_tri_array)
        xz_boolean_array = np.logical_and(x_tri_array, z_tri_array)
        boolean_array = np.logical_and(xy_boolean_array, xz_boolean_array)

        find_where = np.where(boolean_array)    # Indexes of overlapped nodes
        x_index, y_index = find_where
        n_overlaps = len(x_index)

        unique_list = []
        # n_stack = np.stack((x_index, y_index), axis=1)  # Sortirani parovi preklopljenih čvorova
        # print(n_stack)

        # vals, counts = np.unique(n_stack, return_counts=True)
        # print(vals)
        # print(counts)

        # vals, index, counts = np.unique(n_stack, return_index=True, return_counts=True)
        # print(vals)
        # print(index)
        # print(counts)

        # inters = np.intersect1d(x_index, y_index)   # Svi čvorovi koji se preklapaju
        # print(inters)

        # masked_arr = np.ma.masked_where(counts <= 1, counts)    # čvorovi koji se pojavljuju samo jednom
        # print(masked_arr)

        if n_overlaps > 1:
            print("Full model coincident node identification complete. "
                  "Node overlaps detected. Total overlaps:", n_overlaps)
        else:
            print("Full model coincident node identification complete."
                  " No node overlaps found.")

    def merge_coincident_nodes(self):
        coincident_nodes = self.check_node_overlap_np()      # New NumPy check node overlap

        merge_nodes = {}
        delete_list = []
        for nodes in coincident_nodes:
            for node_to_delete in nodes[1]:
                merge_nodes[node_to_delete] = nodes[0]
                delete_list.append(node_to_delete)

        # Change overlapped nodes for all elements
        for element in self.elements.items():
            el_id, el_object = element
            local_node_id = 0
            for node in el_object.nodes:
                if node in delete_list:
                    # print("         ", node.id, ", zamjena sa", merge_nodes[node].id, ", indeks u listi čvorova", local_node_id)
                    el_object.nodes[local_node_id] = merge_nodes[node]
                local_node_id += 1

        # Delete overlapped nodes
        for delete_node in delete_list:
            del self.nodes[delete_node.id]

        print("Node merge complete, deleted", len(delete_list), "nodes.")
        print("Total number of nodes:", self.num_nodes)
        print("Total number of elements before element merge:", self.num_elements)

        # for node in self.initial_node_overlaps.values():
        #     print(node.id)

    def check_element_overlap(self):
        overlap_array = []
        element_dict = self.flange_element_overlaps
        element_combos = itertools.combinations(element_dict.keys(), 2)
        for elements in element_combos:
            element_1_id, element_2_id = elements
            element_1 = element_dict[element_1_id]
            element_2 = element_dict[element_2_id]

            nodes1 = [node.id for node in element_1.nodes]
            nodes2 = [node.id for node in element_2.nodes]

            if set(nodes1) == set(nodes2):
                overlap_array.append(elements)

        return overlap_array

    def merge_coincident_elements(self):
        overlapped_elements = self.check_element_overlap()
        delete_list = []
        for elements in overlapped_elements:
            element_1_id, element_2_id = elements
            element_1 = self.elements[element_1_id]
            element_2 = self.elements[element_2_id]
            tp_e1 = self.properties[element_1.prop_id].tp
            tp_e2 = self.properties[element_2.prop_id].tp

            if tp_e2 > tp_e1:
                delete_list.append(element_1)
            else:
                delete_list.append(element_2)

        # Delete overlapped elements
        for delete_element in delete_list:
            del self.elements[delete_element.id]

        print("Total number of elements:", self.num_elements)

    def identify_plating_nodes(self):
        """
        :return: Dictionary of all plating nodes for pressure load case.
        """
        plating_nodes = {}
        nodes_dict = self.nodes.values()
        z_coords = [node.p[2] for node in nodes_dict]
        id_list = [node.id for node in nodes_dict]

        hw = np.max(z_coords)
        boolean_array = np.isclose(z_coords, hw)
        node_index = np.where(boolean_array)
        node_index = np.concatenate(node_index)

        for index in range(0, len(node_index)):
            node_id = id_list[node_index[index]]
            node = self.nodes[node_id]
            plating_nodes[node.id] = node
        return plating_nodes
