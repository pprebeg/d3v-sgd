import itertools
from femdir.geofem import *


class GeoGrillageFEM (GeoFEM):
    def __init__(self, name=''):
        """
        :param name:
        node_overlaps - Dictionary of nodes with expected overlaps.

        Conversion dictionaries for Grillage model properties [key] into GeoFEM
        properties [value]:

        plate_property_IDs - Model Plate properties into GeoFEM PlateProperty
        stiff_beam_prop_IDs - Model stiffener BeamProperty into GeoFEM Beam property
        web_prop_IDs - Model BeamProperty into GeoFEM PlateProperty
        flange_prop_IDs - Model BeamProperty into GeoFEM PlateProperty
        half_web_property_IDs - Model BeamProperty into GeoFEM PlateProperty,
            for Segments on AOS with half web thickness.
        """
        super().__init__(name)
        self.initial_node_overlaps = {}

        self.plate_property_IDs = {}
        self.stiff_beam_prop_IDs = {}
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
        elem = BeamElement()
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

    def add_FB_beam_property(self, name, hw, tw, mat):
        prop = FB_Profile_BeamProperty()
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

    def add_Bulb_beam_property(self, name, hw_ekv, tw_ekv, bf_ekv, tf_ekv, mat):
        prop = Bulb_Profile_BeamProperty()
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
        for node in self.initial_node_overlaps.values():
            duplicate_coords = False
            if not overlap_array:
                overlap_array.append([node.p, node])
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
                overlap_counter += len(row[1:])
                remaining_node = row[1]         # First node remains
                to_delete_nodes = row[2:]       # Duplicate nodes
                coincident_nodes.append([remaining_node, to_delete_nodes])

        print("Coincident node identification complete. "
              "Total coincident nodes:", overlap_counter,
              ", at", len(overlap_array), "unique coordinates.")

        return coincident_nodes

    def full_model_node_overlap_check(self):
        all_overlaps_dict = {}  # Pairs of overlapping nodes
        all_nodes = self.nodes
        print("Starting full model coincident node check...")

        node_combinations = itertools.combinations(all_nodes.keys(), 2)
        for nodes in node_combinations:
            node1_id, node2_id = nodes
            if np.allclose(all_nodes[node1_id].p, all_nodes[node2_id].p):
                all_overlaps_dict[node1_id] = all_nodes[node2_id]
        n_overlaps = len(all_overlaps_dict)

        if n_overlaps > 1:
            print("Full model coincident node identification complete. "
                  "Node overlaps detected. Total overlaps:", n_overlaps)
        else:
            print("Full model coincident node identification complete."
                  " No node overlaps found.")

    def merge_coincident_nodes(self):
        coincident_nodes = self.check_node_overlap()
        print("Starting coincident node merge...")

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
        print("Total number of elements:", self.num_elements)

        # for node in self.initial_node_overlaps.values():
        #     print(node.id)

        # TEST PRETRAGE
        """
        for item in merge_nodes.items():
            key, val = item
            print("zamjena čvora", key.id, "sa", val.id)
        for node in delete_list:
            print("Brisanje čvora", node.id)
        """
        # TEST Pretraga za brisanje
        """
        for node in self.nodes.items():
            n_id, n_object = node
            if n_object in delete_list:
                print("Brisanje cvora", n_object.id, n_id, self.nodes[n_id].id)
        """
        # TEST Nakon merge:
        """
        for element in self.elements.items():
            el_id, el_object = element
            element_type = el_object.get_type()
            print(el_id, element_type.name, ", element nodes:", [node.id for node in el_object.nodes])
        """

    def check_element_overlap(self):
        delete_list = []
        element_combos = itertools.combinations(self.elements.keys(), 2)
        for elements in element_combos:
            element_1_id, element_2_id = elements
            element_1 = self.elements[element_1_id]
            element_2 = self.elements[element_2_id]

            nodes1 = [node.id for node in element_1.nodes]
            nodes2 = [node.id for node in element_2.nodes]

            e1_type = element_1.get_type().name
            e2_type = element_2.get_type().name
            e1_prop_id = element_1.prop_id
            e2_prop_id = element_2.prop_id

            if set(nodes1) == set(nodes2):
                print("Overlapped elements:", e1_type, element_1.id,
                      e2_type, element_2.id,
                      ", on nodes:", nodes1, nodes2)

                tp_e1 = self.properties[e1_prop_id].tp
                tp_e2 = self.properties[e2_prop_id].tp

                if tp_e2 > tp_e1:
                    # element_2 remains
                    # delete element_1
                    # delete_list.append(element_to_delete)
                    pass
                else:
                    # element_1 remains
                    # delete element_2
                    # delete_list.append(element_to_delete)
                    pass

                print("   ", e1_type, element_1.id, ":", tp_e1, "mm, ",
                      e2_type, element_2.id, ":", tp_e2, "mm")

        return delete_list

    def merge_overlapped_elements(self):
        # overlapped_elements = self.check_element_overlap()
        # delete all in delete list...
        pass
