from femdir.geofem import *


class GeoGrillageFEM (GeoFEM):
    def __init__(self, name=''):
        super().__init__(name)
        self._node_overlaps = {}
        self._id_node_count = 1
        self._id_element_count = 1
        self._id_prop_count = 1

    def add_node(self, node_coords):
        node = Node(self._id_node_count, node_coords)
        self.addNode(node)
        self._id_node_count += 1
        return node

    def add_node_to_node_overlaps(self, node):
        self._node_overlaps[node.id] = node

    def add_element(self, elem: Element, idProp, nodeIds):
        if self._id_element_count == 55:
            print('brakepoint')
        elem.init(self._id_element_count)
        # elem.property = self.getProperty(idProp)
        for idNod in nodeIds:
            node = self.getNode(idNod)
            elem.addNode(node)
        self.addElement(elem)
        self._id_element_count += 1
        return elem

    def add_quad_element(self, idProp, nodeIds):
        elem = QuadElement()
        self.add_element(elem, idProp, nodeIds)
        return elem

    def add_beam_element(self, idProp, nodeIds, vect_orient):
        elem = BeamElement()
        elem.set_beam_orientation(BeamOrientationVector(vect_orient))
        self.add_element(elem, id, idProp, nodeIds)
        return elem

    def add_property(self, prop):
        prop.id = self._id_prop_count
        self.addProperty(prop)
        self._id_prop_count += 1
        return prop
