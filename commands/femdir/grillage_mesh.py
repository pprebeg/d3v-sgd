"""
Tool for Grillage Structure Analysis
University of Zagreb, Faculty of Mechanical Engineering and Naval Architecture
Department of Naval Architecture and Ocean Engineering

Master's thesis project

    Gordan Kos, univ.bacc.ing.nav.arch.
    Dr.sc. Pero Prebeg, dipl.ing.


MODULE FOR GRILLAGE FINITE ELEMENT MESH DEFINITION

Assumptions:
    1.) All stiffeners are continuous over the breadth and length of the hatch cover,
        in compliance with IACS CSR: Chapter 9, Section 5, 2. Arrangements, 2.2.2
    2.) Web height of all primary supporting members is the same for simplicity

"""

from grillage.grillage_model import *

# import femdir.geofementity as gfe
# import femdir.geofem as gfem

# from femdir.geofem import *
from femdir.geofementity import *


# class GrillageMesh(gfem.GeoFEM):
# class GrillageMesh(GeoFEM):
class GrillageMesh:
    def __init__(self, grillage: Grillage):
        # super().__init__()
        self._grillage = grillage
        self._plating_mesh_dim = {}
        self._mesh_dim_x = []               # List of mesh x dimensions in the longitudinal direction
        self._mesh_dim_y = []               # List of mesh y dimensions in the transverse direction
        self._plating_nodes = []        # Preimenovati da se zna da sadrži samo referentni ID grupe čvorova oplate *** WIP
        self._flange_aspect_ratio = 4.0     # Maximum aspect ratio value for primary supporting member flange quad elements
        self._plate_aspect_ratio = 4.0      # Maximum aspect ratio value for plating and primary supporting member web quad elements

    def add_element(self, el_id, el_property, nodes):
        pass

    def plating_mesh_dim(self):
        return self._plating_mesh_dim

    def add_plating_mesh_dim(self, plate_id, value):
        self._plating_mesh_dim[plate_id] = value

    @property
    def mesh_dim_x(self):
        return self._mesh_dim_x

    @property
    def mesh_dim_y(self):
        return self._mesh_dim_y

    def plating_nodes(self):
        return self._plating_nodes

    @property
    def flange_aspect_ratio(self):
        return self._flange_aspect_ratio

    @flange_aspect_ratio.setter
    def flange_aspect_ratio(self, value):
        self._flange_aspect_ratio = value

    @property
    def plate_aspect_ratio(self):
        return self._plate_aspect_ratio

    @plate_aspect_ratio.setter
    def plate_aspect_ratio(self, value):
        self._plate_aspect_ratio = value

    @staticmethod
    def find_closest_divisor(length, value):
        # Method for determining the number of Quad elements along length L, when one dimension (value) of the element is known
        """
        :param length: Length L which should be divided into n equal parts, each with length x.
        :param value: Value to which length x should be closes to.
        :return: Closest divisor of length L, which results in a length x closest to given value.
        Equivalent to the number of Finite Elements n, with dimension x.
        """
        if np.mod(length, value) == 0:
            return length / value
        else:
            i = 1
            res = []
            while i <= length:
                if np.mod(length, i) == 0:
                    res.append(i)
                i += 1

            min_diff = value
            min_div_id = 1
            if not res:     # If input dimensions are decimal, a divisor may not exist
                return np.round(length / value, 0)
            else:           # If input dimensions are integers and a divisor exists
                for i in range(0, len(res)):
                    if min_diff > abs((length / res[i]) - value):
                        min_diff = abs((length / res[i]) - value)
                        min_div_id = i
                n = res[min_div_id]
                return n

    # Metoda za identifikaciju čvora preko koordinata *** WIP
    """
    def IdentifyNode(self, x, y, z):
        # Returns FENode ID located at coordinates x, y ,z
        for i in FEMesh.nodes(self).keys():
            coords = FEMesh.nodes(self)[i].coords
            if coords[0] == x and coords[1] == y and coords[2] == z:
                return FEMesh.nodes(self)[i].id
    """

    def CheckNodeOverlap(self):
        pass

    def element_size_stiffener_spacing(self, plate: Plate):
        # Returns the quad element size, based only on stiffener spacing and assuming one element between stiffeners for a plating zone
        stiff_spacing = Plate.get_stiffener_spacing(plate) * 1000   # Stiffener spacing in [mm]

        if plate.stiff_dir == BeamDirection.LONGITUDINAL:
            L = Plate.plate_longitudinal_dim(plate) * 1000          # Longitudinal plating zone dimension [mm]
            dim_y = stiff_spacing                                   # Distance between nodes in the transverse (y) direction
            dim_x = L / self.find_closest_divisor(L, dim_y)         # Distance between nodes in the longitudinal (x) direction
            return np.array([dim_x, dim_y])

        elif plate.stiff_dir == BeamDirection.TRANSVERSE:
            B = Plate.plate_transverse_dim(plate) * 1000            # Transverse plating zone dimension [mm]
            dim_x = stiff_spacing                                   # Distance between nodes in the longitudinal (x) direction
            dim_y = B / self.find_closest_divisor(B, dim_x)         # Distance between nodes in the transverse (y) direction
            return np.array([dim_x, dim_y])

    def element_size_flange_width(self, grillage, segment: Segment):
        # Returns the maximum quad element length based on flange width and maximum allowed aspect ratio for a segment
        max_aspect_ratio = self._flange_aspect_ratio
        dim_max = 0.0
        beam_property = segment.beam_prop
        bf_net = segment.beam_prop.bf - grillage.corrosion_addition()[1].tc

        if isinstance(beam_property, TBeamProperty) and not \
                (isinstance(beam_property, FBBeamProperty) or isinstance(beam_property, LBeamProperty)):
            dim_max = max_aspect_ratio * 0.5 * bf_net

        elif isinstance(beam_property, LBeamProperty):
            dim_max = max_aspect_ratio * bf_net

        return dim_max

    def element_size_plating_zone(self, grillage, plate: Plate):
        # Returns the quad element size based on stiffener spacing and maximum allowed aspect ratio for a plating zone
        plate_aspect_ratio = self._plate_aspect_ratio
        dim_p = self.element_size_stiffener_spacing(plate)
        dim_xp = dim_p[0]
        dim_yp = dim_p[1]

        dim_xf1 = self.element_size_flange_width(grillage, plate.long_seg1)
        dim_xf2 = self.element_size_flange_width(grillage, plate.long_seg2)
        dim_yf1 = self.element_size_flange_width(grillage, plate.trans_seg1)
        dim_yf2 = self.element_size_flange_width(grillage, plate.trans_seg2)

        dim_xf = np.minimum(dim_xf1, dim_xf2)   # Minimum element x dimension based on flange element aspect ratio
        dim_yf = np.minimum(dim_yf1, dim_yf2)   # Minimum element y dimension based on flange element aspect ratio

        dim_x = dim_xp  # Initial element size along x axis is based on stiffener spacing
        dim_y = dim_yp  # Initial element size along y axis is based on stiffener spacing

        if plate.stiff_dir == BeamDirection.TRANSVERSE:
            if dim_x > dim_xf:                         # If element size based on stiffener spacing along x asis is larger than flange aspect ratio allows for
                div_round_up = np.ceil(dim_x / dim_xf)       # Equal division of elements between transverse stiffeners
                dim_x = dim_x / div_round_up                    # If ordinary stiffeners are oriented transversely, x dimension has to be divided into equal parts

                if dim_y / dim_x > plate_aspect_ratio:    # Checks plating element aspect ratio after reducing dim_x for transverse stiffeners
                    div_round_up = np.ceil(dim_y / dim_x)
                    dim_y = dim_y / div_round_up

            if dim_y > dim_yf:                      # Checks flange element aspect ratio for transverse stiffeners
                div_round_up = np.ceil(dim_y / dim_yf)
                dim_y = dim_y / div_round_up

        if plate.stiff_dir == BeamDirection.LONGITUDINAL:
            if dim_y > dim_yf:                         # If element size based on stiffener spacing along y asis is larger than flange aspect ratio allows for
                div_round_up = np.ceil(dim_y / dim_yf)       # Equal division of elements between longitudinal stiffeners
                dim_y = dim_y / div_round_up                    # If ordinary stiffeners are oriented longitudinally, y dimension has to be divided into equal parts

                if dim_x / dim_y > plate_aspect_ratio:    # Checks plating element aspect ratio after reducing dim_y for longitudinal stiffeners
                    div_round_up = np.ceil(dim_x / dim_y)
                    dim_x = dim_x / div_round_up

            if dim_x > dim_xf:                      # Checks flange element aspect ratio for longitudinal stiffeners
                div_round_up = np.ceil(dim_x / dim_xf)
                dim_x = dim_x / div_round_up

        return np.array((dim_x, dim_y))

    def element_size_mesh(self, grillage):
        if grillage.hc_variant_check() is True:
            n_long = int(grillage.N_transverse - 1)     # Number of plating zones along the longitudinal axis
            n_tran = int(grillage.N_longitudinal - 1)   # Number of plating zones along the transverse axis

            self._mesh_dim_x = np.zeros(n_long)
            self._mesh_dim_y = np.zeros(n_tran)

            plating_mesh_dim_x = {}        # Dictionary of calculated dimension x for all plating zones
            plating_mesh_dim_y = {}        # Dictionary of calculated dimension y for all plating zones

            # Calculate the quad element size based on stiffener spacing and maximum allowed aspect ratio for all plating zones
            #   Dimensions are saved into dictionaries plating_mesh_dim_x, plating_mesh_dim_y
            for plate in grillage.plating().values():
                plate_zones_dim = self.element_size_plating_zone(grillage, plate)
                dim_x = plate_zones_dim[0]
                dim_y = plate_zones_dim[1]

                plating_mesh_dim_x[plate.id] = dim_x
                plating_mesh_dim_y[plate.id] = dim_y

            # Assign dimension y for all plating zones between longitudinal primary supporting members
            for i_long in range(1, len(grillage.longitudinal_members())):
                long1 = grillage.longitudinal_members()[i_long]
                long2 = grillage.longitudinal_members()[i_long + 1]
                plating_zones = grillage.plating_zones_between_psm(long1, long2)    # List of all plating zones between PSM

                restriction_y = False
                dim_y_list = []                     # List of element y dimensions for all plates between PSM

                min_y = 0.0
                for i in range(0, len(plating_zones)):
                    # Minimum y dimension from flange width
                    dim_yf1 = self.element_size_flange_width(grillage, plating_zones[i].trans_seg1)
                    dim_yf2 = self.element_size_flange_width(grillage, plating_zones[i].trans_seg2)
                    min_y = np.minimum(dim_yf1, dim_yf2)

                for i in range(0, len(plating_zones)):
                    plate_id = plating_zones[i].id
                    stiff_dir = plating_zones[i].stiff_dir  # Stiffener direction of plate in list plating_zones
                    dim_y = plating_mesh_dim_y[plate_id]    # Quad element size in the y direction for plate in list plating_zones
                    dim_y_list.append(dim_y)

                    if stiff_dir == BeamDirection.LONGITUDINAL:    # If dimension y is limited by longitudinal stiffener spacing
                        restriction_y = True        # Dimension restriction along y axis exists because there are longitudinal stiffeners
                        self._mesh_dim_y[i_long - 1] = dim_y

                        if dim_y > min_y:                           # If dimension y exceeds the minimum required
                            div_round_up = np.ceil(dim_y / min_y)
                            dim_y = dim_y / div_round_up
                            self._mesh_dim_y[i_long - 1] = dim_y

                if restriction_y is False:          # If there are no longitudinal stiffeners between long1 and long2
                    dim_y = np.amin(dim_y_list)
                    self._mesh_dim_y[i_long - 1] = dim_y

            # Assign dimension x for all plating zones between transverse primary supporting members
            for i_tran in range(1, len(grillage.transverse_members())):
                tran1 = grillage.transverse_members()[i_tran]
                tran2 = grillage.transverse_members()[i_tran + 1]
                plating_zones = grillage.plating_zones_between_psm(tran1, tran2)  # List of all plating zones between PSM

                restriction_x = False
                dim_x_list = []                     # List of element x dimensions for all plates between PSM

                min_x = 0.0
                for i in range(0, len(plating_zones)):
                    # Minimum x dimension from flange width
                    dim_xf1 = self.element_size_flange_width(grillage, plating_zones[i].long_seg1)
                    dim_xf2 = self.element_size_flange_width(grillage, plating_zones[i].long_seg2)
                    min_x = np.minimum(dim_xf1, dim_xf2)

                for i in range(0, len(plating_zones)):
                    plate_id = plating_zones[i].id
                    stiff_dir = plating_zones[i].stiff_dir
                    dim_x = plating_mesh_dim_x[plate_id]
                    dim_x_list.append(dim_x)

                    if stiff_dir == BeamDirection.TRANSVERSE:      # If dimension x is limited by transverse stiffener spacing
                        restriction_x = True        # Dimension restriction along x axis exists because there are transverse stiffeners
                        self._mesh_dim_x[i_tran - 1] = dim_x

                        if dim_x > min_x:                           # If dimension y exceeds the minimum required
                            div_round_up = np.ceil(dim_x / min_x)
                            dim_x = dim_x / div_round_up
                            self._mesh_dim_x[i_tran - 1] = dim_x

                if restriction_x is False:          # If there are no transverse stiffeners between tran1 and tran2
                    dim_x = np.amin(dim_x_list)
                    self._mesh_dim_x[i_tran - 1] = dim_x

    def get_mesh_dim_x(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Characteristic quad element x dimension for any plating zone. Returns the value based on longitudinal segment ID
         from the list of all x dimensions: mesh_dim_x.
        """
        if any(self._mesh_dim_x):
            segment_id = plate.long_seg1.id
            dim_x_id = segment_id - 1    # Trik: Segment ID za svaki nosač počinje sa 1, pa služi za identifikaciju položaja dim_x u mesh_dim_x
            dim_x = self._mesh_dim_x[dim_x_id]
            return dim_x
        else:
            print("ERROR: Mesh x dimensions list is blank! Calculate mesh quad element size first.")

    def get_mesh_dim_y(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Characteristic quad element y dimension for any plating zone. Returns the value based on transverse segment ID
         from the list of all y dimensions: mesh_dim_y.
        """
        if any(self._mesh_dim_y):
            segment_id = plate.trans_seg1.id
            dim_y_id = segment_id - 1
            dim_y = self._mesh_dim_y[dim_y_id]
            return dim_y
        else:
            print("ERROR: Mesh y dimensions list is blank! Calculate mesh quad element size first.")

    # Metoda za određivanje broja i pozicioniranje elemenata veličine dim_x, dim_y na polju oplate *** WIP
    # def get_number_of_quads(self, grillgeo, plate: Plate):
    #     """
    #     :param grillgeo: Selected grillgeo variant.
    #     :param plate: Selected plating zone.
    #     :return: Number of quad elements for any plating zone between Primary Supporting Members.
    #     """
    #
    #     bf_LS1 = plate.long_seg1.beam_prop.bf - grillgeo.corrosion_addition()[1].tc     # Net flange width of longitudinal segment 1 - L
    #     bf_LS2 = plate.long_seg2.beam_prop.bf - grillgeo.corrosion_addition()[1].tc     # Net flange width of longitudinal segment 2 - T
    #     bf_TS1 = plate.trans_seg1.beam_prop.bf - grillgeo.corrosion_addition()[1].tc    # Net flange width of transverse segment 1 - T
    #     bf_TS2 = plate.trans_seg2.beam_prop.bf - grillgeo.corrosion_addition()[1].tc    # Net flange width of transverse segment 2 - T
    #
    #     # IZRAZI VRIJEDE KADA JE JEDAN UZDUŽNI NOSAČ RUBI L PROFIL, A OSTALI T!
    #     # POSTAVITI ISPITIVANJE KOJI JE KAKAV TIP U KONAČNOM KODU
    #     rem_dist_x = plate.plate_longitudinal_dim() * 1000 - (bf_TS1 / 2) - (bf_TS2 / 2)  # Remaining distance along x axis with flange width excluded
    #     rem_dist_y = plate.plate_transverse_dim() * 1000 - bf_LS1 - (bf_LS2 / 2)
    #
    #     dim_x = self.get_mesh_dim_x(plate)
    #     dim_y = self.get_mesh_dim_y(plate)
    #
    #     # Kada bi se ostatak jednoliko rasporedio - ali ne može zbog pozicija ukrepa!
    #     """
    #     rem_dim_x = plate.plate_longitudinal_dim() * 1000 - bf_TS1 / 2 - bf_TS2 / 2     # Remaining distance along x axis with flange width excluded
    #     rem_dim_y = plate.plate_transverse_dim() * 1000 - bf_LS1 - bf_LS2 / 2           # Remaining distance along y axis with flange width excluded
    #
    #     floor_x = np.floor_divide(rem_dim_x, dim_x)     # Number of elements with dimension dim_x that fit in length rem_dim_x
    #     floor_y = np.floor_divide(rem_dim_y, dim_y)     # Number of elements with dimension dim_y that fit in length rem_dim_y
    #     mod_x = np.mod(rem_dim_x, dim_x)    # Remainder distance along x axis, [mm] - half of this distance is added to the first and last element
    #     mod_y = np.mod(rem_dim_y, dim_y)    # Remainder distance along y axis, [mm] - half of this distance is added to the first and last element
    #     print("mod_x", mod_x, ", mod_y", mod_y)
    #     print("floor_x =", floor_x, ", floor_y =", floor_y)
    #     """
    #
    #     stiffener_offset = np.round(plate.get_equal_stiffener_offset() * 1000, 4)
    #
    #     if plate.stiff_dir == BeamOrientation.TRANSVERSE:
    #         rem_dim_x1 = np.round(stiffener_offset - bf_TS1 / 2, 4)     # Element x dimension next to transverse segment 1  (left)
    #         rem_dim_x2 = np.round(stiffener_offset - bf_TS2 / 2, 4)     # Element x dimension next to transverse segment 2  (right)
    #
    #         n_elem_plate_y = plate.plate_transverse_dim() * 1000 / dim_y    # Number of elements with dim_y that fit in transverse plate dimension
    #         n_elem_dim_y = n_elem_plate_y - 2                               # Number of plate elements that will have dimension dim_y
    #
    #         l_zone_dim_y = n_elem_dim_y * dim_y
    #         rem_dim_y = (rem_dist_y - l_zone_dim_y) / 2
    #
    #         print(rem_dim_y)
    #         # print(rem_dim_x1, rem_dim_x2)
    #         # print(rem_dim_y1, rem_dim_y2)
    #
    #     # print(bf_LS1, bf_LS2, bf_TS1, bf_TS2)
    #
    #     # print("dim_x =", dim_x, ", dim_y =", dim_y)
    #     # print("rem_dist_x =", rem_dist_x, ", rem_dist_y =", rem_dist_y)
    #     # print("Stiffener offset od jakih nosača", stiffener_offset)


# Generacija čvorova i elemenata - metoda bez preklapanja *** WIP (razrađeno samo za prve dvije zone oplate)
"""
    def GenerateNodes(self, grillgeo):
        plate_id = 1
        plate_zone = grillgeo.plating()[plate_id]
        dim_x = self.plating_mesh_dim()[1][0]                                 # Element dimension along x axis, [mm]
        dim_y = self.plating_mesh_dim()[1][1]                                 # Element dimension along y axis, [mm]
        n_x = int((Plate.plate_longitudinal_dim(plate_zone) * 1000 / dim_x) + 1)    # Number of nodes along x axis
        n_y = int((Plate.plate_transverse_dim(plate_zone) * 1000 / dim_y) + 1)      # Number of nodes along y axis

        # ******** PLATING ZONE NODE ARRAYS ********
        #   2D array of all nodes on the first plating zone - reference for matching nodes and elements

        # PRVA ZONA OPLATE (POČETNA)
        plate_zone1 = np.zeros((n_y, n_x))
        node_id = 1
        for j in range(0, n_y):
            for i in range(0, n_x):
                plate_zone1[j, i] = node_id
                node_id += 1
        self._plating_nodes.append(plate_zone1)  # All plate zone arrays are saved in plate_nodes list

        # DRUGA ZONA OPLATE - KORISTI ZAJEDNIČKI POČETNI ČVOR (redak i = 0) - NEKAKO AUTOMATIZIRATI ZA OSTALE ZONE - ZAJEDNICKI SEGMENT?
        plate_zone2 = np.zeros((n_y, n_x))
        for j in range(0, n_y):
            for i in range(0, n_x):
                plate_zone2[j, 0] = plate_zone1[j, n_x - 1]
                plate_zone2[j, i] = node_id
                if i != 0:
                    node_id += 1
        self._plating_nodes.append(plate_zone2)  # All plate zone arrays are saved in plate_nodes list

        # POSTOJE 4 RAZLICITA SLUCAJA:
        #   1.) Nema zajednickih cvorova (segmenata) prilikom generacije - npr. polje 1
        #   2.) Postoje zajednicki cvorovi duz poprecnog segmenta - npr. polje 2, 3
        #   3.) Postoje zajednicki cvorovi duz uzduznog segmenta - npr. polje 4, 7, 10
        #   4.) Postoje zajednicki cvorovi duz uzduznog i poprecnog segmenta - npr. polje 5, 6, 8, 9...

        # Generate plating zone nodes
        node_id = 1
        for plate_zone in range(1, len(self._plating_nodes) + 1):
            segment1 = grillgeo.plating()[plate_zone].long_seg1     # Reference segment
            segment1_node1 = Segment.get_segment_node1(segment1)    # Reference node for generating plating zone FENodes
            dim_x = self.plating_mesh_dim()[plate_zone][0]    # Element dimension along x axis, [mm]
            dim_y = self.plating_mesh_dim()[plate_zone][1]    # Element dimension along y axis, [mm]

            for y in range(0, n_y):
                for x in range(0, n_x):
                    x_coord = segment1_node1[0] * 1000 + dim_x * x
                    y_coord = segment1_node1[1] * 1000 + dim_y * y
                    z_coord = segment1_node1[2] * 1000
                    node = FENode(node_id, x_coord, y_coord, z_coord)
                    self.add_nodes(node)
                    node_id += 1

    def GenerateElements(self, grillgeo):
        # First plating zone ID = 1:
        dim_x = self.plating_mesh_dim()[1][0]               # Element dimension along x axis, [mm]
        dim_y = self.plating_mesh_dim()[1][1]               # Element dimension along y axis, [mm]

        # Generate elements starting with element ID = 1
        element_id = 1
        for node_array in range(0, len(self._plating_nodes)):
            plate_zone = grillgeo.plating()[node_array + 1]
            n_x = int((Plate.plate_longitudinal_dim(plate_zone) * 1000 / dim_x) + 1)    # Number of nodes along x axis
            n_y = int((Plate.plate_transverse_dim(plate_zone) * 1000 / dim_y) + 1)      # Number of nodes along y axis

            plate_property_id = plate_zone.plate_prop.id

            for j in range(0, n_y - 1):
                for i in range(0, n_x - 1):
                    node1_id = self._plating_nodes[node_array][j, i]
                    node2_id = self._plating_nodes[node_array][j, i + 1]
                    node3_id = self._plating_nodes[node_array][j + 1, i + 1]
                    node4_id = self._plating_nodes[node_array][j + 1, i]

                    node1 = self.nodes()[node1_id]
                    node2 = self.nodes()[node2_id]
                    node3 = self.nodes()[node3_id]
                    node4 = self.nodes()[node4_id]

                    element = QuadElement(element_id, plate_property_id, node1, node2, node3, node4)
                    self.add_elements(element)
                    element_id += 1


class MeshData:
    def __init__(self, filename: str):
        self._filename = filename

    # Save FEM model to a file
    def write_file(self, mesh: FEMesh):
        with open(self._filename, "w") as f:
            f.write("# Node record" + "\n")
            for i in mesh.nodes():
                node = mesh.nodes()[i]
                coords = node.coords
                x = coords[0]
                y = coords[1]
                z = coords[2]

                line = str("node")
                line += ' ' + str(node.id)
                line += ' ' + str("coords")
                line += ' ' + str("6")
                line += ' ' + str(x)
                line += ' ' + str(y)
                line += ' ' + str(z) + "\n"
                f.write(line)

            f.write("# Element record" + "\n")
            for i in mesh.elements():
                element = mesh.elements()[i]
                line = str("element")
                line += ' ' + str(element.id)
                line += ' ' + str("nodes")
                line += ' ' + str("4")
                line += ' ' + str(element.node1.id)
                line += ' ' + str(element.node2.id)
                line += ' ' + str(element.node3.id)
                line += ' ' + str(element.node4.id) + "\n"
                f.write(line)
"""
