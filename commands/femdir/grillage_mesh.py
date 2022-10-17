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


class MeshSolution(Enum):
    # Meshing solution variants:
    V1 = 1  # Variant with primary supporting member flange mesh being reflected onto plating mesh
    V2 = 2  # Variant with ??


class UniquePlateProperty:
    def __init__(self, id_, tp, mat):
        self._id_ = id_
        self._tp = tp
        self._mat = mat
        self._plate_prop = []       # List of PlateProperty objects used in the grillage model with the same properties
        self._beam_prop = []        # List of BeamProperty objects used in the grillage model with the same properties

    @property
    def id(self):
        return self._id_

    @property
    def tp(self):
        return self._tp

    @tp.setter
    def tp(self, value):
        self._tp = value

    @property
    def mat(self):
        return self._mat

    @mat.setter
    def mat(self, value):
        self._mat = value

    @property
    def plate_prop(self):
        return self._plate_prop

    @plate_prop.setter
    def plate_prop(self, value):
        self._plate_prop = value

    @property
    def beam_prop(self):
        return self._beam_prop

    @beam_prop.setter
    def beam_prop(self, value):
        self._beam_prop = value


# class GrillageMesh(gfem.GeoFEM):
# class GrillageMesh(GeoFEM):
class GrillageMesh:
    def __init__(self, grillage: Grillage, variant: MeshSolution):
        # super().__init__()
        self._grillage = grillage
        self._variant = variant
        self._mesh_dim_x = []               # List of base mesh x dimensions (dim_x) in the longitudinal direction
        self._mesh_dim_y = []               # List of base mesh y dimensions (dim_y) in the transverse direction
        self._tr_el_dim_x = []              # List of transition element x dimensions
        self._tr_el_dim_y = []              # List of transition element y dimensions
        self._min_num_ebs = 1               # Minimum number of elements between stiffeners according to Rules; default = 1
        self._flange_aspect_ratio = 8       # Maximum aspect ratio value for primary supporting member flange quad elements; default = 8
        self._plate_aspect_ratio = 3        # Maximum aspect ratio value for plating and primary supporting member web quad elements; default = 3
        self._unique_properties = {}        # Dictionary of unique plate thickness and material combinations used in the grillage model

    def add_element(self, el_id, el_property, nodes):
        pass

    @property
    def mesh_dim_x(self):
        return self._mesh_dim_x

    @property
    def mesh_dim_y(self):
        return self._mesh_dim_y

    @property
    def min_num_ebs(self):
        return self._min_num_ebs

    @min_num_ebs.setter
    def min_num_ebs(self, value):
        self._min_num_ebs = value

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

    @property
    def unique_properties(self):
        return self._unique_properties

    def identify_unique_property(self, grillage):
        # Identify unique plate thickness and material properties used in plating and beams of the grillage model
        upp_id = 1
        duplicate_tp = False  # Plating thickness and material duplicate
        duplicate_tw = False  # Primary supporting member web thickness and material duplicate
        duplicate_tf = False  # Primary supporting member flange thickness and material duplicate

        # Identify unique plate properties used for plating in the grillage model
        for plate in grillage.plate_props().values():
            tp = plate.tp_net(grillage.corrosion_addition()[1], plate.tp)  # Current plate net thickness
            mat = plate.plate_mat  # Current plate material

            if self._unique_properties:  # If the dictionary is not empty
                for item in self._unique_properties.values():  # Check if values in dict have the same properties
                    if item.tp == tp and item.mat == mat:  # Duplicates have the same plate thickness and material
                        duplicate_tp = True
                        item.plate_prop.append(plate)      # Save duplicate plate object to the list
                        break
                    else:
                        duplicate_tp = False

                if duplicate_tp is False:  # Create unique property if current plate is not a duplicate
                    curr_upp = UniquePlateProperty(upp_id, tp, mat)
                    curr_upp.plate_prop.append(plate)
                    self._unique_properties[upp_id] = curr_upp
                    upp_id += 1

            else:  # Create a unique property if the dictionary is empty
                curr_upp = UniquePlateProperty(upp_id, tp, mat)
                curr_upp.plate_prop.append(plate)
                self._unique_properties[upp_id] = curr_upp
                upp_id += 1

        # Identify unique web and flange properties used for beams in the grillage model
        for beam in grillage.beam_props().values():
            if beam.beam_type == "T" or beam.beam_type == "L":
                tw = beam.tw_net(grillage.corrosion_addition()[1])  # Net web thickness
                tf = beam.tf_net(grillage.corrosion_addition()[1])  # Net flange thickness
                mat = beam.mat

                for item in self._unique_properties.values():  # Check if values in dict have the same properties as psm web
                    if item.tp == tw and item.mat == mat:
                        duplicate_tw = True
                        item.beam_prop.append(beam)  # Save duplicate plate object to the list
                        break
                    else:
                        duplicate_tw = False

                for item in self._unique_properties.values():  # Check if values in dict have the same properties as psm flange
                    if item.tp == tf and item.mat == mat:
                        duplicate_tf = True
                        item.beam_prop.append(beam)  # Save duplicate plate object to the list
                        break
                    else:
                        duplicate_tf = False

                if duplicate_tw is False:  # Create unique plate property if web plate is not a duplicate
                    curr_upp = UniquePlateProperty(upp_id, tw, mat)
                    curr_upp.beam_prop.append(beam)
                    self._unique_properties[upp_id] = curr_upp
                    upp_id += 1

                if duplicate_tf is False:  # Create unique plate property if flange plate is not a duplicate
                    curr_upp = UniquePlateProperty(upp_id, tf, mat)
                    curr_upp.beam_prop.append(beam)
                    self._unique_properties[upp_id] = curr_upp
                    upp_id += 1

            # if beam.beam_type == "FB":
            #     tw = beam.tw
            #     pass

    @staticmethod
    def find_closest_divisor(length, spacing):
        # Method for determining the number of Quad elements along length L, when one dimension (value) of the element is known
        """
        :param length: Length L which should be divided into n equal parts, each with length x.
        :param spacing: Value to which length x should be closes to.
        :return: Closest divisor of length L, which results in a length x closest to given value.
        Equivalent to the number of Finite Elements n, with dimension x.
        """
        if np.mod(length, spacing) == 0:
            return length / spacing
        else:
            i = 1
            res = []
            while i <= length:
                if np.mod(length, i) == 0:
                    res.append(i)
                i += 1

            min_diff = spacing
            min_div_id = 1
            if not res:     # If input dimensions are decimal, a divisor may not exist
                return np.round(length / spacing, 0)
            else:           # If input dimensions are integers and a divisor exists
                for i in range(0, len(res)):
                    if min_diff > abs((length / res[i]) - spacing):
                        min_diff = abs((length / res[i]) - spacing)
                        min_div_id = i
                n = res[min_div_id]
                return n

    @staticmethod
    def get_flange_el_width(segment: Segment):
        """
        :param segment: Segment of a primary supporting member.
        :return: Flange quad element dimension across the width of the flange of the selected segment.
        """
        grillage = segment.primary_supp_mem.grillage
        bf = segment.beam_prop.bf - grillage.corrosion_addition()[1].tc  # Net flange width
        if segment.beam_prop.beam_type == "T":
            dim = bf / 2                            # T profile flange is represented with 2 elements across the flange width
            return dim
        elif segment.beam_prop.beam_type == "L":
            dim = bf                                # L profile flange is represented with 1 element across the flange width
            return dim
        elif segment.beam_prop.beam_type == "FB":
            dim = 0
            return dim

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
        # koristiti np.isclose()
        pass

    def element_size_stiffener_spacing(self, plate: Plate):
        # Returns the quad element size, based only on stiffener spacing for any plating zone
        # Method is valid for determining the coarsest plating mesh allowed by the Rules for any meshing variant

        stiff_spacing = Plate.get_stiffener_spacing(plate) * 1000           # Stiffener spacing in [mm]
        n_elem_bs = self._min_num_ebs                                       # Number of elements between stiffeners
        if plate.stiff_dir == BeamDirection.LONGITUDINAL:
            L = Plate.plate_longitudinal_dim(plate) * 1000                  # Longitudinal plating zone dimension [mm]
            dim_y = stiff_spacing / n_elem_bs                               # Distance between nodes in the transverse (y) direction
            if self._variant == MeshSolution.V1:                            # Use reduced plating zone dimension for MeshSolution V1
                dim_xf1 = self.get_flange_el_width(plate.trans_seg1)        # Longitudinal segment 1 flange element x dimension
                dim_xf2 = self.get_flange_el_width(plate.trans_seg2)        # Longitudinal segment 2 flange element x dimension
                L_red = L - dim_xf1 - dim_xf2                               # Reduced longitudinal plating zone dimension
                dim_x = L_red / self.find_closest_divisor(L_red, dim_y)     # Quad element x dimension
            else:                                                           # Use full plating zone dimension for other MeshSolution
                dim_x = L / self.find_closest_divisor(L, dim_y)
            return np.array([dim_x, dim_y])

        elif plate.stiff_dir == BeamDirection.TRANSVERSE:
            B = Plate.plate_transverse_dim(plate) * 1000                    # Transverse plating zone dimension [mm]
            dim_x = stiff_spacing / n_elem_bs                               # Distance between nodes in the longitudinal (x) direction
            if self._variant == MeshSolution.V1:                            # Use reduced plating zone dimension for MeshSolution V1
                dim_yf1 = self.get_flange_el_width(plate.long_seg1)         # Transverse segment 1 flange element y dimension
                dim_yf2 = self.get_flange_el_width(plate.long_seg2)         # Transverse segment 2 flange element y dimension
                B_red = B - dim_yf1 - dim_yf2                               # Reduced transverse plating zone dimension
                dim_y = B_red / self.find_closest_divisor(B_red, dim_x)     # Quad element y dimension
            else:                                                           # Use full plating zone dimension for other MeshSolution
                dim_y = B / self.find_closest_divisor(B, dim_x)
            return np.array([dim_x, dim_y])

    def element_size_flange_width(self, segment: Segment):
        # Returns the maximum quad element length based on net flange width and maximum allowed aspect ratio for any segment
        # Method is valid for determining the coarsest flange element dimension allowed by the Rules for any meshing variant
        if self.get_flange_el_width(segment) != 0:
            max_aspect_ratio = self._flange_aspect_ratio
            dim_max = max_aspect_ratio * self.get_flange_el_width(segment)
            return dim_max
        else:
            return 0

    def element_size_plating_zone(self, plate: Plate):
        # Local consideration of base mesh dimensions dim_x and dim_y, for each plating zone individually
        # Method is limited to meshing solution V1 with psm flange mesh being reflected onto plating mesh

        plate_aspect_ratio = self._plate_aspect_ratio
        dim_p = self.element_size_stiffener_spacing(plate)          # Plate element dimensions based on stififener spacing
        dim_x = dim_p[0]                                            # Initial element x dimension is based on stiffener spacing
        dim_y = dim_p[1]                                            # Initial element y dimension is based on stiffener spacing

        dim_xf1 = self.element_size_flange_width(plate.long_seg1)   # Longitudinal segment 1 flange element x dimension
        dim_xf2 = self.element_size_flange_width(plate.long_seg2)   # Longitudinal segment 2 flange element x dimension
        dim_yf1 = self.element_size_flange_width(plate.trans_seg1)  # Transverse segment 1 flange element y dimension
        dim_yf2 = self.element_size_flange_width(plate.trans_seg2)  # Transverse segment 2 flange element y dimension

        dim_xf = np.minimum(dim_xf1, dim_xf2)           # Minimum element x dimension based on flange element aspect ratio
        dim_yf = np.minimum(dim_yf1, dim_yf2)           # Minimum element y dimension based on flange element aspect ratio

        if dim_x > dim_xf:                              # If element size based on stiffener spacing along x asis exceeds maximum flange dimension
            div_round_up = np.ceil(dim_x / dim_xf)      # Equal division of elements between transverse stiffeners
            dim_x = dim_x / div_round_up                # Base mesh dimension x refinement

            if dim_y / dim_x > plate_aspect_ratio:      # Check plating element aspect ratio after refining dim_x
                div_round_up = np.ceil(dim_y / dim_x)
                dim_y = dim_y / div_round_up            # Base mesh dimension y refinement

        if dim_y > dim_yf:                              # If element size based on stiffener spacing along y asis exceeds maximum flange dimension
            div_round_up = np.ceil(dim_y / dim_yf)      # Equal division of elements between longitudinal stiffeners
            dim_y = dim_y / div_round_up                # Base mesh dimension y refinement

            if dim_x / dim_y > plate_aspect_ratio:      # Check plating element aspect ratio after refining dim_y
                div_round_up = np.ceil(dim_x / dim_y)
                dim_x = dim_x / div_round_up            # Base mesh dimension x refinement

        return np.array((dim_x, dim_y))

    def element_size_mesh(self, grillage):
        # Global consideration of base mesh dimensions dim_x and dim_y
        # Method is limited to meshing solution V1 with psm flange mesh being reflected onto plating mesh

        if grillage.hc_variant_check() is True:         # Calculate element size only if grillage model passes hc_variant_check
            n_long = int(grillage.N_transverse - 1)     # Number of plating zones along the longitudinal axis
            n_tran = int(grillage.N_longitudinal - 1)   # Number of plating zones along the transverse axis

            self._mesh_dim_x = np.zeros(n_long)
            self._mesh_dim_y = np.zeros(n_tran)

            plating_mesh_dim_x = {}        # Dictionary of calculated dimension x for all plating zones
            plating_mesh_dim_y = {}        # Dictionary of calculated dimension y for all plating zones

            # Calculate the quad element size based on stiffener spacing and maximum allowed aspect ratio for all plating zones.
            #   Base mesh dimensions dim_x and dim_y are saved into dictionaries plating_mesh_dim_x, plating_mesh_dim_y.
            for plate in grillage.plating().values():
                plate_zones_dim = self.element_size_plating_zone(plate)
                dim_x = plate_zones_dim[0]
                dim_y = plate_zones_dim[1]

                plating_mesh_dim_x[plate.id] = dim_x
                plating_mesh_dim_y[plate.id] = dim_y

            # Assign dimension y for all plating zones between longitudinal primary supporting members
            for i_long in range(1, len(grillage.longitudinal_members())):
                long1 = grillage.longitudinal_members()[i_long]
                long2 = grillage.longitudinal_members()[i_long + 1]
                plating_zones = grillage.plating_zones_between_psm(long1, long2)    # List of all plating zones between PSM

                restriction_y = False   # Dimension dim_y is restricted by existance of longitudinal stiffeners between longitudinal psm
                dim_y_list = []         # List of element y dimensions for all plates between PSM

                min_y = 0.0
                for i in range(0, len(plating_zones)):
                    # Minimum y dimension based on flange width
                    dim_yf1 = self.element_size_flange_width(plating_zones[i].trans_seg1)
                    dim_yf2 = self.element_size_flange_width(plating_zones[i].trans_seg2)
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

                restriction_x = False   # Dimension dim_x is restricted by existance of transverse stiffeners between transverse psm
                dim_x_list = []         # List of element x dimensions for all plates between PSM

                min_x = 0.0
                for i in range(0, len(plating_zones)):
                    # Minimum x dimension based on flange width
                    dim_xf1 = self.element_size_flange_width(plating_zones[i].long_seg1)
                    dim_xf2 = self.element_size_flange_width(plating_zones[i].long_seg2)
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

    def get_base_dim_x(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Characteristic quad element x dimension for any plating zone. Returns the value based on longitudinal segment ID
                from the list of all x dimensions: mesh_dim_x. Meshing solution V1.
        """
        if any(self._mesh_dim_x):
            segment_id = plate.long_seg1.id
            dim_x_id = segment_id - 1    # Trik: Segment ID za svaki nosač počinje sa 1, pa služi za identifikaciju položaja dim_x u mesh_dim_x
            dim_x = self._mesh_dim_x[dim_x_id]
            return dim_x
        else:
            print("ERROR: Mesh x dimensions list is blank! Calculate mesh quad element size first.")

    def get_base_dim_y(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Characteristic quad element y dimension for any plating zone. Returns the value based on transverse segment ID
                from the list of all y dimensions: mesh_dim_y. Meshing solution V1.
        """
        if any(self._mesh_dim_y):
            segment_id = plate.trans_seg1.id
            dim_y_id = segment_id - 1
            dim_y = self._mesh_dim_y[dim_y_id]
            return dim_y
        else:
            print("ERROR: Mesh y dimensions list is blank! Calculate mesh quad element size first.")

    def element_size_transition_x(self, plate: Plate):
        # Local consideration of transition element mesh dimension tr_dim_x
        """
        :param plate: Selected plating zone.
        :return: Transition element dimensions x, located next to psm flange elements. First value (tr_elem_dim_x1)
                represents the element closest to the first transverse segment (trans_seg_1) and the second
                (tr_elem_dim_x2) represents the element closest to the second transverse segment (trans_seg_2)
                that define the plating zone. Meshing solution V1.
        """

        dim_x = self.get_base_dim_x(plate)
        fl_dim_x1 = self.get_flange_el_width(plate.trans_seg1)      # Flange element dimension x across the flange width of transverse segment 1
        fl_dim_x2 = self.get_flange_el_width(plate.trans_seg2)      # Flange element dimension x across the flange width of transverse segment 2

        stiff_offset = plate.get_equal_stiffener_offset() * 1000    # Stiffener offset of stiffners on the plating zone

        tr_elem_dim_x1 = 0.0
        tr_elem_dim_x2 = 0.0

        # If stiffener direction is transverse, transition element dimension x is defined by the
        #   remaining distance between stiffener offset and flange element dimension x
        if plate.stiff_dir == BeamDirection.TRANSVERSE:
            remaining_dist1 = stiff_offset - fl_dim_x1     # Remaining distance between stiffener offset and flange x dimension of transverse segment 1
            remaining_dist2 = stiff_offset - fl_dim_x2     # Remaining distance between stiffener offset and flange x dimension of transverse segment 2

            n_elem1 = np.floor(remaining_dist1 / dim_x)    # Number of elements with dimension dim_x that fit inside the remaining distance 1
            n_elem2 = np.floor(remaining_dist2 / dim_x)    # Number of elements with dimension dim_x that fit inside the remaining distance 2

            tr_elem_dim_x1 = stiff_offset - n_elem1 * dim_x - fl_dim_x1     # Transition element dimension x next to transverse segment 1
            tr_elem_dim_x2 = stiff_offset - n_elem2 * dim_x - fl_dim_x2     # Transition element dimension x next to transverse segment 2

        # If stiffener direction is longitudinal, transition element x dimension does not exist
        elif plate.stiff_dir == BeamDirection.LONGITUDINAL:
            tr_elem_dim_x1 = 0
            tr_elem_dim_x2 = 0

        # Aspect ratio check
        if tr_elem_dim_x1 != 0:
            ar1 = dim_x / tr_elem_dim_x1
            if ar1 > self._plate_aspect_ratio:
                tr_elem_dim_x1 += dim_x

        if tr_elem_dim_x2 != 0:
            ar2 = dim_x / tr_elem_dim_x2
            if ar2 > self._plate_aspect_ratio:
                tr_elem_dim_x2 += dim_x

        return np.array((tr_elem_dim_x1, tr_elem_dim_x2))

    def element_size_transition_y(self, plate: Plate):
        # Local consideration of transition element mesh dimension tr_dim_y
        """
        :param plate: Selected plating zone.
        :return: Transition element dimensions y, located next to psm flange elements. First value (tr_elem_dim_y1)
                 represents the element closest to the first longitudinal segment (long_seg1) and the second
                (tr_elem_dim_y2) represents the element closest to the second longitudinal segment (long_seg2)
                that define the plating zone. Meshing solution V1.
        """
        dim_y = self.get_base_dim_y(plate)
        fl_dim_y1 = self.get_flange_el_width(plate.long_seg1)      # Flange element dimension y across the flange width of longitudinal segment 1
        fl_dim_y2 = self.get_flange_el_width(plate.long_seg2)      # Flange element dimension y across the flange width of longitudinal segment 2

        stiff_offset = plate.get_equal_stiffener_offset() * 1000    # Stiffener offset of stiffners on the plating zone

        tr_elem_dim_y1 = 0.0
        tr_elem_dim_y2 = 0.0

        if plate.stiff_dir == BeamDirection.LONGITUDINAL:
            remaining_dist1 = stiff_offset - fl_dim_y1     # Between stiffener offset and flange y dimension of longitudinal segment 1
            remaining_dist2 = stiff_offset - fl_dim_y2     # Between stiffener offset and flange y dimension of longitudinal segment 2

            n_elem1 = np.floor(remaining_dist1 / dim_y)    # Number of elements with dimension dim_y that fit inside the remaining distance 1
            n_elem2 = np.floor(remaining_dist2 / dim_y)    # Number of elements with dimension dim_y that fit inside the remaining distance 2

            tr_elem_dim_y1 = stiff_offset - n_elem1 * dim_y - fl_dim_y1     # Transition element dimension x next to longitudinal segment 1
            tr_elem_dim_y2 = stiff_offset - n_elem2 * dim_y - fl_dim_y2     # Transition element dimension x next to longitudinal segment 2

        # If stiffener direction is transverse, transition element y dimension does not exist
        elif plate.stiff_dir == BeamDirection.TRANSVERSE:
            tr_elem_dim_y1 = 0
            tr_elem_dim_y2 = 0

        # Aspect ratio check
        if tr_elem_dim_y1 != 0:
            ar1 = dim_y / tr_elem_dim_y1
            if ar1 > self._plate_aspect_ratio:
                tr_elem_dim_y1 += dim_y

        if tr_elem_dim_y2 != 0:
            ar2 = dim_y / tr_elem_dim_y2
            if ar2 > self._plate_aspect_ratio:
                tr_elem_dim_y2 += dim_y

        return np.array((tr_elem_dim_y1, tr_elem_dim_y2))

    def transition_element_mesh(self, grillage):
        # Global consideration of transition element mesh dimensions tr_dim_x and tr_dim_y
        # Method is limited to meshing solution V1 with psm flange mesh being reflected onto plating mesh

        transition_element_dim_x = {}               # Dictionary of calculated dimension x for all plating zones
        transition_element_dim_y = {}               # Dictionary of calculated dimension y for all plating zones
        n_long = int(grillage.N_transverse - 1)     # Number of plating zones along the longitudinal axis
        n_tran = int(grillage.N_longitudinal - 1)   # Number of plating zones along the transverse axis

        self._tr_el_dim_x = np.zeros((2, n_long))
        self._tr_el_dim_y = np.zeros((2, n_tran))

        # Calculate the transition quad element size for all plating zones
        for plate in grillage.plating().values():
            tr_elem_dim_x = self.element_size_transition_x(plate)
            tr_elem_dim_y = self.element_size_transition_y(plate)
            transition_element_dim_x[plate.id] = tr_elem_dim_x
            transition_element_dim_y[plate.id] = tr_elem_dim_y

        # Assign transition element y dimensions between pairs of longitudinal members
        for i_long in range(1, len(grillage.longitudinal_members())):
            long1 = grillage.longitudinal_members()[i_long]
            long2 = grillage.longitudinal_members()[i_long + 1]
            plating_zones = grillage.plating_zones_between_psm(long1, long2)
            for i in range(0, len(plating_zones)):
                tr_elem_dim_y = transition_element_dim_y[plating_zones[i].id]
                stiff_dir = plating_zones[i].stiff_dir

                if stiff_dir == BeamDirection.LONGITUDINAL:
                    self._tr_el_dim_y[0, i_long - 1] = tr_elem_dim_y[0]
                    self._tr_el_dim_y[1, i_long - 1] = tr_elem_dim_y[1]
                    break
                else:
                    self._tr_el_dim_y[0, i_long - 1] = tr_elem_dim_y[0]
                    self._tr_el_dim_y[1, i_long - 1] = tr_elem_dim_y[1]

        # Assign transition element x dimensions between pairs of transverse members
        for i_tran in range(1, len(grillage.transverse_members())):
            tran1 = grillage.transverse_members()[i_tran]
            tran2 = grillage.transverse_members()[i_tran + 1]
            plating_zones = grillage.plating_zones_between_psm(tran1, tran2)

            for i in range(0, len(plating_zones)):
                tr_elem_dim_x = transition_element_dim_x[plating_zones[i].id]
                stiff_dir = plating_zones[i].stiff_dir

                if stiff_dir == BeamDirection.TRANSVERSE:
                    self._tr_el_dim_x[0, i_tran - 1] = tr_elem_dim_x[0]
                    self._tr_el_dim_x[1, i_tran - 1] = tr_elem_dim_x[1]
                    break
                else:
                    self._tr_el_dim_x[0, i_tran - 1] = tr_elem_dim_x[0]
                    self._tr_el_dim_x[1, i_tran - 1] = tr_elem_dim_x[1]

    def get_tr_dim_x(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Transition quad element x dimensions for any plating zone. Returns the value based on longitudinal segment ID.
                First value (dim_1) represents the element closest to the first transverse segment (trans_seg_1) and the second
                (dim_2) represents the element closest to the second transverse segment (trans_seg_2) that define the plating zone.
        """
        segment_id = plate.long_seg1.id
        dim_id = segment_id - 1
        dim = self._tr_el_dim_x
        dim_1 = dim[0][dim_id]
        dim_2 = dim[1][dim_id]
        return dim_1, dim_2

    def get_tr_dim_y(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Transition quad element x dimensions for any plating zone. Returns the value based on transverse segment ID.
                First value (dim_1) represents the element closest to the first longitudinal segment (long_seg_1) and the second
                (dim_2) represents the element closest to the second longitudinal segment (long_seg_2) that define the plating zone.
        """
        segment_id = plate.trans_seg1.id
        dim_id = segment_id - 1
        dim = self._tr_el_dim_y
        dim_1 = dim[0][dim_id]
        dim_2 = dim[1][dim_id]
        return dim_1, dim_2

    def get_tr_element_num(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Number of transition elements along x and y dimension of the plating zone.
                Method is limited to meshing solution V1 with psm flange mesh being reflected onto plating mesh.
        """

        n_tr_el_x = 2   # Number of transitional elements on the plating zone along x axis
        n_tr_el_y = 2   # Number of transitional elements on the plating zone along y axis
        if self.get_tr_dim_x(plate)[0] == 0:
            n_tr_el_x = 0
        if self.get_tr_dim_y(plate)[0] == 0:
            n_tr_el_y = 0
        return n_tr_el_x, n_tr_el_y

    def get_flange_element_num(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Number of flange elements along x and y dimension reflected onto the plating zone.
                Method is limited to meshing solution V1 with psm flange mesh being reflected onto plating mesh.
        """

        fl_dim_x1 = self.get_flange_el_width(plate.trans_seg1)
        fl_dim_x2 = self.get_flange_el_width(plate.trans_seg2)
        fl_dim_y1 = self.get_flange_el_width(plate.long_seg1)
        fl_dim_y2 = self.get_flange_el_width(plate.long_seg2)

        n_fl_el_x = 2
        n_fl_el_y = 2

        if fl_dim_x1 == 0:
            n_fl_el_x -= 1
        if fl_dim_x2 == 0:
            n_fl_el_x -= 1

        if fl_dim_y1 == 0:
            n_fl_el_y -= 1
        if fl_dim_y2 == 0:
            n_fl_el_y -= 1

        return n_fl_el_x, n_fl_el_y

    def get_element_number(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Number of all elements along x and y dimension of the plating zone.
                Method is limited to meshing solution V1 with psm flange mesh being reflected onto plating mesh.
        """

        L = plate.plate_longitudinal_dim() * 1000  # Longitudinal plating zone dimension [mm]
        B = plate.plate_transverse_dim() * 1000  # Transverse plating zone dimension [mm]

        dim_x = self.get_base_dim_x(plate)
        dim_y = self.get_base_dim_y(plate)

        tr_el_dim_x1 = self.get_tr_dim_x(plate)[0]
        tr_el_dim_x2 = self.get_tr_dim_x(plate)[1]
        tr_el_dim_y1 = self.get_tr_dim_y(plate)[0]
        tr_el_dim_y2 = self.get_tr_dim_y(plate)[1]

        fl_dim_x1 = self.get_flange_el_width(plate.trans_seg1)
        fl_dim_x2 = self.get_flange_el_width(plate.trans_seg2)
        fl_dim_y1 = self.get_flange_el_width(plate.long_seg1)
        fl_dim_y2 = self.get_flange_el_width(plate.long_seg2)

        n_tr_el_x = self.get_tr_element_num(plate)[0]
        n_tr_el_y = self.get_tr_element_num(plate)[1]

        n_fl_el_x = self.get_flange_element_num(plate)[0]
        n_fl_el_y = self.get_flange_element_num(plate)[1]

        n_el_dim_x = np.round((L - tr_el_dim_x1 - tr_el_dim_x2 - fl_dim_x1 - fl_dim_x2) / dim_x)    # Number of elements with dim_x along x axis
        n_el_dim_y = np.round((B - tr_el_dim_y1 - tr_el_dim_y2 - fl_dim_y1 - fl_dim_y2) / dim_y)    # Number of elements with dim_y along y axis

        n_el_x = n_el_dim_x + n_tr_el_x + n_fl_el_x
        n_el_y = n_el_dim_y + n_tr_el_y + n_fl_el_y
        return n_el_x, n_el_y

    def get_mesh_dim_x(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Dimensions of all quad elements along x axis, in order, for the selected plating zone.
        """
        tr_dim_x1 = self.get_tr_dim_x(plate)[0]
        tr_dim_x2 = self.get_tr_dim_x(plate)[1]
        fl_dim_x1 = self.get_flange_el_width(plate.trans_seg1)
        fl_dim_x2 = self.get_flange_el_width(plate.trans_seg2)
        base_dim_x = self.get_base_dim_x(plate)

        n_elem_tr = int(self.get_tr_element_num(plate)[0])
        n_elem_flange = int(self.get_flange_element_num(plate)[0])
        n_elem_x = int(self.get_element_number(plate)[0]) - n_elem_flange - n_elem_tr

        element_dim = {}
        element_id = 1
        start = element_id
        end = int(start + n_elem_flange / 2)    # Nedostaje provjera postoji li uopće prirubnica na ovom nosaču
        for element in range(start, end):       # pretpostavka je da postoji (zasada)
            element_dim[element_id] = fl_dim_x1
            element_id += 1

        start = end
        end = int(start + n_elem_tr / 2)
        for element in range(start, end):
            element_dim[element_id] = tr_dim_x1
            element_id += 1

        start = end
        end = int(start + n_elem_x)
        for element in range(start, end):
            element_dim[element_id] = base_dim_x
            element_id += 1

        start = end
        end = int(start + n_elem_tr / 2)
        for element in range(start, end):
            element_dim[element_id] = tr_dim_x2
            element_id += 1

        start = end
        end = int(start + n_elem_flange / 2)
        for element in range(start, end):
            element_dim[element_id] = fl_dim_x2
            element_id += 1

        return element_dim

    def get_mesh_dim_y(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Dimensions of all quad elements along y axis, in order, for the selected plating zone.
        """
        tr_dim_y1 = self.get_tr_dim_y(plate)[0]
        tr_dim_y2 = self.get_tr_dim_y(plate)[1]
        fl_dim_y1 = self.get_flange_el_width(plate.long_seg1)
        fl_dim_y2 = self.get_flange_el_width(plate.long_seg2)
        base_dim_y = self.get_base_dim_y(plate)

        n_elem_tr = int(self.get_tr_element_num(plate)[1])
        n_elem_flange = int(self.get_flange_element_num(plate)[1])
        n_elem_x = int(self.get_element_number(plate)[1]) - n_elem_flange - n_elem_tr

        element_dim = {}
        element_id = 1
        start = element_id
        end = int(start + n_elem_flange / 2)    # Nedostaje provjera postoji li uopće prirubnica na ovom nosaču
        for element in range(start, end):       # pretpostavka je da postoji (zasada)
            element_dim[element_id] = fl_dim_y1
            element_id += 1

        start = end
        end = int(start + n_elem_tr / 2)
        for element in range(start, end):
            element_dim[element_id] = tr_dim_y1
            element_id += 1

        start = end
        end = int(start + n_elem_x)
        for element in range(start, end):
            element_dim[element_id] = base_dim_y
            element_id += 1

        start = end
        end = int(start + n_elem_tr / 2)
        for element in range(start, end):
            element_dim[element_id] = tr_dim_y2
            element_id += 1

        start = end
        end = int(start + n_elem_flange / 2)
        for element in range(start, end):
            element_dim[element_id] = fl_dim_y2
            element_id += 1

        return element_dim

    def get_plate_node_coords(self, plate: Plate):
        # Node coordinates calculation for any plating zone
        # Method is limited to meshing solution V1 with psm flange mesh being reflected onto plating mesh

        # Druga verzija - vektorski sa numpy poljima
        mesh_dim_x = self.get_mesh_dim_x(plate)
        mesh_dim_y = self.get_mesh_dim_y(plate)

        ref_node1 = Segment.get_segment_node1(plate.long_seg1)              # Reference node 2
        ref_node2 = Segment.get_segment_node2(plate.long_seg1)              # Reference node 2
        ref_vector = np.subtract(ref_node2, ref_node1)                      # Reference vector in the direction of the reference segment
        unit_ref_vector = ref_vector / np.linalg.norm(ref_vector)           # Unit reference vector
        normal_vector = np.array((0, 0, 1))                                 # Vector normal to the plating surface
        perpendicular_vector = np.cross(normal_vector, unit_ref_vector)     # Unit perpendicular vector

        node_id = 1                     # Local node ID on the plating zone
        dim_y_index = 1
        spacing_vector_x = np.zeros(3)
        spacing_vector_y = np.zeros(3)

        for y_node in range(0, len(mesh_dim_y) + 1):        # Nodes along y axis
            dim_x_index = 1
            if y_node > 0:
                spacing_vector_y += mesh_dim_y[dim_y_index] * perpendicular_vector
                dim_y_index += 1
            else:
                spacing_vector_y = np.zeros(3)

            for x_node in range(0, len(mesh_dim_x) + 1):    # Nodes along x axis
                if x_node > 0:
                    spacing_vector_x += mesh_dim_x[dim_x_index] * unit_ref_vector
                    dim_x_index += 1
                else:
                    spacing_vector_x = np.zeros(3)

                spacing_vector = spacing_vector_x + spacing_vector_y    # Node position vector in the local coordinate system
                node = spacing_vector + ref_node1 * 1000   # Node coordinats in the global coordinate system, ref_node1 = position vector
                print("Node ID:", node_id, node)
                node_id += 1

        # Prva verzija koda, direktan izračun koordinata x,y,z:    Radi OK
        """
        mesh_dim_x = self.get_mesh_dim_x(plate)
        mesh_dim_y = self.get_mesh_dim_y(plate)

        segment = plate.long_seg1
        grillage = segment.primary_supp_mem.grillage
        ref_node = grillage.get_segment_nodes(segment)[0]

        node_id = 1
        x_coord = 0
        y_coord = 0

        j = 1
        for y in range(0, len(mesh_dim_y) + 1):
            i = 1
            if y == 0:
                y_coord = ref_node[1] * 1000
            else:
                y_coord += mesh_dim_y[j]
                j += 1

            for x in range(0, len(mesh_dim_x) + 1):
                if x == 0:
                    x_coord = ref_node[0] * 1000
                else:
                    x_coord += mesh_dim_x[i]
                    i += 1
                print("Node ID:", node_id, ", koordinate:", x_coord, y_coord)
                node_id += 1
        """

    def get_all_node_coords(self, grillage):
        # Get all plating node coordinates

        # Ovo nije dobra implementacija za generaciju čvorova polovičnog modela!
        #   Ovakva jednostavna petlja kroz pola zona ne funkcionira za paran broj nosača, gdje os simetrije siječe centralnu zonu oplate na pola
        # Za polovični ili četvrtinski model je potrebna provjera ima li poklopac paran ili neparan broj nosača oko osi simetrije
        # Što se događa s elementima na osi simetrije?
        # Kako napraviti novi beam property za ukrepu na osi simetrije ili novi jaki nosač sa pola širine struka na osi simetrije?
        #   - riješiti na razini mreže, dodjelom drugačijeg beam property? Postoji li već implementirano rješenje za ovo?

        # Preko prve verzije:
        """
        node_id = 1
        n_plate = int(len(grillage.plating().keys()) / 2)
        print(n_plate)

        for plate_id in range(1, n_plate + 1):
            plate = grillage.plating()[plate_id]
            mesh_dim_x = self.get_mesh_dim_x(plate)
            mesh_dim_y = self.get_mesh_dim_y(plate)

            segment = plate.long_seg1
            grillage = segment.primary_supp_mem.grillage
            ref_node = grillage.get_segment_nodes(segment)[0]

            x_coord = 0
            y_coord = 0
            z_coord = segment.beam_prop.hw

            j = 1

            for y in range(0, len(mesh_dim_y) + 1):
                i = 1
                if y == 0:
                    y_coord = ref_node[1] * 1000
                else:
                    y_coord += mesh_dim_y[j]
                    j += 1

                for x in range(0, len(mesh_dim_x) + 1):
                    if x == 0:
                        x_coord = ref_node[0] * 1000
                    else:
                        x_coord += mesh_dim_x[i]
                        i += 1
                    coords = np.array((x_coord, y_coord, z_coord))
                    print("Plate ID:", plate.id, "Node ID:", node_id, ", koordinate:", coords)
                    node_id += 1

        """

        # Preko vektora:
        node_id = 1
        n_plate = int(len(grillage.plating().keys()) / 2)
        print(n_plate)

        for plate_id in range(1, n_plate + 1):
            plate = grillage.plating()[plate_id]
            mesh_dim_x = self.get_mesh_dim_x(plate)
            mesh_dim_y = self.get_mesh_dim_y(plate)

            ref_node1 = Segment.get_segment_node1(plate.long_seg1)  # Reference node 2
            ref_node2 = Segment.get_segment_node2(plate.long_seg1)  # Reference node 2
            ref_vector = np.subtract(ref_node2, ref_node1)  # Reference vector in the direction of the reference segment
            unit_ref_vector = ref_vector / np.linalg.norm(ref_vector)  # Unit reference vector
            normal_vector = np.array((0, 0, 1))  # Vector normal to the plating surface
            perpendicular_vector = np.cross(normal_vector, unit_ref_vector)  # Unit perpendicular vector

            dim_y_index = 1
            spacing_vector_x = np.zeros(3)
            spacing_vector_y = np.zeros(3)

            for y_node in range(0, len(mesh_dim_y) + 1):  # Nodes along y axis
                dim_x_index = 1
                if y_node > 0:
                    spacing_vector_y += mesh_dim_y[dim_y_index] * perpendicular_vector
                    dim_y_index += 1
                else:
                    spacing_vector_y = np.zeros(3)

                for x_node in range(0, len(mesh_dim_x) + 1):  # Nodes along x axis
                    if x_node > 0:
                        spacing_vector_x += mesh_dim_x[dim_x_index] * unit_ref_vector
                        dim_x_index += 1
                    else:
                        spacing_vector_x = np.zeros(3)

                    spacing_vector = spacing_vector_x + spacing_vector_y  # Node position vector in the local coordinate system
                    node = spacing_vector + ref_node1 * 1000  # Node coordinats in the global coordinate system, ref_node1 = position vector
                    print("Node ID:", node_id, node)
                    node_id += 1


# Generacija čvorova i elemenata - metoda bez preklapanja *** WIP (razrađeno samo za prve dvije zone oplate)
"""
    def GenerateNodes(self, grillage):
        plate_id = 1
        plate_zone = grillage.plating()[plate_id]
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

        # POSTOJE 4 RAZLICITA SLUCAJA: - puno if provjera
        #   1.) Nema zajednickih cvorova (segmenata) prilikom generacije - npr. polje 1
        #   2.) Postoje zajednicki cvorovi duz poprecnog segmenta - npr. polje 2, 3
        #   3.) Postoje zajednicki cvorovi duz uzduznog segmenta - npr. polje 4, 7, 10
        #   4.) Postoje zajednicki cvorovi duz uzduznog i poprecnog segmenta - npr. polje 5, 6, 8, 9...

        # Generate plating zone nodes
        node_id = 1
        for plate_zone in range(1, len(self._plating_nodes) + 1):
            segment1 = grillage.plating()[plate_zone].long_seg1     # Reference segment
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

    def GenerateElements(self, grillage):
        # First plating zone ID = 1:
        dim_x = self.plating_mesh_dim()[1][0]               # Element dimension along x axis, [mm]
        dim_y = self.plating_mesh_dim()[1][1]               # Element dimension along y axis, [mm]

        # Generate elements starting with element ID = 1
        element_id = 1
        for node_array in range(0, len(self._plating_nodes)):
            plate_zone = grillage.plating()[node_array + 1]
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
"""
