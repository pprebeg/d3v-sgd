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
# from femdir.geofementity import *


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
        self._min_num_eweb = 3              # Minimum number of elements representing the web of a psm along its height; default = 3
        self._num_eaf = 1                   # Number of elements across primary supporting member flange; default = 1
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
    def min_num_eweb(self):
        return self._min_num_eweb

    @min_num_eweb.setter
    def min_num_eweb(self, value):
        self._min_num_eweb = value

    @property
    def num_eaf(self):
        return self._num_eaf

    @num_eaf.setter
    def num_eaf(self, value):
        self._num_eaf = value

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

    # VIŠE SE NE KORISTI - Zapravo nije bilo potrebe za cjelobrojnom podjelom dimenzija preko ove metode,
    #  dozvoljene su decimalne vrijednosti dimenzija elemenata! Metoda je korisna kao zaseban kalkulator

    # @staticmethod
    # def find_closest_divisor(length, spacing):
    #     # Method for determining the number of Quad elements along length L, when one dimension (value) of the element is known
    #     """
    #     :param length: Length L which should be divided into n equal parts, each with length x.
    #     :param spacing: Value to which length x should be closes to.
    #     :return: Closest divisor of length L, which results in a length x closest to given value.
    #     Equivalent to the number of Finite Elements n, with dimension x.
    #     """
    #     if np.mod(length, spacing) == 0:
    #         return length / spacing
    #     else:
    #         i = 1
    #         res = []
    #         while i <= length:
    #             if np.mod(length, i) == 0:
    #                 res.append(i)
    #             i += 1
    #
    #         min_diff = spacing
    #         min_div_id = 1
    #         if not res:     # If input dimensions are decimal, a divisor may not exist
    #             return np.round(length / spacing, 0)
    #         else:           # If input dimensions are integers and a divisor exists
    #             for i in range(0, len(res)):
    #                 if min_diff > abs((length / res[i]) - spacing):
    #                     min_diff = abs((length / res[i]) - spacing)
    #                     min_div_id = i
    #             n = res[min_div_id]
    #             return n

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

    @staticmethod
    def element_aspect_ratio(dim_x, dim_y):
        """
        :param dim_x: Dimension of the quad element along x axis.
        :param dim_y: Dimension of the quad element along y axis.
        :return: Aspect ratio of the quad element.
        """
        if dim_x > dim_y:
            ar = dim_x / dim_y
        elif dim_x < dim_y:
            ar = dim_y / dim_x
        else:
            ar = 1
        return ar

    def refine_plate_element(self, length, dim_limit, perpendicular_dim):
        """
        :param length: Dimension which is to be equally divided with elements of size dim.
        :param dim_limit: Maximum dimension allowed for the element along the given length.
        :param perpendicular_dim: Element dimension perpendicular to the given length which determines the aspect ratio.
        :return: Element dimension that is less than the maximum allowed and results in plate element dimension within
                allowed plating aspect ratio.
        """
        n_elements = np.ceil(length / dim_limit)    # Minimum number of elements that equally divide the given length, with dim < dim_limit
        dim = length / n_elements
        ar = self.element_aspect_ratio(dim, perpendicular_dim)

        if perpendicular_dim / self._plate_aspect_ratio >= dim_limit:   # Impossible condition, skip while loop and return dim
            return dim

        else:           # Increase the number of elements until aspect ratio falls under the maximum allowed for plate elements
            while ar > self.plate_aspect_ratio:
                n_elements += 1
                dim = length / n_elements
                ar = self.element_aspect_ratio(dim, perpendicular_dim)
            return dim

    def element_size_perp_to_stiffeners(self, plate: Plate):
        # Returns quad element size perpendicular to stiffener direction
        # Method is valid for determining the maximum element dimension allowed by the Rules for any meshing variant
        stiff_spacing = Plate.get_stiffener_spacing(plate) * 1000   # Stiffener spacing in [mm]
        n_elem_bs = self._min_num_ebs                               # Number of elements between stiffeners
        perpendicular_dim = stiff_spacing / n_elem_bs                             # Quad element dimension between stiffeners
        return perpendicular_dim

    def get_reduced_plate_dim(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Reduced plate dimensions based on plate stiffener orientation for MeshSolution.V1.
        """
        if plate.stiff_dir == BeamDirection.LONGITUDINAL:
            L = plate.plate_longitudinal_dim() * 1000          # Longitudinal plating zone dimension [mm]
            dim_xf1 = self.get_flange_el_width(plate.trans_seg1)    # Transverse segment 1 flange element x dimension
            dim_xf2 = self.get_flange_el_width(plate.trans_seg2)    # Transverse segment 2 flange element x dimension
            L_red = L - dim_xf1 - dim_xf2                           # Reduced longitudinal plating zone dimension for MeshSolution.V1
            return L_red

        elif plate.stiff_dir == BeamDirection.TRANSVERSE:
            B = plate.plate_transverse_dim() * 1000            # Transverse plating zone dimension [mm]
            dim_yf1 = self.get_flange_el_width(plate.long_seg1)     # Longitudinal segment 1 flange element y dimension
            dim_yf2 = self.get_flange_el_width(plate.long_seg2)     # Longitudinal segment 2 flange element y dimension
            B_red = B - dim_yf1 - dim_yf2                           # Reduced transverse plating zone dimension for MeshSolution.V1
            return B_red

    def get_flange_el_width(self, segment: Segment):
        """
        :param segment: Segment of a primary supporting member.
        :return: Flange quad element dimension across the width of the flange of the selected segment (perpendicular to the
                segment direction). For longitudinal segments this method returns dimension dim_yf, and for transverse segments
                returns dimension dim_xf. Method is valid for determining the coarsest flange element dimension allowed by the
                Rules for any meshing variant
        """
        grillage = segment.primary_supp_mem.grillage
        bf = segment.beam_prop.bf - grillage.corrosion_addition()[1].tc  # Net flange width
        if segment.beam_prop.beam_type == "T":
            dim = bf / (self._num_eaf * 2)      # T profile flange is represented with 2 elements across the flange width by default
            return dim
        elif segment.beam_prop.beam_type == "L":
            dim = bf / self._num_eaf            # L profile flange is represented with 1 element across the flange width by default
            return dim
        elif segment.beam_prop.beam_type == "FB":
            dim = 0
            return dim

    def get_flange_el_length(self, segment: Segment):
        """
        :param segment: Segment of a primary supporting member.
        :return: Maximum flange quad element length based on flange element width and maximum flange aspect ratio (parallel to the
                segment direction). For longitudinal segments this method returns dimension dim_xf, and for transverse segments
                returns dimension dim_yf. Method is valid for determining the coarsest flange element dimension allowed by the
                Rules for any meshing variant
        """
        if self.get_flange_el_width(segment) != 0:
            max_aspect_ratio = self._flange_aspect_ratio
            dim_max = max_aspect_ratio * self.get_flange_el_width(segment)
            return dim_max
        else:
            return 0

    def get_web_el_height(self, segment: Segment):
        """
        :param segment: Segment of a primary supporting member.
        :return: Quad element dimension along the height of a primary supporting member.
        """
        hw = segment.beam_prop.hw
        dim = hw / self._min_num_eweb
        return dim

    def element_size_stiffener_spacing(self, plate: Plate):
        # Returns the quad element size, based only on stiffener spacing with aspect ratio close to 1
        # Method is valid for determining the coarsest plating mesh allowed by the Rules for any meshing variant
        if plate.stiff_dir == BeamDirection.LONGITUDINAL:
            L = plate.plate_longitudinal_dim() * 1000          # Longitudinal plating zone dimension [mm]
            dim_y = self.element_size_perp_to_stiffeners(plate)     # Distance between nodes in the transverse (y) direction

            if self._variant == MeshSolution.V1:                    # Use reduced plating zone dimension for MeshSolution V1
                L = self.get_reduced_plate_dim(plate)               # Reduced longitudinal plating zone dimension

            div_round = np.round(L / dim_y)
            dim_x = L / div_round                                   # Quad element x dimension
            return np.array([dim_x, dim_y])

        elif plate.stiff_dir == BeamDirection.TRANSVERSE:
            B = plate.plate_transverse_dim() * 1000            # Transverse plating zone dimension [mm]
            dim_x = self.element_size_perp_to_stiffeners(plate)     # Distance between nodes in the longitudinal (x) direction

            if self._variant == MeshSolution.V1:                    # Use reduced plating zone dimension for MeshSolution V1
                B = self.get_reduced_plate_dim(plate)               # Reduced transverse plating zone dimension

            div_round = np.round(B / dim_x)
            dim_y = B / div_round                                   # Quad element y dimension
            return np.array([dim_x, dim_y])

    def element_size_plating_zone(self, plate: Plate):
        # Local consideration of base mesh dimensions dim_x and dim_y, for each plating zone individually
        # Method is limited to meshing solution V1 with psm flange mesh being reflected onto plating mesh

        plate_aspect_ratio = self._plate_aspect_ratio
        dim_p = self.element_size_stiffener_spacing(plate)      # Plate element dimensions based on stififener spacing
        dim_x = dim_p[0]                                        # Initial element x dimension is based on stiffener spacing
        dim_y = dim_p[1]                                        # Initial element y dimension is based on stiffener spacing

        dim_xf1 = self.get_flange_el_length(plate.long_seg1)    # Longitudinal segment 1 flange element x dimension
        dim_xf2 = self.get_flange_el_length(plate.long_seg2)    # Longitudinal segment 2 flange element x dimension
        dim_yf1 = self.get_flange_el_length(plate.trans_seg1)   # Transverse segment 1 flange element y dimension
        dim_yf2 = self.get_flange_el_length(plate.trans_seg2)   # Transverse segment 2 flange element y dimension

        dim_xf = np.minimum(dim_xf1, dim_xf2)                   # Maximum element x dimension based on flange element aspect ratio
        dim_yf = np.minimum(dim_yf1, dim_yf2)                   # Maximum element y dimension based on flange element aspect ratio

        if dim_x > dim_xf:
            if plate.stiff_dir == BeamDirection.LONGITUDINAL:
                dim_x = self.refine_plate_element(self.get_reduced_plate_dim(plate), dim_xf, dim_y)     # Refine along L_red, max dimension dim_xf, aspect ratio dim_y
                if self.element_aspect_ratio(dim_x, dim_y) > plate_aspect_ratio:                        # Check plating element aspect ratio after refining dim_x
                    dim_y = self.refine_plate_element(plate.get_stiffener_spacing() * 1000, dim_yf, dim_x)

            elif plate.stiff_dir == BeamDirection.TRANSVERSE:
                dim_x = self.refine_plate_element(plate.get_stiffener_spacing() * 1000, dim_xf, dim_y)  # Refine along stiff spacing, max dim_xf, aspect ratio dim_y
                if self.element_aspect_ratio(dim_x, dim_y) > plate_aspect_ratio:                        # Check plating element aspect ratio after refining dim_x
                    dim_y = self.refine_plate_element(self.get_reduced_plate_dim(plate), dim_yf, dim_x)

        if dim_y > dim_yf:
            if plate.stiff_dir == BeamDirection.LONGITUDINAL:
                dim_y = self.refine_plate_element(plate.get_stiffener_spacing() * 1000, dim_yf, dim_x)
                if self.element_aspect_ratio(dim_x, dim_y) > plate_aspect_ratio:
                    dim_x = self.refine_plate_element(self.get_reduced_plate_dim(plate), dim_xf, dim_y)

            elif plate.stiff_dir == BeamDirection.TRANSVERSE:
                dim_y = self.refine_plate_element(self.get_reduced_plate_dim(plate), dim_yf, dim_x)
                if self.element_aspect_ratio(dim_x, dim_y) > plate_aspect_ratio:
                    dim_x = self.refine_plate_element(plate.get_stiffener_spacing() * 1000, dim_x, dim_y)

        # STARI KOD:
        """
        if dim_x > dim_xf:                              # If element size based on stiffener spacing along x asis exceeds maximum flange dimension
            # **** Moguće poboljšati kod da se dobije nešto grublja mreža provjerom orijentacije ukrepa! ****
            # Ako je orijentacija BeamDirection.LONGITUDINAL dimx se može odrediti prema broju elemenata dim_xf koji stane između el prirubnica
            # Ako je orijentacija BeamDirection.TRANSVERSE dimx mora biti višekratnik razmaka između ukrepa - postojeći kod

            div_round_up = np.ceil(dim_x / dim_xf)      # Equal division of elements between transverse stiffeners
            dim_x = dim_x / div_round_up                # Base mesh dimension x refinement

            ar = self.element_aspect_ratio(dim_x, dim_y)
            if ar > plate_aspect_ratio:      # Check plating element aspect ratio after refining dim_x
                div_round_up = np.ceil(dim_y / dim_x)
                dim_y = dim_y / div_round_up            # Base mesh dimension y refinement

        if dim_y > dim_yf:                              # If element size based on stiffener spacing along y asis exceeds maximum flange dimension
            div_round_up = np.ceil(dim_y / dim_yf)      # Equal division of elements between longitudinal stiffeners
            dim_y = dim_y / div_round_up                # Base mesh dimension y refinement

            ar = self.element_aspect_ratio(dim_x, dim_y)
            if ar > plate_aspect_ratio:      # Check plating element aspect ratio after refining dim_y
                div_round_up = np.ceil(dim_x / dim_y)
                dim_x = dim_x / div_round_up            # Base mesh dimension x refinement
        """
        return np.array((dim_x, dim_y))

    def element_size_mesh(self, grillage):
        # Global consideration of base mesh dimensions dim_x and dim_y
        # Method is limited to meshing solution V1 with psm flange mesh being reflected onto plating mesh

        if grillage.hc_variant_check() is True:         # Calculate element size only if grillage model passes hc_variant_check
            n_long = int(grillage.N_transverse - 1)     # Number of plating zones along the longitudinal axis
            n_tran = int(grillage.N_longitudinal - 1)   # Number of plating zones along the transverse axis

            self._mesh_dim_x = np.zeros(n_long)         # Final base mesh dimension x between transverse primary supporting members
            self._mesh_dim_y = np.zeros(n_tran)         # Final base mesh dimension y between longitudinal primary supporting members

            plating_mesh_dim_x = {}                     # Dimension x for all plating zones, based on element_size_plating_zone
            plating_mesh_dim_y = {}                     # Dimension y for all plating zones, based on element_size_plating_zone

            # Calculate the quad element size based on stiffener spacing and maximum allowed aspect ratio for all plating zones
            for plate in grillage.plating().values():
                plate_zones_dim = self.element_size_plating_zone(plate)
                plating_mesh_dim_x[plate.id] = plate_zones_dim[0]
                plating_mesh_dim_y[plate.id] = plate_zones_dim[1]

            # Assign dimension y for all plating zones between longitudinal primary supporting members
            for i_long in range(1, len(grillage.longitudinal_members())):
                long1 = grillage.longitudinal_members()[i_long]
                long2 = grillage.longitudinal_members()[i_long + 1]
                plating_zones = grillage.plating_zones_between_psm(long1, long2)    # List of all plating zones between PSM

                if any(plate.stiff_dir == BeamDirection.LONGITUDINAL for plate in plating_zones):
                    for plate in plating_zones:
                        if plate.stiff_dir == BeamDirection.LONGITUDINAL:   # If dimension y is limited by longitudinal stiffener spacing

                            # Find maximum allowed value of base mesh dim_y
                            max_y = np.Infinity
                            segments_list = grillage.segments_between_psm(long1, long2)  # All transverse segments between long1 and long2
                            for segment in segments_list:
                                dim_yf = self.get_flange_el_length(segment)  # Maximum quad element y dimension based on net flange width
                                if dim_yf < max_y:
                                    max_y = dim_yf

                            dim_y = plating_mesh_dim_y[plate.id]  # Base element size in the y direction for plate in list plating_zones
                            if dim_y > max_y:                               # If dimension y exceeds the maximum allowed
                                dim_y = self.refine_plate_element(plate.get_stiffener_spacing() * 1000, max_y, plating_mesh_dim_x[plate.id])
                                self._mesh_dim_y[i_long - 1] = dim_y        # Save value of dim_y
                            else:                                           # If dim_y does not exceed the maximum max_y
                                self._mesh_dim_y[i_long - 1] = dim_y
                            break                                           # Stop after finding the first zone with longitudinal stiffeners
                else:
                    # Duga verzija:
                    """
                    dim_y_list = []  # List of element y dimensions for all plates between adjacent longitudinal PSM
                    for plate in plating_zones:  # ** Postoji prostor za optimizaciju! -eliminacija provjere duplikatne vrijednosti segmenata
                        dim_y = plating_mesh_dim_y[plate.id]  # Quad element size in the y direction for plate in list plating_zones
                        dim_y_list.append(dim_y)
                    """
                    dim_y_list = [plating_mesh_dim_y[plate.id] for plate in plating_zones]
                    dim_y = np.amin(dim_y_list)  # Use minimum value of all saved dim_y for plating zones between longitudinal psm
                    self._mesh_dim_y[i_long - 1] = dim_y

            # Assign dimension x for all plating zones between transverse primary supporting members
            for i_tran in range(1, len(grillage.transverse_members())):
                tran1 = grillage.transverse_members()[i_tran]
                tran2 = grillage.transverse_members()[i_tran + 1]
                plating_zones = grillage.plating_zones_between_psm(tran1, tran2)  # List of all plating zones between PSM

                if any(plate.stiff_dir == BeamDirection.TRANSVERSE for plate in plating_zones):
                    for plate in plating_zones:
                        if plate.stiff_dir == BeamDirection.TRANSVERSE:     # If dimension x is limited by transverse stiffener spacing

                            # Find maximum allowed value of base mesh dim_x
                            max_x = np.Infinity
                            segments_list = grillage.segments_between_psm(tran1, tran2)  # All longitudinal segments between tran1 and tran2
                            for segment in segments_list:
                                dim_xf = self.get_flange_el_length(segment)  # Maximum quad element x dimension based on net flange width
                                if dim_xf < max_x:
                                    max_x = dim_xf

                            dim_x = plating_mesh_dim_x[plate.id]            # Quad element size in the x direction for plate in list plating_zones
                            if dim_x > max_x:                               # If dimension x exceeds the maximum allowed
                                dim_x = self.refine_plate_element(plate.get_stiffener_spacing() * 1000, max_x, plating_mesh_dim_y[plate.id])
                                self._mesh_dim_x[i_tran - 1] = dim_x        # Save value of dim_x
                            else:
                                self._mesh_dim_x[i_tran - 1] = dim_x
                            break                                           # Stop after finding the first zone with transverse stiffeners
                else:
                    dim_x_list = [plating_mesh_dim_x[plate.id] for plate in plating_zones]
                    dim_x = np.amin(dim_x_list)  # Use minimum value of all saved dim_x for plating zones between transverse psm
                    self._mesh_dim_x[i_tran - 1] = dim_x

    def get_base_dim_x(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Base quad element x dimension for any plating zone. Returns the value based on longitudinal segment ID
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
        :return: Base quad element y dimension for any plating zone. Returns the value based on transverse segment ID
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
        dim_y = self.get_base_dim_y(plate)
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
            ar1 = dim_y / tr_elem_dim_x1
            if ar1 > self._plate_aspect_ratio:
                tr_elem_dim_x1 += dim_x

        if tr_elem_dim_x2 != 0:
            ar2 = dim_y / tr_elem_dim_x2
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
        dim_x = self.get_base_dim_x(plate)
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
            ar1 = dim_x / tr_elem_dim_y1
            if ar1 > self._plate_aspect_ratio:
                tr_elem_dim_y1 += dim_y

        if tr_elem_dim_y2 != 0:
            ar2 = dim_x / tr_elem_dim_y2
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

    def generate_plating_zone_elements(self, plate: Plate, start_node_id, start_element_id, split_along=AOS.NONE):
        """
        :param plate: Selected plate for calculating node locations, generating Noode and Element objects.
        :param start_node_id: Starting node ID which allows continued numeration after other methods.
        :param start_element_id: Starting element ID which allows continued numeration after other methods.
        :param split_along: Variable determines the meshing limits of the selected plating zone based on Axis Of Symmetry.
        :return: Determines node coordinates and generates finite element Node and Element objects on the selected plating zone.
                Returns last node and element ID, to continue node and element numbering on the next plating zone or segment.
                Method is valid for any meshing solution variant MeshSolution using quad elements,
                given different input dictionaries of dimensions of all elements along x and y axis.
        """

        print("Čvorovi zone oplate", plate.id)

        mesh_dim_x = self.get_mesh_dim_x(plate)     # Input dictionary of dimensions of all elements along x axis for MeshSolution.V1
        mesh_dim_y = self.get_mesh_dim_y(plate)     # Input dictionary of dimensions of all elements along y axis for MeshSolution.V1

        row_limit = len(mesh_dim_y) + 1     # Number of node rows on the entire plating zone
        column_limit = len(mesh_dim_x) + 1  # Number of node columns on the entire plating zone

        # **** DODATAK **** WIP
        # Ako zonu oplate presjeca os simetrije:

        # Modify row and column limits based on split_along Axis Of Symmetry input
        if split_along == AOS.LONGITUDINAL:     # Longitudinal axis of symmetry splits the plating zone
            pass
        elif split_along == AOS.TRANSVERSE:     # Transverse axis of symmetry splits the plating zone
            pass
        elif split_along == AOS.BOTH:           # Both longitudinal and transverse axis of symmetry splits the plating zone
            pass

        # Node coordinates calculation
        node_id = start_node_id
        ref_node1 = Segment.get_segment_node1(plate.long_seg1)              # Reference node 1 coordinates in [mm]
        ref_node2 = Segment.get_segment_node2(plate.long_seg1)              # Reference node 2 coordinates in [mm]
        ref_vector = np.subtract(ref_node2, ref_node1)                      # Reference vector in the direction of the reference segment
        unit_ref_vector = ref_vector / np.linalg.norm(ref_vector)           # Unit reference vector
        normal_vector = np.array((0, 0, 1))                                 # Vector normal to the plating surface
        perpendicular_vector = np.cross(normal_vector, unit_ref_vector)     # Unit perpendicular vector

        spacing_vector_x = np.zeros(3)
        spacing_vector_y = np.zeros(3)
        dim_y_index = 1
        for row in range(0, row_limit):        # Row of nodes along x axis
            dim_x_index = 1
            if row > 0:
                spacing_vector_y += mesh_dim_y[dim_y_index] * perpendicular_vector
                dim_y_index += 1
            else:
                spacing_vector_y = np.zeros(3)

            for column in range(0, column_limit):    # Column of nodes along y axis
                if column > 0:
                    spacing_vector_x += mesh_dim_x[dim_x_index] * unit_ref_vector
                    dim_x_index += 1
                else:
                    spacing_vector_x = np.zeros(3)

                spacing_vector = spacing_vector_x + spacing_vector_y    # Node position vector in the local coordinate system
                node = spacing_vector + ref_node1   # Node coordinats in the global coordinate system, ref_node1 = position vector
                # Ovdje instanciraj objekt čvora
                print("Node ID:", node_id, ", koordinate:", node)
                node_id += 1
        last_node_id = node_id      # Return last node ID for continued node numeration

        # Node ID array - reference for quad element generation
        plate_id_array = np.zeros((row_limit, column_limit))
        node_id = start_node_id
        for row in range(0, row_limit):
            for column in range(0, column_limit):
                plate_id_array[row, column] = node_id
                node_id += 1
        # print(plate_id_array)

        # Identify UniquePlateProperty
        plate_prop = plate.plate_prop
        for unique_prop in self._unique_properties.values():
            for grillage_prop in unique_prop.plate_prop:
                if grillage_prop is plate_prop:
                    # print("Zoni oplate pripada jedinstveni property", unique_prop.id, ", čiji je ekvivalent u grillage modelu", plate_prop.id)
                    pass

        print(" \n Elementi zone oplate", plate.id)

        # Generate elements
        element_id = start_element_id
        for row in range(0, row_limit - 1):
            for column in range(0, column_limit - 1):
                node1_id = plate_id_array[row, column]
                node2_id = plate_id_array[row, column + 1]
                node3_id = plate_id_array[row + 1, column + 1]
                node4_id = plate_id_array[row + 1, column]
                # Ovdje instanciraj objekt quad elementa, dodaj čvorove, svojstvo elementu, itd.
                # Quad elementu dodijeliti svojstvo identificirano u self._unique_properties
                print("Quad element ID", element_id, ", Node IDs:", node1_id, node2_id, node3_id, node4_id)
                element_id += 1

        return last_node_id, element_id

    def generate_stiffener_elements(self):
        pass
        # WIP
        # Nova ideja: preko koordinata čvora identificiraj samo prvi čvor na retku ili stupcu čvorova gdje se nalazi ukrepa, a
        #   ostali čvorovi pojedine ukrepe su poznati preko smjera generacije čvorova i orijentacije ukrepa
        #   ubaciti ovo u metodu za izradu elemenata oplate? - već postoji referentna matrica sa ID čvorova
        """
        for i_stiff in range(1, int(plate.get_stiffener_number()) + 1):
            stiffener_end_coords = plate.get_stiff_coords(i_stiff)
            ref_node1 = stiffener_end_coords[0]
            print("Početni čvor grednog elementa ukrepe:", ref_node1)

            # Identify node
            #   pretraživanje instanciranih čvorova, pronaći koji čvor ima koordinate ref_node1
        """

        """
        # Izračun koordinata svih čvorova svih ukrepa na zoni oplate 
        
        element_length = None
        if plate.stiff_dir == BeamDirection.LONGITUDINAL:
            element_length = self.get_mesh_dim_x(plate)
        elif plate.stiff_dir == BeamDirection.TRANSVERSE:
            element_length = self.get_mesh_dim_y(plate)

        spacing_vector = np.zeros(3)
        for i_stiff in range(1, int(plate.get_stiffener_number()) + 1):
            stiffener_end_coords = plate.get_stiff_coords(i_stiff)
            ref_node1 = stiffener_end_coords[0]
            ref_node2 = stiffener_end_coords[1]
            ref_vector = np.subtract(ref_node2, ref_node1)  # Reference vector in the direction of the reference segment
            unit_ref_vector = ref_vector / np.linalg.norm(ref_vector)  # Unit reference vector

            dim_index = 1
            for node in range(0, len(element_length)):
                if node > 0:
                    spacing_vector += element_length[dim_index] * unit_ref_vector
                    second_spacing_vector = element_length[dim_index + 1] * unit_ref_vector
                    dim_index += 1
                else:
                    spacing_vector = np.zeros(3)
                    second_spacing_vector = element_length[1] * unit_ref_vector

                node1_coords = spacing_vector + ref_node1
                node2_coords = second_spacing_vector + node1_coords
                print("Čvorovi grednog elementa ukrepe:", node1_coords, node2_coords)
        """

    def generate_segment_web_elements(self, segment: Segment, start_node_id, start_element_id, split=False):
        """
        :param segment: Selected segment for calculating node locations, generating Noode and Element objects.
        :param start_node_id: Starting node ID which allows continued numeration after other methods.
        :param start_element_id: Starting element ID which allows continued numeration after other methods.
        :param split: Variable determines the meshing limits of the selected segment based on Axis Of Symmetry.
        :return: Determines node coordinates and generates finite element Node and Element objects on the selected segment.
                Returns last node and element ID, to continue node and element numbering on the next segment.
                Method is limited to meshing solutions with uniform web mesh and no deformed quad elements.
        """

        print("Čvorovi segmenta", segment.id, "jakog nosača", segment.primary_supp_mem.id)

        direction = segment.primary_supp_mem.direction
        grillage = segment.primary_supp_mem.grillage
        dim_z = self.get_web_el_height(segment)                 # Vertical dimension of every segment web element

        mesh_long_dim = None
        for plate in grillage.plating().values():               # Identify which plating zone the segment belongs to
            segment_defines_plate = plate.test_plate_segment(segment)
            if segment_defines_plate:
                if direction == BeamDirection.LONGITUDINAL:
                    mesh_long_dim = self.get_mesh_dim_x(plate)  # Mesh dimensions along the entire length of the segment
                elif direction == BeamDirection.TRANSVERSE:
                    mesh_long_dim = self.get_mesh_dim_y(plate)
                break                                           # Stop after finding the first plating zone segment belongs to

        # **** DODATAK **** WIP
        # Ako segment presjeca os simetrije:

        # Modify column limits based on split input
        column_limit = len(mesh_long_dim) + 1                   # Number of node columns on the entire segment
        if split:
            pass

        # Node coordinates calculation
        node_id = start_node_id
        ref_node1 = Segment.get_segment_node1(segment)          # Reference node 1 coordinates in [mm], origin of the local csy
        ref_node2 = Segment.get_segment_node2(segment)          # Reference node 2 coordinates in [mm]
        ref_vector = np.subtract(ref_node2, ref_node1)          # Reference vector in the direction of the reference segment
        unit_ref_vector = ref_vector / np.linalg.norm(ref_vector)  # Unit reference vector
        perpendicular_vector = np.array((0, 0, -1))             # Vector in the direction of psm flange, opposite of global z axis direction
        long_spacing_vector = np.zeros(3)                       # Lonngitudinal spacing vector in the direction of segment psm

        for row in range(0, self._min_num_eweb + 1):            # Total number of rows of web element nodes is equal to min_num_ewb + 1
            vertical_spacing_vector = perpendicular_vector * dim_z * row
            long_dim_index = 1
            for column in range(0, column_limit):
                if column > 0:
                    long_spacing_vector += mesh_long_dim[long_dim_index] * unit_ref_vector
                    long_dim_index += 1
                else:
                    long_spacing_vector = np.zeros(3)

                position_vector = long_spacing_vector + vertical_spacing_vector    # Node position vector in the local coordinate system
                node_coords = position_vector + ref_node1
                print("Node ID:", node_id, node_coords)
                node_id += 1

        # Node ID array - reference for quad element generation
        plate_id_array = np.zeros((self._min_num_eweb + 1, column_limit))
        node_id = start_node_id
        for row in range(0, self._min_num_eweb + 1):
            for column in range(0, column_limit):
                plate_id_array[row, column] = node_id
                node_id += 1
        # print(plate_id_array)

        print(" \n Elementi segmenta", segment.id, "jakog nosača", segment.primary_supp_mem.id)

        # Generate elements
        element_id = start_element_id
        for row in range(0, self._min_num_eweb):
            for column in range(0, column_limit - 1):
                node1_id = plate_id_array[row, column]
                node2_id = plate_id_array[row, column + 1]
                node3_id = plate_id_array[row + 1, column + 1]
                node4_id = plate_id_array[row + 1, column]
                # Ovdje instanciraj objekt quad elementa, dodaj čvorove, svojstvo elementu, itd.
                # Quad elementu dodijeliti svojstvo identificirano u self._unique_properties
                print("Quad element ID", element_id, ", Node IDs:", node1_id, node2_id, node3_id, node4_id)
                element_id += 1

        return node_id, element_id

    def generate_segment_flange_elements(self, segment: Segment, start_node_id):
        # WIP
        # return node_id, element_id
        pass

    def generate_mesh(self, grillage, AxisOfSymm: AOS):
        plate_limit = len(grillage.plating())               # Generate plating mesh on all plating zones
        long_member_limit = grillage.N_longitudinal         # Generate mesh on all longitudinal members
        tran_member_limit = grillage.N_transverse           # Generate mesh on all transverse members
        long_segment_limit = grillage.N_transverse - 1      # Generate mesh on all longitudinal segments
        tran_segment_limit = grillage.N_longitudinal - 1    # Generate mesh on all transverse segments

        # **** WIP - Potrebno smisliti logiku za postavljanje limita prema odabranoj osi simetrije
        #   odgovoriti na pitanje na kojim zonama oplate i na kojim segmentima će se raditi mreža ako je odabrana neka os simetrije
        if AxisOfSymm == AOS.LONGITUDINAL:  # Grillage has longitudinal axis of symmetry
            pass
        elif AxisOfSymm == AOS.TRANSVERSE:  # Grillage has transverse axis of symmetry
            pass
        elif AxisOfSymm == AOS.BOTH:        # Grillage has both longitudinal and transverse axis of symmetry
            pass
        elif AxisOfSymm == AOS.NONE:        # Grillage no axis of symmetry
            pass

        starting_node_id = 1
        starting_element_id = 1

        # Generate all plating quad elements and stiffener beam elements
        for plate_id in range(1, plate_limit + 1):
            plate = grillage.plating()[plate_id]
            last_IDs = self.generate_plating_zone_elements(plate, starting_node_id, starting_element_id)
            starting_node_id = last_IDs[0]
            starting_element_id = last_IDs[1]
            # self.generate_stiffener_elements(plate)

        # Generate all longitudinal segment quad elements
        for long_id in range(1, long_member_limit + 1):
            longitudinal = grillage.longitudinal_members()[long_id]
            for segment_id in range(0, long_segment_limit):
                segment = longitudinal.segments[segment_id]
                last_IDs = self.generate_segment_web_elements(segment, starting_node_id, starting_element_id)
                starting_node_id = last_IDs[0]
                starting_element_id = last_IDs[1]
                # last_IDs = self.generate_segment_flange_elements(segment, starting_node_id)
                # starting_node_id = last_IDs[0]
                # starting_element_id = last_IDs[1]

        # Generate all transverse segment quad elements
        for tran_id in range(1, tran_member_limit + 1):
            transverse = grillage.transverse_members()[tran_id]
            for segment_id in range(0, tran_segment_limit):
                segment = transverse.segments[segment_id]
                last_IDs = self.generate_segment_web_elements(segment, starting_node_id, starting_element_id)
                starting_node_id = last_IDs[0]
                starting_element_id = last_IDs[1]
                # last_IDs = self.generate_segment_flange_elements(segment, starting_node_id)
                # starting_node_id = last_IDs[0]
                # starting_element_id = last_IDs[1]
