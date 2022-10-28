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


# class MeshSize(GeoFEM):
class MeshSize:
    """
    ISPRAVITI:
    Određivanje dimenzija bi trebalo biti sukladno odabranoj osi simetrije! Ako je odabrana izrada četvrtinskog modela,
    nema smisla provjeravati i usklađivati dimenzije na cijelom poklopcu. Na ovaj način nestaje problem presjecanja elementa sa osi simetrije.
    Izraditi metode za identifikaciju kakva mreža treba biti na kojoj zoni oplate.
    """
    def __init__(self, grillage: Grillage, axis_of_symm=AOS.NONE):
        # super().__init__()
        self._grillage = grillage           # Input grillage model
        self._axis_of_symm = axis_of_symm   # Input axis of symmetry of the grillage model

        self._min_num_ebs = 1               # Minimum number of elements between stiffeners; default = 1
        self._min_num_eweb = 3              # Minimum number of elements representing the web of a psm along its height; default = 3
        self._num_eaf = 1                   # Number of elements across primary supporting member flange; default = 1
        self._flange_aspect_ratio = 8       # Maximum aspect ratio value for primary supporting member flange quad elements; default = 8
        self._plate_aspect_ratio = 3        # Maximum aspect ratio value for plating and primary supporting member web quad elements; default = 3
        self._des_plate_aspect_ratio = 2    # Desirable plating aspect ratio value less than the maximum; default = 2
        self._unique_properties = {}        # Dictionary of unique plate thickness and material combinations used in the grillage model

        self.full_plate_zones = {}          # Plating zones to be fully meshed
        self.long_half_plate_zones = {}     # Plating zones to be split with a longitudinal axis of symmetry
        self.tran_half_plate_zones = {}     # Plating zones to be split with a transverse axis of symmetry
        self.quarter_plate_zone = {}        # Plating zone to be split with both axis of symmetry

    @property
    def axis_of_symm(self):
        return self._axis_of_symm

    @axis_of_symm.setter
    def axis_of_symm(self, value):
        self._axis_of_symm = value

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
    def des_plate_aspect_ratio(self):
        return self._des_plate_aspect_ratio

    @des_plate_aspect_ratio.setter
    def des_plate_aspect_ratio(self, value):
        self._des_plate_aspect_ratio = value
        if value > self.plate_aspect_ratio:
            raise Exception("Desired plate element aspect ratio can not be greater than the maximum plate aspect ratio!")

    def plate_zone_reference_ID_array(self):
        """
        :return: 2D array of plating zone IDs arranged to represent relative placement of zones on the grillage model.
        """
        total_rows = self._grillage.N_longitudinal - 1        # Total number of plating zone rows on the grillage model
        total_columns = self._grillage.N_transverse - 1       # Total number of plating zone columns on the grillage model

        plating_zone_array = np.zeros((total_rows, total_columns))
        zone_id = 1

        for row in range(0, total_rows):
            for column in range(0, total_columns):
                plating_zone_array[row, column] = zone_id
                zone_id += 1

        return plating_zone_array

    def longitudinal_psm_extent(self):
        """
        :return: Dictionary of longitudinal primary supporting members to be considered for mesh generation,
                based on input Axis of Symmetry.
        """
        longitudinals = {}

        if np.mod(self._grillage.N_longitudinal, 2) == 0:
            n_long = self._grillage.N_longitudinal / 2
        else:
            n_long = np.ceil(self._grillage.N_longitudinal / 2)

        if self.axis_of_symm == AOS.LONGITUDINAL or self.axis_of_symm == AOS.BOTH:
            for i in range(1, int(n_long) + 1):
                longitudinals[i] = self._grillage.longitudinal_members()[i]
            return longitudinals
        else:
            return self._grillage.longitudinal_members()

    def transverse_psm_extent(self):
        """
        :return: Dictionary of transverse primary supporting members to be considered for mesh generation,
                based on input Axis of Symmetry.
        """
        transversals = {}

        if np.mod(self._grillage.N_transverse, 2) == 0:
            n_tran = self._grillage.N_transverse / 2
        else:
            n_tran = np.ceil(self._grillage.N_transverse / 2)

        if self.axis_of_symm == AOS.TRANSVERSE or self.axis_of_symm == AOS.BOTH:
            for i in range(1, int(n_tran) + 1):
                transversals[i] = self._grillage.transverse_members()[i]
            return transversals
        else:
            return self._grillage.transverse_members()

    def longitudinal_segment_extent(self):
        pass

    def transverse_segment_extent(self):
        pass

    def identify_long_full_plate_zones(self):
        """
        :return: Identifies Plate objects for full mesh generation on the entire plating zone; grillage with a longitudinal axis of symm.
                Stores identified zones in full_plate_zones dictionary.
        """
        total_rows = self._grillage.N_longitudinal - 1    # Total number of plating zone rows on the grillage model
        plating_zone_array = self.plate_zone_reference_ID_array()

        n_full_rows = int(np.floor(total_rows / 2))
        zone_arr_split = plating_zone_array[0:n_full_rows, :]
        for i in zone_arr_split:
            for j in i:
                self.full_plate_zones[j] = self._grillage.plating()[j]

    def identify_tran_full_plate_zones(self):
        """
        :return: Identifies Plate objects for full mesh generation on the entire plating zone; grillage with a transverse axis of symm.
                Stores identified zones in full_plate_zones dictionary.
        """
        total_columns = self._grillage.N_transverse - 1   # Total number of plating zone columns on the grillage model
        plating_zone_array = self.plate_zone_reference_ID_array()

        n_full_columns = int(np.floor(total_columns / 2))
        zone_arr_split = plating_zone_array[:, 0:n_full_columns]
        for i in zone_arr_split:
            for j in i:
                self.full_plate_zones[j] = self._grillage.plating()[j]

    def identify_both_full_plate_zones(self):
        """
        :return: Identifies Plate objects for full mesh generation on the entire plating zone; grillage with both axis of symm.
                Stores identified zones in full_plate_zones dictionary.
        """
        total_rows = self._grillage.N_longitudinal - 1    # Total number of plating zone rows on the grillage model
        total_columns = self._grillage.N_transverse - 1   # Total number of plating zone columns on the grillage model
        plating_zone_array = self.plate_zone_reference_ID_array()

        middle_row_index = int(np.floor(total_rows / 2))
        middle_column_index = int(np.floor(total_columns / 2))
        zone_arr_split = plating_zone_array[0:middle_row_index, 0:middle_column_index]
        for i in zone_arr_split:
            for j in i:
                self.full_plate_zones[j] = self._grillage.plating()[j]

    def identify_long_half_plate_zones(self):
        """
        :return: Identifies Plate objects for half mesh generation; plating zones split with a longitudinal axis of symmetry.
                Stores identified zones in long_half_plate_zones dictionary.
        """
        total_rows = self._grillage.N_longitudinal - 1              # Total number of plating zone rows on the grillage model
        plating_zone_array = self.plate_zone_reference_ID_array()

        if np.mod(self._grillage.N_longitudinal, 2) == 0:           # Odd number of plating zone rows
            middle_row_index = int(np.floor(total_rows / 2))
            middle_row = plating_zone_array[middle_row_index, :]    # Find middle row of plating zone index array
            for i in middle_row:
                self.long_half_plate_zones[i] = self._grillage.plating()[i]

    def identify_tran_half_plate_zones(self):
        """
        :return: Identifies Plate objects for half mesh generation; plating zones split with a transverse axis of symmetry.
                Stores identified zones in tran_half_plate_zones dictionary.
        """
        total_columns = self._grillage.N_transverse - 1   # Total number of plating zone columns on the grillage model
        plating_zone_array = self.plate_zone_reference_ID_array()

        if np.mod(self._grillage.N_transverse, 2) == 0:                 # Odd number of plating zone columns
            middle_column_index = int(np.floor(total_columns / 2))
            middle_column = plating_zone_array[:, middle_column_index]  # Find middle column of plating zone index array
            for i in middle_column:
                self.tran_half_plate_zones[i] = self._grillage.plating()[i]

    def identify_both_half_plate_zones(self):
        """
        :return: Identifies Plate objects for half mesh generation; grillage with both axis of symm.
                Stores identified zones in long_half_plate_zones and tran_half_plate_zones dictionary.
        """
        total_rows = self._grillage.N_longitudinal - 1              # Total number of plating zone rows on the grillage model
        total_columns = self._grillage.N_transverse - 1   # Total number of plating zone columns on the grillage model
        plating_zone_array = self.plate_zone_reference_ID_array()
        middle_row_index = int(np.floor(total_rows / 2))
        middle_column_index = int(np.floor(total_columns / 2))

        if np.mod(self._grillage.N_longitudinal, 2) == 0:           # Odd number of plating zone rows
            middle_row = plating_zone_array[middle_row_index, 0:middle_column_index]        # Find middle row of plating zone index array
            for i in middle_row:
                self.long_half_plate_zones[i] = self._grillage.plating()[i]

        if np.mod(self._grillage.N_transverse, 2) == 0:                 # Odd number of plating zone columns
            middle_column = plating_zone_array[0:middle_row_index, middle_column_index]     # Find middle column of plating zone index array
            for i in middle_column:
                self.tran_half_plate_zones[i] = self._grillage.plating()[i]

    def identify_quarter_plate_zone(self):
        """
        :return: Identifies Plate object for quarter mesh generation; plating zone split with both axis of symmetry.
                 Only possible if grillage has both axis of symmetry, with even number of longitudinal and transverse
                 primary supporting members. There can be only one plating zone split with both axis of symmetry.
                 Identified zone is saved in quarter_plate_zone dictionary for consistency with other methods.
        """
        if np.mod(self._grillage.N_longitudinal, 2) == 0 and np.mod(self._grillage.N_transverse, 2) == 0:
            total_rows = self._grillage.N_longitudinal - 1      # Total number of plating zone rows on the grillage model
            total_columns = self._grillage.N_transverse - 1     # Total number of plating zone columns on the grillage model
            plating_zone_array = self.plate_zone_reference_ID_array()
            middle_row_index = int(np.floor(total_rows / 2))
            middle_column_index = int(np.floor(total_columns / 2))
            middle_plate_ID = plating_zone_array[middle_row_index, middle_column_index]
            self.quarter_plate_zone[1] = self._grillage.plating()[middle_plate_ID]

    def grillage_plate_extent(self):
        """
        :return: Determines limits of plate mesh generation based on input Axis of Symmetry value.
                Calls specific methods for identifying which plating zones will be fully or partially meshed.
                If grillage has no axis of symmetry, all plating zones inside grillage.plating() will be meshed.
        """
        if self.axis_of_symm == AOS.LONGITUDINAL:
            self.identify_long_full_plate_zones()
            self.identify_long_half_plate_zones()

        elif self.axis_of_symm == AOS.TRANSVERSE:
            self.identify_tran_full_plate_zones()
            self.identify_tran_half_plate_zones()

        elif self.axis_of_symm == AOS.BOTH:
            self.identify_both_full_plate_zones()
            self.identify_both_half_plate_zones()
            self.identify_quarter_plate_zone()

        else:
            self.full_plate_zones = self._grillage.plating()

    @property
    def unique_properties(self):
        return self._unique_properties

    def identify_unique_property(self):
        # Identify unique plate thickness and material properties used in plating and beams of the grillage model
        upp_id = 1
        duplicate_tp = False  # Plating thickness and material duplicate
        duplicate_tw = False  # Primary supporting member web thickness and material duplicate
        duplicate_tf = False  # Primary supporting member flange thickness and material duplicate

        # Identify unique plate properties used for plating in the grillage model
        for plate in self._grillage.plate_props().values():
            tp = plate.tp_net(self._grillage.corrosion_addition()[1], plate.tp)  # Current plate net thickness
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
        for beam in self._grillage.beam_props().values():
            if beam.beam_type == "T" or beam.beam_type == "L":
                tw = beam.tw_net(self._grillage.corrosion_addition()[1])  # Net web thickness
                tf = beam.tf_net(self._grillage.corrosion_addition()[1])  # Net flange thickness
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
        """
        :param length: Length L which should be divided into n equal parts, each with length x.
        :param spacing: Value to which length x should be closes to.
        :return: Closest divisor of length L, which results in a length x closest to given spacing value.
                If a divisor does not exist, this method returns None.
        """
        if np.mod(length, spacing) == 0:
            return length / spacing

        else:
            i = 1
            res = []
            while i <= length:
                if np.mod(length, i) == 0:
                    res.append(i)
                    # print(i)    # Provjera djelitelja
                i += 1

            min_diff = spacing
            min_div_id = 1

            if not res:     # If input length is decimal and a divisor does not exist
                return None
            else:           # If input dimensions are integers and a divisor exists
                for i in range(0, len(res)):
                    if min_diff > abs((length / res[i]) - spacing):
                        min_diff = abs((length / res[i]) - spacing)
                        min_div_id = i
                n = res[min_div_id]
                return n

    @staticmethod
    def find_largest_divisor(length, max_val):
        """
        :param length: Length L which should be divided into n equal parts, each with length x.
        :param max_val: Maximum value which length x can not exceed.
        :return: Largest divisor of length L, which results in a length x smaller or equal to the given maximum value.
                If a divisor does not exist, such as when length L is decimal, this method returns None.
        """
        i = 1
        divisors = []
        while i <= length:
            if np.mod(length, i) == 0:
                divisors.append(i)
            i += 1

        curr_x_max = 0                  # Current maximum value of x
        n = None

        if not divisors:                # If input length is decimal and a divisor does not exist
            return None
        else:                           # If a divisor exists
            for i in divisors:
                curr_x = length / i     # Current dimension x value
                if curr_x_max < curr_x <= max_val:
                    curr_x_max = curr_x
                    n = i
            return n

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

    @staticmethod
    def refine_plate_element(length, dim_limit):
        """
        :param length: Dimension which is to be equally divided with elements of size dim.
        :param dim_limit: Maximum dimension allowed for the element along the given length, defined by aspect ratio or other limits.
        :return: Element dimension value that is less than the maximum allowed and equally divides given length.
        """
        n_elements = np.ceil(length / dim_limit)  # Maximum number of elements that equally divide the given length, with dim < dim_limit
        dim = length / n_elements
        return dim

    def CheckNodeOverlap(self):
        # koristiti np.isclose()
        pass

    def element_size_perp_to_stiffeners(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Quad element dimension in [mm] perpendicular to stiffener direction based on number of elements between stiffeners.
        """
        stiff_spacing = plate.get_stiffener_spacing() * 1000    # Stiffener spacing in [mm]
        n_elem_bs = self._min_num_ebs                           # Number of elements between stiffeners
        perp_dim = stiff_spacing / n_elem_bs                    # Quad element dimension between stiffeners
        return perp_dim

    def element_size_para_to_stiffeners(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Quad element dimension parlallel to stiffener direction, based only on quad element dimension
                 perpendicular to stiffener direction, desired and maximum aspect ratio. Method attempts to find
                 a divisor of plating zone dimension parallel to stiffeners that results in quad element size not
                 greater than the desired element dimension parallel to stiffener direction. If no divisors are found
                 or the divisor is the plating zone dimension itself, method refine_plate_element is used.
        """
        L = plate.plate_dim_parallel_to_stiffeners()                # Reduced plating zone dimension parallel to stiffener direction
        y = self.element_size_perp_to_stiffeners(plate)             # Element dimension perpendicular to stiffener direction

        des_x_val = y * self.des_plate_aspect_ratio                 # Desired element dimension parallel to stiffener direction
        n_elem = self.find_largest_divisor(L, des_x_val)            # Number of elements with desired dimension des_x_val
        if n_elem is not None:
            x = L / n_elem                                          # Element dimension parallel to stiffener direction
            ar = self.element_aspect_ratio(x, y)
            if ar > self.plate_aspect_ratio:                        # Check aspect ratio in case L is a prime number
                x = self.refine_plate_element(L, des_x_val)
        else:                                                       # There is no divisor in case length L is not a integer
            x = self.refine_plate_element(L, des_x_val)
        return x

    def get_flange_el_width(self, segment: Segment):
        """
        :param segment: Segment of a primary supporting member.
        :return: Flange quad element dimension across the width of the flange of the selected segment (perpendicular to the
                segment direction) based on the number of elements across the flange. For longitudinal segments this method
                returns dimension dim_yf and for transverse segments returns dimension dim_xf.
        """
        bf_net = TBeamProperty.bf_net(segment.beam_prop, self._grillage.corrosion_addition()[1])

        if segment.beam_prop.beam_type == "T":
            return bf_net / (self._num_eaf * 2)  # T profile flange is represented with 2 elements across the flange width by default

        elif segment.beam_prop.beam_type == "L":
            return bf_net / self._num_eaf        # L profile flange is represented with 1 element across the flange width by default

        elif segment.beam_prop.beam_type == "FB":
            return 0

    def get_flange_el_length(self, segment: Segment):
        """
        :param segment: Segment of a primary supporting member.
        :return: Maximum flange quad element length based on flange element width and maximum flange aspect ratio (parallel to the
                segment direction). For longitudinal segments this method returns dimension dim_xf, and for transverse segments
                returns dimension dim_yf.
        """
        if self.get_flange_el_width(segment) != 0:
            return self._flange_aspect_ratio * self.get_flange_el_width(segment)
        else:
            return 0

    def get_min_flange_el_length(self, segment1: Segment, segment2: Segment):
        """
        :param segment1: First selected segment.
        :param segment2: Second selected segment.
        :return: Minimum value of two segment flange element lengths. If both segments do not have a flange, method returns 0.
        """
        dim_f1 = self.get_flange_el_length(segment1)
        dim_f2 = self.get_flange_el_length(segment2)

        if dim_f1 != 0 and dim_f2 != 0:
            min_dim = np.minimum(dim_f1, dim_f2)
        elif dim_f1 != 0 and dim_f2 == 0:
            min_dim = dim_f1
        elif dim_f1 == 0 and dim_f2 != 0:
            min_dim = dim_f2
        else:
            min_dim = 0

        return min_dim

    def get_min_flange_el_length_between_psm(self, member1: PrimarySuppMem, member2: PrimarySuppMem):
        """
        :param member1: First selected Primary supporting member.
        :param member2: Second selected Primary supporting member.
        :return: Method identifies all segments between two adjacent primary supporting members and returns the minimum flange
                element length of all identified segments. If all flange element length values are 0, the el_length_list will
                be empty and this method will return 0.
        """
        segments_list = self._grillage.segments_between_psm(member1, member2)  # All transverse segments between member1 and member2
        el_length_list = [self.get_flange_el_length(segment) for segment in segments_list if self.get_flange_el_length(segment) != 0]
        if not el_length_list:
            return 0
        else:
            return np.amin(el_length_list)

    def get_web_el_height(self, segment: Segment):
        """
        :param segment: Segment of a primary supporting member.
        :return: Quad element dimension along the height of a primary supporting member.
        """
        hw = segment.beam_prop.hw
        dim = hw / self._min_num_eweb
        return dim


class MeshV1(MeshSize):
    """
    *** Ovo bi trebala biti klasa za određivanje dimenzija mreže specifičnih za varijantu V1 sa preslikavanjem prirubnica ***
    Dobivene dimenzije osnovnih (base) i prijelaznih (transition) elemenata su ulazni podaci za izradu svih čvorova i elemenata.

    Razmisliti da te dimenzije ne budu strogo dimenzije elemenata, već položaji čvorova duž ruba elemenata gdje čvorovi moraju biti.
    Što se događa s mrežom unutar pojedinog elementa je do odabira varijante izrade mreže tog dijela, ali ti rubni čvorovi moraju postojati.

    Mesh Variant V1 reflects primary supporting member flange elements onto plating and has the following limitations:
        1.) Flange width of a primary supporting member has to be the same on all of its segments.
        2.) All primary supporting members need to have the same web height.
        3.) Flange element overlap has to be in the same plane.
        4.) Grillage plating can not be defined with any camber.
    """
    def __init__(self, grillage: Grillage, axis_of_symm=AOS.NONE):
        super().__init__(grillage, axis_of_symm)
        self._mesh_dim_x = []               # List of base mesh x dimensions (dim_x) in the longitudinal direction
        self._mesh_dim_y = []               # List of base mesh y dimensions (dim_y) in the transverse direction
        self._tr_el_dim_x = []              # List of transition element x dimensions
        self._tr_el_dim_y = []              # List of transition element y dimensions

    @property
    def mesh_dim_x(self):
        return self._mesh_dim_x

    @property
    def mesh_dim_y(self):
        return self._mesh_dim_y

    def get_reduced_plate_dim(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Reduced plate dimensions based on plate stiffener orientation for MeshVariantV1.
        """
        if plate.stiff_dir == BeamDirection.LONGITUDINAL:
            L = plate.plate_longitudinal_dim() * 1000               # Longitudinal plating zone dimension [mm]
            dim_xf1 = self.get_flange_el_width(plate.trans_seg1)    # Transverse segment 1 flange element x dimension
            dim_xf2 = self.get_flange_el_width(plate.trans_seg2)    # Transverse segment 2 flange element x dimension
            L_red = L - dim_xf1 - dim_xf2                           # Reduced longitudinal plating zone dimension
            return L_red

        elif plate.stiff_dir == BeamDirection.TRANSVERSE:
            B = plate.plate_transverse_dim() * 1000                 # Transverse plating zone dimension [mm]
            dim_yf1 = self.get_flange_el_width(plate.long_seg1)     # Longitudinal segment 1 flange element y dimension
            dim_yf2 = self.get_flange_el_width(plate.long_seg2)     # Longitudinal segment 2 flange element y dimension
            B_red = B - dim_yf1 - dim_yf2                           # Reduced transverse plating zone dimension
            return B_red

    def element_size_para_to_stiffeners(self, plate: Plate):
        """
        :param plate:
        :return: Quad element dimension parlallel to stiffener direction, based only on quad element dimension
                 perpendicular to stiffener direction, desired and maximum aspect ratio. Method attempts to find
                 a divisor of plating zone dimension parallel to stiffeners that results in quad element size not
                 greater than the desired element dimension parallel to stiffener direction. If no divisors are found
                 or the divisor is the plating zone dimension itself, method refine_plate_element is used.
        """
        L = self.get_reduced_plate_dim(plate)                 # Reduced plating zone dimension parallel to stiffener direction
        y = self.element_size_perp_to_stiffeners(plate)       # Element dimension perpendicular to stiffener direction
        # L = 4133
        # L = 2090
        # y = 935

        des_x_val = y * self.des_plate_aspect_ratio                 # Desired element dimension parallel to stiffener direction
        n_elem = self.find_largest_divisor(L, des_x_val)            # Number of elements with desired dimension des_x_val
        if n_elem is not None:
            x = L / n_elem                                          # Element dimension parallel to stiffener direction
            ar = self.element_aspect_ratio(x, y)
            if ar > self.plate_aspect_ratio:                        # Check aspect ratio in case L is a prime number
                x = self.refine_plate_element(L, des_x_val)
        else:                                                       # There is no divisor in case length L is not a integer
            x = self.refine_plate_element(L, des_x_val)
        return x

    def element_size_plating_zone(self, plate: Plate):
        # Local consideration of base mesh dimensions dim_x and dim_y, for each plating zone individually
        if plate.stiff_dir == BeamDirection.LONGITUDINAL:
            dim_x = self.element_size_para_to_stiffeners(plate)
            dim_y = self.element_size_perp_to_stiffeners(plate)
        else:
            dim_x = self.element_size_perp_to_stiffeners(plate)
            dim_y = self.element_size_para_to_stiffeners(plate)

        dim_xf = self.get_min_flange_el_length(plate.long_seg1, plate.long_seg2)
        dim_yf = self.get_min_flange_el_length(plate.trans_seg1, plate.trans_seg2)

        if dim_x > dim_xf != 0:
            if plate.stiff_dir == BeamDirection.LONGITUDINAL:
                dim_x = self.refine_plate_element(self.get_reduced_plate_dim(plate), dim_xf)     # Refine along L_red, max dimension dim_xf, aspect ratio dim_y
                if self.element_aspect_ratio(dim_x, dim_y) > self._plate_aspect_ratio:                        # Check plating element aspect ratio after refining dim_x
                    dim_y_limit = np.minimum(dim_yf, self.element_size_perp_to_stiffeners(plate))           # Maximum allowed dim_y dimension
                    dim_y = self.refine_plate_element(plate.get_stiffener_spacing() * 1000, dim_y_limit)

            elif plate.stiff_dir == BeamDirection.TRANSVERSE:
                dim_x_limit = np.minimum(dim_xf, self.element_size_perp_to_stiffeners(plate))               # Maximum allowed dim_x dimension
                dim_x = self.refine_plate_element(plate.get_stiffener_spacing() * 1000, dim_x_limit)  # Refine along stiff spacing, max dim_xf, aspect ratio dim_y
                if self.element_aspect_ratio(dim_x, dim_y) > self._plate_aspect_ratio:                        # Check plating element aspect ratio after refining dim_x
                    dim_y = self.refine_plate_element(self.get_reduced_plate_dim(plate), dim_yf)

        if dim_y > dim_yf != 0:
            if plate.stiff_dir == BeamDirection.LONGITUDINAL:
                dim_y_limit = np.minimum(dim_yf, self.element_size_perp_to_stiffeners(plate))  # Maximum allowed dim_y dimension
                dim_y = self.refine_plate_element(plate.get_stiffener_spacing() * 1000, dim_y_limit)
                if self.element_aspect_ratio(dim_x, dim_y) > self._plate_aspect_ratio:
                    dim_x = self.refine_plate_element(self.get_reduced_plate_dim(plate), dim_xf)

            elif plate.stiff_dir == BeamDirection.TRANSVERSE:
                dim_y = self.refine_plate_element(self.get_reduced_plate_dim(plate), dim_yf)
                if self.element_aspect_ratio(dim_x, dim_y) > self._plate_aspect_ratio:
                    dim_x_limit = np.minimum(dim_xf, self.element_size_perp_to_stiffeners(plate))  # Maximum allowed dim_x dimension
                    dim_x = self.refine_plate_element(plate.get_stiffener_spacing() * 1000, dim_x_limit)

        return dim_x, dim_y

    def element_base_size_mesh(self):
        # Global consideration of base mesh dimensions dim_x and dim_y

        # ISPRAVITI:
        # ********** Implementirati provjeru osi simetrije i sukladno odabrati osnovne dimenzije ********
        # Ne petlja kroz sve zone oplate i nosače cijelog modela nego samo dijela za koji se radi mreža!.

        if self._grillage.hc_variant_check() is True:         # Calculate element size only if grillage model passes hc_variant_check
            n_long = int(self._grillage.N_transverse - 1)     # Number of plating zones along the longitudinal axis
            n_tran = int(self._grillage.N_longitudinal - 1)   # Number of plating zones along the transverse axis

            self._mesh_dim_x = np.zeros(n_long)         # Final base mesh dimension x between transverse primary supporting members
            self._mesh_dim_y = np.zeros(n_tran)         # Final base mesh dimension y between longitudinal primary supporting members

            plating_mesh_dim_x = {}                     # Dimension x for all plating zones, based on element_size_plating_zone
            plating_mesh_dim_y = {}                     # Dimension y for all plating zones, based on element_size_plating_zone

            # Calculate the quad element size based on stiffener spacing and maximum allowed aspect ratio for all plating zones
            for plate in self._grillage.plating().values():
                dim_x, dim_y = self.element_size_plating_zone(plate)
                plating_mesh_dim_x[plate.id] = dim_x
                plating_mesh_dim_y[plate.id] = dim_y

            # ******** Možda razbiti na više metoda koje zasebno usklađuju osnovne x i y dimenzije mreže *************

            # Assign dimension y for all plating zones between longitudinal primary supporting members
            for i_long in range(1, len(self._grillage.longitudinal_members())):
                long1 = self._grillage.longitudinal_members()[i_long]
                long2 = self._grillage.longitudinal_members()[i_long + 1]
                plating_zones = self._grillage.plating_zones_between_psm(long1, long2)    # List of all plating zones between PSM

                if any(plate.stiff_dir == BeamDirection.LONGITUDINAL for plate in plating_zones):
                    for plate in plating_zones:
                        if plate.stiff_dir == BeamDirection.LONGITUDINAL:   # If dimension y is limited by longitudinal stiffener spacing
                            max_y = self.get_min_flange_el_length_between_psm(long1, long2)
                            dim_y = plating_mesh_dim_y[plate.id]  # Base element size in the y direction for plate in list plating_zones
                            if dim_y > max_y != 0:         # If dimension y exceeds the maximum allowed
                                dim_y = self.refine_plate_element(plate.get_stiffener_spacing() * 1000, max_y)
                                self._mesh_dim_y[i_long - 1] = dim_y        # Save value of dim_y
                            else:                                           # If dim_y does not exceed the maximum max_y
                                self._mesh_dim_y[i_long - 1] = dim_y
                            break                                           # Stop after finding the first zone with longitudinal stiffeners
                else:
                    dim_y_list = [plating_mesh_dim_y[plate.id] for plate in plating_zones]
                    dim_y = np.amin(dim_y_list)  # Use minimum value of all saved dim_y for plating zones between longitudinal psm
                    self._mesh_dim_y[i_long - 1] = dim_y

            # Assign dimension x for all plating zones between transverse primary supporting members
            for i_tran in range(1, len(self._grillage.transverse_members())):
                tran1 = self._grillage.transverse_members()[i_tran]
                tran2 = self._grillage.transverse_members()[i_tran + 1]
                plating_zones = self._grillage.plating_zones_between_psm(tran1, tran2)  # List of all plating zones between PSM

                if any(plate.stiff_dir == BeamDirection.TRANSVERSE for plate in plating_zones):
                    for plate in plating_zones:
                        if plate.stiff_dir == BeamDirection.TRANSVERSE:     # If dimension x is limited by transverse stiffener spacing
                            max_x = self.get_min_flange_el_length_between_psm(tran1, tran2)
                            dim_x = plating_mesh_dim_x[plate.id]            # Quad element size in the x direction for plate in list plating_zones
                            if dim_x > max_x != 0:         # If dimension x exceeds the maximum allowed
                                dim_x = self.refine_plate_element(plate.get_stiffener_spacing() * 1000, max_x)
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
                from the list of all x dimensions: mesh_dim_x.
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
                from the list of all y dimensions: mesh_dim_y.
        """
        if any(self._mesh_dim_y):
            segment_id = plate.trans_seg1.id
            dim_y_id = segment_id - 1
            dim_y = self._mesh_dim_y[dim_y_id]
            return dim_y
        else:
            print("ERROR: Mesh y dimensions list is blank! Calculate mesh quad element size first.")

    def element_size_transition(self, plate: Plate, segment_id):
        dim_x = self.get_base_dim_x(plate)
        dim_y = self.get_base_dim_y(plate)
        stiff_offset = plate.get_equal_stiffener_offset() * 1000    # Stiffener offset of stiffners on the plating zone

        # If stiffener direction is longitudinal, transition elements next to transverse segments do not exist
        if plate.stiff_dir == BeamDirection.LONGITUDINAL:
            plate_segments = {1: plate.long_seg1, 2: plate.long_seg2}
            flange_width = self.get_flange_el_width(plate_segments[segment_id]) * self.num_eaf
            remaining_dist = stiff_offset - flange_width    # Between stiffener offset and flange width
            n_elem = np.floor(remaining_dist / dim_y)       # Number of elements with dimension dim_y that fit inside the remaining distance
            dim_tr_y = stiff_offset - n_elem * dim_y - flange_width  # Transition element dimension x next to longitudinal segment 1

            if remaining_dist < 0:
                raise Exception("Transition element dimension y on plating zone", plate.id, "is negative!",
                                "Primary supporting member flange overlaps a stiffener, check flange width value.")
            if dim_tr_y != 0:
                ar = self.element_aspect_ratio(dim_tr_y, dim_x)
                if ar > self._plate_aspect_ratio:
                    dim_tr_y += dim_y
                return 0, dim_tr_y

        # If stiffener direction is transverse, transition elements next to longitudinal segments do not exist
        else:
            plate_segments = {1: plate.trans_seg1, 2: plate.trans_seg2}
            flange_width = self.get_flange_el_width(plate_segments[segment_id]) * self.num_eaf
            remaining_dist = stiff_offset - flange_width
            n_elem = np.floor(remaining_dist / dim_x)
            dim_tr_x = stiff_offset - n_elem * dim_x - flange_width

            if remaining_dist < 0:
                raise Exception("Transition element dimension x on plating zone", plate.id, "is negative!",
                                "Primary supporting member flange overlaps a stiffener, check flange width value.")

            if dim_tr_x != 0:
                ar = self.element_aspect_ratio(dim_tr_x, dim_y)
                if ar > self._plate_aspect_ratio:
                    dim_tr_x += dim_x
                return dim_tr_x, 0

    def transition_element_mesh(self):
        # Global consideration of transition element mesh dimensions tr_dim_x and tr_dim_y
        n_long = int(self._grillage.N_transverse - 1)     # Number of plating zones along the longitudinal axis
        n_tran = int(self._grillage.N_longitudinal - 1)   # Number of plating zones along the transverse axis
        self._tr_el_dim_x = np.zeros((2, n_long))
        self._tr_el_dim_y = np.zeros((2, n_tran))

        # Assign transition element y dimensions between pairs of longitudinal members
        for i_long in range(1, len(self._grillage.longitudinal_members())):
            long1 = self._grillage.longitudinal_members()[i_long]
            long2 = self._grillage.longitudinal_members()[i_long + 1]
            plating_zones = self._grillage.plating_zones_between_psm(long1, long2)

            for plate in plating_zones:
                tr_dim_y1 = self.element_size_transition(plate, 1)[1]
                tr_dim_y2 = self.element_size_transition(plate, 2)[1]

                if plate.stiff_dir == BeamDirection.LONGITUDINAL:
                    self._tr_el_dim_y[0, i_long - 1] = tr_dim_y1
                    self._tr_el_dim_y[1, i_long - 1] = tr_dim_y2
                    break
                else:
                    self._tr_el_dim_y[0, i_long - 1] = tr_dim_y1
                    self._tr_el_dim_y[1, i_long - 1] = tr_dim_y2

        # Assign transition element x dimensions between pairs of transverse members
        for i_tran in range(1, len(self._grillage.transverse_members())):
            tran1 = self._grillage.transverse_members()[i_tran]
            tran2 = self._grillage.transverse_members()[i_tran + 1]
            plating_zones = self._grillage.plating_zones_between_psm(tran1, tran2)

            for plate in plating_zones:
                tr_dim_x1 = self.element_size_transition(plate, 1)[0]
                tr_dim_x2 = self.element_size_transition(plate, 2)[0]

                if plate.stiff_dir == BeamDirection.TRANSVERSE:
                    self._tr_el_dim_x[0, i_tran - 1] = tr_dim_x1
                    self._tr_el_dim_x[1, i_tran - 1] = tr_dim_x2
                    break
                else:
                    self._tr_el_dim_x[0, i_tran - 1] = tr_dim_x1
                    self._tr_el_dim_x[1, i_tran - 1] = tr_dim_x2

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
        """
        tr_el_transverse = 2        # Number of transitional elements on the plating zone along the transverse segment
        tr_el_longitudinal = 2      # Number of transitional elements on the plating zone along the longitudinal segment

        if self.get_tr_dim_x(plate)[0] == 0:
            tr_el_transverse = 0
        if self.get_tr_dim_y(plate)[0] == 0:
            tr_el_longitudinal = 0

        return tr_el_transverse, tr_el_longitudinal

    def get_flange_element_num(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Number of flange elements along x and y dimension reflected onto the plating zone.
        """
        flange_element_number_x = self.num_eaf * 2
        flange_element_number_y = self.num_eaf * 2

        if self.get_flange_el_width(plate.trans_seg1) == 0:
            flange_element_number_x -= self.num_eaf
        if self.get_flange_el_width(plate.trans_seg2) == 0:
            flange_element_number_x -= self.num_eaf

        if self.get_flange_el_width(plate.long_seg1) == 0:
            flange_element_number_y -= self.num_eaf
        if self.get_flange_el_width(plate.long_seg2) == 0:
            flange_element_number_y -= self.num_eaf
        return flange_element_number_x, flange_element_number_y

    def get_element_number(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Number of all elements along x and y dimension of the plating zone.
        """

        L = plate.plate_longitudinal_dim() * 1000  # Longitudinal plating zone dimension [mm]
        B = plate.plate_transverse_dim() * 1000  # Transverse plating zone dimension [mm]

        dim_x = self.get_base_dim_x(plate)
        dim_y = self.get_base_dim_y(plate)

        tr_el_dim_x1, tr_el_dim_x2 = self.get_tr_dim_x(plate)
        tr_el_dim_y1, tr_el_dim_y2 = self.get_tr_dim_y(plate)

        fl_dim_x1 = self.get_flange_el_width(plate.trans_seg1)
        fl_dim_x2 = self.get_flange_el_width(plate.trans_seg2)
        fl_dim_y1 = self.get_flange_el_width(plate.long_seg1)
        fl_dim_y2 = self.get_flange_el_width(plate.long_seg2)

        n_tr_el_x, n_tr_el_y = self.get_tr_element_num(plate)
        n_fl_el_x, n_fl_el_y = self.get_flange_element_num(plate)

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

        tr_dim_x1, tr_dim_x2 = self.get_tr_dim_x(plate)
        fl_dim_x1 = self.get_flange_el_width(plate.trans_seg1)
        fl_dim_x2 = self.get_flange_el_width(plate.trans_seg2)
        base_dim_x = self.get_base_dim_x(plate)

        n_elem_tr = int(self.get_tr_element_num(plate)[0])
        n_elem_flange = int(self.get_flange_element_num(plate)[0])
        n_elem_x = int(self.get_element_number(plate)[0]) - n_elem_flange - n_elem_tr

        element_dim = {}    # Dictionary of all quad element dimensions along the x axis of the plating zone
        element_id = 1

        # Flange elements of transverse segment 1
        start = element_id
        end = int(start + n_elem_flange / 2)
        for element in range(start, end):
            element_dim[element_id] = fl_dim_x1
            element_id += 1

        # Transition elements next to transverse segment 1
        start = end
        end = int(start + n_elem_tr / 2)
        for element in range(start, end):
            element_dim[element_id] = tr_dim_x1
            element_id += 1

        # Base mesh elements
        start = end
        end = int(start + n_elem_x)
        for element in range(start, end):
            element_dim[element_id] = base_dim_x
            element_id += 1

        # Transition elements next to transverse segment 2
        start = end
        end = int(start + n_elem_tr / 2)
        for element in range(start, end):
            element_dim[element_id] = tr_dim_x2
            element_id += 1

        # Flange elements of transverse segment 2
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
        tr_dim_y1, tr_dim_y2 = self.get_tr_dim_y(plate)
        fl_dim_y1 = self.get_flange_el_width(plate.long_seg1)
        fl_dim_y2 = self.get_flange_el_width(plate.long_seg2)
        base_dim_y = self.get_base_dim_y(plate)

        n_elem_tr = int(self.get_tr_element_num(plate)[1])
        n_elem_flange = int(self.get_flange_element_num(plate)[1])
        n_elem_x = int(self.get_element_number(plate)[1]) - n_elem_flange - n_elem_tr

        element_dim = {}    # Dictionary of all quad element dimensions along the x axis of the plating zone
        element_id = 1

        start = element_id
        end = int(start + n_elem_flange / 2)        # If there is no flange, n_elem_flange = 0 and end = start, so the for loop skips
        for element in range(start, end):
            element_dim[element_id] = fl_dim_y1
            element_id += 1

        start = end
        end = int(start + n_elem_tr / 2)        # If there is no transition element, n_elem_tr = 0 and end = start, so the for loop skips
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

    #  *************** Ovdje bi bio kraj ove klase; Slijedeće metode raspodijeliti na novu strukturu klasa *******************

    def generate_plating_zone_elements(self, plate: Plate, start_n_id, start_e_id, split_along=AOS.NONE):
        """
        :param plate: Selected plate for calculating node locations, generating Noode and Element objects.
        :param start_n_id: Starting node ID which allows continued numeration after other methods.
        :param start_e_id: Starting element ID which allows continued numeration after other methods.
        :param split_along: Variable determines the meshing limits of the selected plating zone based on Axis Of Symmetry.
        :return: Determines node coordinates and generates finite element Node and Element objects on the selected plating zone.
                Returns last node and element ID, to continue node and element numbering on the next plating zone or segment.
        """

        print("Čvorovi zone oplate", plate.id)

        mesh_dim_x = self.get_mesh_dim_x(plate)  # Input dictionary of dimensions of all elements along x axis
        mesh_dim_y = self.get_mesh_dim_y(plate)  # Input dictionary of dimensions of all elements along y axis

        row_limit = len(mesh_dim_y) + 1  # Number of node rows on the entire plating zone
        column_limit = len(mesh_dim_x) + 1  # Number of node columns on the entire plating zone

        # **** DODATAK **** WIP
        # Ako zonu oplate presjeca os simetrije:

        # Modify row and column limits based on split_along Axis Of Symmetry input
        if split_along == AOS.LONGITUDINAL:  # Longitudinal axis of symmetry splits the plating zone
            pass
        elif split_along == AOS.TRANSVERSE:  # Transverse axis of symmetry splits the plating zone
            pass
        elif split_along == AOS.BOTH:  # Both longitudinal and transverse axis of symmetry splits the plating zone
            pass

        # Node coordinates calculation
        node_id = start_n_id
        ref_node1 = Segment.get_segment_node1(plate.long_seg1)  # Reference node 1 coordinates in [mm]
        ref_node2 = Segment.get_segment_node2(plate.long_seg1)  # Reference node 2 coordinates in [mm]
        ref_vector = np.subtract(ref_node2, ref_node1)  # Reference vector in the direction of the reference segment
        unit_ref_vector = ref_vector / np.linalg.norm(ref_vector)  # Unit reference vector
        normal_vector = np.array((0, 0, 1))  # Vector normal to the plating surface
        perpendicular_vector = np.cross(normal_vector, unit_ref_vector)  # Unit perpendicular vector

        spacing_vector_x = np.zeros(3)
        spacing_vector_y = np.zeros(3)
        dim_y_index = 1
        for row in range(0, row_limit):  # Row of nodes along x axis
            dim_x_index = 1
            if row > 0:
                spacing_vector_y += mesh_dim_y[dim_y_index] * perpendicular_vector
                dim_y_index += 1
            else:
                spacing_vector_y = np.zeros(3)

            for column in range(0, column_limit):  # Column of nodes along y axis
                if column > 0:
                    spacing_vector_x += mesh_dim_x[dim_x_index] * unit_ref_vector
                    dim_x_index += 1
                else:
                    spacing_vector_x = np.zeros(3)

                spacing_vector = spacing_vector_x + spacing_vector_y  # Node position vector in the local coordinate system
                node = spacing_vector + ref_node1  # Node coordinats in the global coordinate system, ref_node1 = position vector
                # Ovdje instanciraj objekt čvora
                print("Node ID:", node_id, ", koordinate:", node)
                node_id += 1
        last_node_id = node_id  # Return last node ID for continued node numeration

        # Node ID array - reference for quad element generation
        plate_id_array = np.zeros((row_limit, column_limit))
        node_id = start_n_id
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
        element_id = start_e_id
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

    def generate_segment_web_elements(self, segment: Segment, start_n_id, start_e_id, split=False):
        """
        :param segment: Selected segment for calculating node locations, generating Noode and Element objects.
        :param start_n_id: Starting node ID which allows continued numeration after other methods.
        :param start_e_id: Starting element ID which allows continued numeration after other methods.
        :param split: Variable determines the meshing limits of the selected segment based on Axis Of Symmetry.
        :return: Determines node coordinates and generates finite element Node and Element objects on the selected segment.
                Returns last node and element ID, to continue node and element numbering on the next segment.
        """

        print("Čvorovi struka segmenta", segment.id, "jakog nosača", segment.primary_supp_mem.id)

        direction = segment.primary_supp_mem.direction
        dim_z = self.get_web_el_height(segment)  # Vertical dimension of every segment web element

        mesh_long_dim = None
        for plate in self._grillage.plating().values():  # Identify which plating zone the segment belongs to
            segment_defines_plate = plate.test_plate_segment(segment)
            if segment_defines_plate:
                if direction == BeamDirection.LONGITUDINAL:
                    mesh_long_dim = self.get_mesh_dim_x(plate)  # Mesh dimensions along the entire length of the segment
                elif direction == BeamDirection.TRANSVERSE:
                    mesh_long_dim = self.get_mesh_dim_y(plate)
                break  # Stop after finding the first plating zone segment belongs to

        # **** DODATAK **** WIP
        # Ako segment presjeca os simetrije:

        # Modify column limits based on split input
        column_limit = len(mesh_long_dim) + 1  # Number of node columns on the entire segment
        if split:
            pass

        # Node coordinates calculation
        node_id = start_n_id
        ref_node1 = Segment.get_segment_node1(segment)  # Reference node 1 coordinates in [mm], origin of the local csy
        ref_node2 = Segment.get_segment_node2(segment)  # Reference node 2 coordinates in [mm]
        ref_vector = np.subtract(ref_node2, ref_node1)  # Reference vector in the direction of the reference segment
        unit_ref_vector = ref_vector / np.linalg.norm(ref_vector)  # Unit reference vector
        perpendicular_vector = np.array((0, 0, -1))  # Vector in the direction of psm flange, opposite of global z axis direction
        long_spacing_vector = np.zeros(3)  # Lonngitudinal spacing vector in the direction of segment psm

        for row in range(0, self._min_num_eweb + 1):  # Total number of rows of web element nodes is equal to min_num_ewb + 1
            vertical_spacing_vector = perpendicular_vector * dim_z * row
            long_dim_index = 1
            for column in range(0, column_limit):
                if column > 0:
                    long_spacing_vector += mesh_long_dim[long_dim_index] * unit_ref_vector
                    long_dim_index += 1
                else:
                    long_spacing_vector = np.zeros(3)

                position_vector = long_spacing_vector + vertical_spacing_vector  # Node position vector in the local coordinate system
                node_coords = position_vector + ref_node1
                print("Node ID:", node_id, node_coords)
                node_id += 1

        # Node ID array - reference for quad element generation
        plate_id_array = np.zeros((self._min_num_eweb + 1, column_limit))
        node_id = start_n_id
        for row in range(0, self._min_num_eweb + 1):
            for column in range(0, column_limit):
                plate_id_array[row, column] = node_id
                node_id += 1
        # print(plate_id_array)

        print(" \n Elementi segmenta", segment.id, "jakog nosača", segment.primary_supp_mem.id)

        # Generate elements
        element_id = start_e_id
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

    def generate_segment_flange_elements(self, segment: Segment, start_n_id):
        # WIP
        # return node_id, element_id
        pass

    def generate_mesh(self):
        if not any(self._mesh_dim_x) or not any(self._mesh_dim_y):
            raise Exception("ERROR: Calculate all base mesh dimensions first!")

        plate_limit = len(self._grillage.plating())  # Generate plating mesh on all plating zones
        long_member_limit = self._grillage.N_longitudinal  # Generate mesh on all longitudinal members
        tran_member_limit = self._grillage.N_transverse  # Generate mesh on all transverse members
        long_segment_limit = self._grillage.N_transverse - 1  # Generate mesh on all longitudinal segments
        tran_segment_limit = self._grillage.N_longitudinal - 1  # Generate mesh on all transverse segments

        starting_node_id = 1
        starting_element_id = 1

        # Generate all plating quad elements beam elements
        for plate_id in range(1, plate_limit + 1):
            plate = self._grillage.plating()[plate_id]
            last_IDs = self.generate_plating_zone_elements(plate, starting_node_id, starting_element_id)
            starting_node_id = last_IDs[0]
            starting_element_id = last_IDs[1]
            # self.generate_stiffener_elements(plate)

        # Generate all longitudinal segment quad elements
        for long_id in range(1, long_member_limit + 1):
            longitudinal = self._grillage.longitudinal_members()[long_id]
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
            transverse = self._grillage.transverse_members()[tran_id]
            for segment_id in range(0, tran_segment_limit):
                segment = transverse.segments[segment_id]
                last_IDs = self.generate_segment_web_elements(segment, starting_node_id, starting_element_id)
                starting_node_id = last_IDs[0]
                starting_element_id = last_IDs[1]
                # last_IDs = self.generate_segment_flange_elements(segment, starting_node_id)
                # starting_node_id = last_IDs[0]
                # starting_element_id = last_IDs[1]


# NOVA STRUKTURA KLASA
class PlatingZoneMesh:
    def __init__(self, plate: Plate, start_n_id, start_e_id, split_along=AOS.NONE):
        self._plate = plate
        self._start_node_id = start_n_id        # Starting node ID
        self._start_element_id = start_e_id     # Starting element ID
        self._split_along = split_along         # Optional argument: if axis of symmetry splits the plating zone
        self._edge_nodes_x = {}                 # Input: Distance between nodes, in order for all elements, along x axis
        self._edge_nodes_y = {}                 # Input: Distance between nodes, in order for all elements, along y axis


class PlateMesh:
    def __init__(self, grillage: Grillage):
        self._grillage = grillage
        self._axis_of_symm = MeshSize.axis_of_symm   # Input axis of symmetry for the grillage model, determines the type of FEM model


class SegmentMesh:
    def __init__(self, segment: Segment, start_n_id, start_e_id, split=False):
        self._segment = segment
        self._start_node_id = start_n_id        # Starting node ID
        self._start_element_id = start_e_id     # Starting element ID
        self._split = split                     # Optional argument: set as True if any axis of symmetry splits the segment in half
        self._edge_nodes_x = {}                 # Input: Distance between nodes, in order for all elements, along x axis
        self._edge_nodes_y = {}                 # Input: Distance between nodes, in order for all elements, along y axis


class TBeamMesh(SegmentMesh):
    def __init__(self, segment: Segment, start_n_id, start_e_id):
        super().__init__(segment, start_n_id, start_e_id)


class LBeamMesh(TBeamMesh):
    def __init__(self, segment: Segment, start_n_id, start_e_id):
        super().__init__(segment, start_n_id, start_e_id)


class FBBeamMesh(TBeamMesh):
    def __init__(self, segment: Segment, start_n_id, start_e_id):
        super().__init__(segment, start_n_id, start_e_id)


class TMeshV1(TBeamMesh):
    def __init__(self, segment: Segment, start_n_id, start_e_id):
        super().__init__(segment, start_n_id, start_e_id)


class LMeshV1(LBeamMesh):
    def __init__(self, segment: Segment, start_n_id, start_e_id):
        super().__init__(segment, start_n_id, start_e_id)


class FBMeshV1(FBBeamMesh):
    def __init__(self, segment: Segment, start_n_id, start_e_id):
        super().__init__(segment, start_n_id, start_e_id)


class PrimarySuppMemMesh:
    def __init__(self, psm: PrimarySuppMem):
        self._psm = psm
        self._axis_of_symm = MeshSize.axis_of_symm


class GrillageMesh:
    def __int__(self,  grillage: Grillage):
        self._grillage = grillage
        self._start_node_id = 1
        self._start_element_id = 1

    def generate_mesh_V1(self):
        PlateMesh(self._grillage)
        
        for member in self._grillage.longitudinal_members().values():
            for segment in member.segments:
                beam_type = segment.beam_prop.beam_type

                if beam_type == "T":
                    return TMeshV1(segment, self._start_node_id, self._start_element_id)
                elif beam_type == "L":
                    return LMeshV1(segment, self._start_node_id, self._start_element_id)
                elif beam_type == "FB":
                    return FBMeshV1(segment, self._start_node_id, self._start_element_id)
