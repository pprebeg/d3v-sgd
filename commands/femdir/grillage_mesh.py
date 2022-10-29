"""
Tool for Grillage Structure Analysis
University of Zagreb, Faculty of Mechanical Engineering and Naval Architecture
Department of Naval Architecture and Ocean Engineering

Master's thesis project

    Gordan Kos, univ.bacc.ing.nav.arch.
    Dr.sc. Pero Prebeg, dipl.ing.


MODULE FOR GRILLAGE FINITE ELEMENT MESH DEFINITION

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
    def __init__(self, grillage: Grillage, axis_of_symm=AOS.NONE):
        """
        Class for calculating mesh dimensions on the selected grillage model.

        Final result of the following methods are distances between edge nodes of all structural elements,
        along both x and y axis, which will be used for all node and element generation.

        :param grillage: Input grillage model.
        :param axis_of_symm: Optional argument: global Axis of Symmetry of the grillage model.
        """
        # super().__init__()
        self._grillage = grillage
        self._axis_of_symm = axis_of_symm

        self._min_num_ebs = 1               # Minimum number of elements between stiffeners; default = 1
        self._min_num_eweb = 3              # Minimum number of elements representing the web of a psm along its height; default = 3
        self._num_eaf = 1                   # Number of elements across primary supporting member flange; default = 1
        self._flange_aspect_ratio = 8       # Maximum aspect ratio value for primary supporting member flange quad elements; default = 8
        self._plate_aspect_ratio = 3        # Maximum aspect ratio value for plating and primary supporting member web quad elements; default = 3
        self._des_plate_aspect_ratio = 2    # Desirable plating aspect ratio value less than the maximum; default = 2
        self._unique_properties = {}        # Dictionary of unique plate thickness and material combinations used in the grillage model

        self.plating_zones_ref_array = []   # 2D array of plating zone IDs to be meshed, arranged to represent their relative placement
        self.plating_zones = {}             # All plating zones included in mesh generation
        self.full_plate_zones = {}          # Plating zones to be fully meshed
        self.long_half_plate_zones = {}     # Plating zones to be split with a longitudinal axis of symmetry
        self.tran_half_plate_zones = {}     # Plating zones to be split with a transverse axis of symmetry
        self.quarter_plate_zone = {}        # Plating zone to be split with both axis of symmetry

        self.full_long_segments = {}        # Longitudinal segments to be fully meshed
        self.mod_full_long_segments = {}    # Longitudinal segments on the longitudinal AOS, to be fully meshed with modified plate property
        self.half_long_segments = {}        # Longitudinal segments split in half by transverse AOS
        self.mod_half_long_segments = {}    # Longitudinal segments split by transverse AOS, to be fully meshed with modified plate property
        self.full_tran_segments = {}        # Transverse segments to be fully meshed
        self.mod_full_tran_segments = {}    # Transverse segments on the transverse AOS, to be fully meshed with modified plate property
        self.half_tran_segments = {}        # Transverse segments split in half by longitudinal AOS
        self.mod_half_tran_segments = {}    # Transverse segments split by longitudinal AOS, to be fully meshed with modified plate property

        self._plate_edge_node_x = {}        # Distance between all nodes along x axis, in order, for all meshed plating zones
        self._plate_edge_node_y = {}        # Distance between all nodes along y axis, in order, for all meshed plating zones

    """
    Razmisliti da ove dimenzije ne budu strogo dimenzije elemenata, već položaji čvorova duž ruba elemenata gdje čvorovi moraju biti.
    Što se događa s mrežom unutar pojedinog elementa je do odabira varijante izrade mreže tog dijela, ali rubni čvorovi moraju postojati.
    """

    @property
    def grillage(self):
        return self._grillage

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

    @property
    def plate_edge_node_x(self):
        return self._plate_edge_node_x

    @property
    def plate_edge_node_y(self):
        return self._plate_edge_node_y

    @des_plate_aspect_ratio.setter
    def des_plate_aspect_ratio(self, value):
        self._des_plate_aspect_ratio = value
        if value > self.plate_aspect_ratio:
            raise Exception("Desired plate element aspect ratio can not be greater than the maximum plate aspect ratio!")

    def hc_plate_zone_reference_ID_array(self):
        """
        :return: 2D array of all plating zone IDs arranged to represent relative placement of zones on the entire grillage model.
        """
        total_rows = self._grillage.N_longitudinal - 1        # Total number of plating zone rows on the grillage model
        total_columns = self._grillage.N_transverse - 1       # Total number of plating zone columns on the grillage model
        total_zones = total_rows * total_columns

        id_list = np.arange(1, total_zones + 1, 1)
        plating_zone_array = np.reshape(id_list, [total_rows, total_columns])
        return plating_zone_array

    # def long_segment_reference_ID_array(self):
    #     """
    #     :return: 1D array of longitudinal segment IDs arranged to represent segments along a primary supporting member.
    #     """
    #     n_tran = self._grillage.N_transverse
    #     segment_ID_array = np.arange(1, n_tran, 1)
    #     return segment_ID_array

    # def tran_segment_reference_ID_array(self):
    #     """
    #     :return: 1D array of transverse segment IDs arranged to represent segments along a primary supporting member.
    #     """
    #     n_long = self._grillage.N_longitudinal
    #     segment_ID_array = np.arange(1, n_long, 1)
    #     return segment_ID_array

    def longitudinal_psm_extent(self):
        """
        :return: Dictionary of longitudinal primary supporting members to be considered for mesh generation,
                based on input Axis of Symmetry.
        """
        longitudinals = {}

        if self.axis_of_symm == AOS.LONGITUDINAL or self.axis_of_symm == AOS.BOTH:
            n_long = np.ceil(self._grillage.N_longitudinal / 2)
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

        if self.axis_of_symm == AOS.TRANSVERSE or self.axis_of_symm == AOS.BOTH:
            n_tran = np.ceil(self._grillage.N_transverse / 2)
            for i in range(1, int(n_tran) + 1):
                transversals[i] = self._grillage.transverse_members()[i]
            return transversals
        else:
            return self._grillage.transverse_members()

    def identify_long_full_segments(self):
        """
        :return: Identifies longitudinal Segment objects for full mesh generation.
                Stores identified segments in full_long_segments dictionary.
        """
        # n_long = self._grillage.N_longitudinal              # Number of longitudinal primary supporting members
        # n_tran = self._grillage.N_transverse                # Number of transverse primary supporting members
        n_long = 5
        n_tran = 5
        n_full_segments = int(np.floor((n_tran - 1) / 2))   # Number of longitudinal segments to be fully meshed
        print("Broj punih uzdužnih segmenata", n_full_segments)

        if np.mod(n_long, 2) == 0:
            n = 1
            for member in self.longitudinal_psm_extent().values():
                for segment_id in range(0, n_full_segments):
                    print("Nosač:", member.id, ", segment:", segment_id)
                    self.full_long_segments[n] = member.segments[segment_id]
                    n += 1
        else:
            for member_id in range(1, len(self.longitudinal_psm_extent())):
                print(member_id)

    def identify_tran_full_segments(self):
        # WIP
        pass

    def identify_long_half_segments(self):
        n_tran = self._grillage.N_transverse
        n_full_segments = int(np.floor((n_tran - 1) / 2))  # Number of segmenets to be fully meshed

        if np.mod(n_tran, 2) == 0:  # Center segment is split in half by transverse axis of symmetry
            middle_segment_ID = n_full_segments
            n = 1
            for member in self._grillage.longitudinal_members().values():
                self.half_long_segments[n] = member.segments[middle_segment_ID]
                n += 1

    def identify_tran_half_segments(self):
        # WIP
        pass

    def grillage_segment_extent(self):
        """
        :return: Determines limits of segment mesh generation based on input Axis of Symmetry value.
                Calls specific methods for identifying which segments will be fully or partially meshed.
                If grillage has no axis of symmetry, all segments on the grillage model will be meshed.
        """

        # **** NA KRAJU SPOJITI OVU METODU SA grillage_plate_extent U JEDNU: grillage_mesh_extents **********

        if self.axis_of_symm == AOS.LONGITUDINAL:
            self.identify_long_full_segments()
            self.identify_tran_full_segments()
            self.identify_tran_half_segments()

        elif self.axis_of_symm == AOS.TRANSVERSE:
            self.identify_tran_full_segments()
            self.identify_long_full_segments()
            self.identify_long_half_segments()

        elif self.axis_of_symm == AOS.BOTH:
            self.identify_long_full_segments()
            self.identify_long_half_segments()
            self.identify_tran_full_segments()
            self.identify_tran_half_segments()
        else:
            self.identify_long_full_segments()
            self.identify_tran_full_segments()

    def identify_long_full_plate_zones(self):
        """
        :return: Identifies Plate objects for full mesh generation on the entire plating zone; grillage with a longitudinal axis of symm.
                Stores identified zones in full_plate_zones dictionary.
        """
        total_rows = self._grillage.N_longitudinal - 1    # Total number of plating zone rows on the grillage model
        plating_zone_array = self.hc_plate_zone_reference_ID_array()

        n_full_rows = int(np.floor(total_rows / 2))
        zone_arr_split = plating_zone_array[0:n_full_rows, :]
        for i in zone_arr_split:
            for j in i:
                self.full_plate_zones[j] = self._grillage.plating()[j]
                self.plating_zones[j] = self._grillage.plating()[j]

    def identify_tran_full_plate_zones(self):
        """
        :return: Identifies Plate objects for full mesh generation on the entire plating zone; grillage with a transverse axis of symm.
                Stores identified zones in full_plate_zones dictionary.
        """
        total_columns = self._grillage.N_transverse - 1   # Total number of plating zone columns on the grillage model
        plating_zone_array = self.hc_plate_zone_reference_ID_array()

        n_full_columns = int(np.floor(total_columns / 2))
        zone_arr_split = plating_zone_array[:, 0:n_full_columns]
        for i in zone_arr_split:
            for j in i:
                self.full_plate_zones[j] = self._grillage.plating()[j]
                self.plating_zones[j] = self._grillage.plating()[j]

    def identify_both_full_plate_zones(self):
        """
        :return: Identifies Plate objects for full mesh generation on the entire plating zone; grillage with both axis of symm.
                Stores identified zones in full_plate_zones dictionary.
        """
        total_rows = self._grillage.N_longitudinal - 1    # Total number of plating zone rows on the grillage model
        total_columns = self._grillage.N_transverse - 1   # Total number of plating zone columns on the grillage model
        plating_zone_array = self.hc_plate_zone_reference_ID_array()

        middle_row_index = int(np.floor(total_rows / 2))
        middle_column_index = int(np.floor(total_columns / 2))
        zone_arr_split = plating_zone_array[0:middle_row_index, 0:middle_column_index]
        for i in zone_arr_split:
            for j in i:
                self.full_plate_zones[j] = self._grillage.plating()[j]
                self.plating_zones[j] = self._grillage.plating()[j]

    def identify_long_half_plate_zones(self):
        """
        :return: Identifies Plate objects for half mesh generation; plating zones split with a longitudinal axis of symmetry.
                Stores identified zones in long_half_plate_zones dictionary.
        """
        total_rows = self._grillage.N_longitudinal - 1              # Total number of plating zone rows on the grillage model
        plating_zone_array = self.hc_plate_zone_reference_ID_array()

        if np.mod(self._grillage.N_longitudinal, 2) == 0:           # Odd number of plating zone rows
            middle_row_index = int(np.floor(total_rows / 2))
            middle_row = plating_zone_array[middle_row_index, :]    # Find middle row of plating zone index array
            for i in middle_row:
                self.long_half_plate_zones[i] = self._grillage.plating()[i]
                self.plating_zones[i] = self._grillage.plating()[i]

    def identify_tran_half_plate_zones(self):
        """
        :return: Identifies Plate objects for half mesh generation; plating zones split with a transverse axis of symmetry.
                Stores identified zones in tran_half_plate_zones dictionary.
        """
        total_columns = self._grillage.N_transverse - 1   # Total number of plating zone columns on the grillage model
        plating_zone_array = self.hc_plate_zone_reference_ID_array()

        if np.mod(self._grillage.N_transverse, 2) == 0:                 # Odd number of plating zone columns
            middle_column_index = int(np.floor(total_columns / 2))
            middle_column = plating_zone_array[:, middle_column_index]  # Find middle column of plating zone index array
            for i in middle_column:
                self.tran_half_plate_zones[i] = self._grillage.plating()[i]
                self.plating_zones[i] = self._grillage.plating()[i]

    def identify_both_half_plate_zones(self):
        """
        :return: Identifies Plate objects for half mesh generation; grillage with both axis of symm.
                Stores identified zones in long_half_plate_zones and tran_half_plate_zones dictionary.
        """
        total_rows = self._grillage.N_longitudinal - 1              # Total number of plating zone rows on the grillage model
        total_columns = self._grillage.N_transverse - 1   # Total number of plating zone columns on the grillage model
        plating_zone_array = self.hc_plate_zone_reference_ID_array()
        middle_row_index = int(np.floor(total_rows / 2))
        middle_column_index = int(np.floor(total_columns / 2))

        if np.mod(self._grillage.N_longitudinal, 2) == 0:           # Odd number of plating zone rows
            middle_row = plating_zone_array[middle_row_index, 0:middle_column_index]        # Find middle row of plating zone index array
            for i in middle_row:
                self.long_half_plate_zones[i] = self._grillage.plating()[i]
                self.plating_zones[i] = self._grillage.plating()[i]

        if np.mod(self._grillage.N_transverse, 2) == 0:                 # Odd number of plating zone columns
            middle_column = plating_zone_array[0:middle_row_index, middle_column_index]     # Find middle column of plating zone index array
            for i in middle_column:
                self.tran_half_plate_zones[i] = self._grillage.plating()[i]
                self.plating_zones[i] = self._grillage.plating()[i]

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
            plating_zone_array = self.hc_plate_zone_reference_ID_array()
            middle_row_index = int(np.floor(total_rows / 2))
            middle_column_index = int(np.floor(total_columns / 2))
            middle_plate_ID = plating_zone_array[middle_row_index, middle_column_index]
            self.quarter_plate_zone[middle_plate_ID] = self._grillage.plating()[middle_plate_ID]
            self.plating_zones[middle_plate_ID] = self._grillage.plating()[middle_plate_ID]

    def longitudinal_symm_plate_ref_array(self):
        """
        :return: 2D array of plating zone IDs to be meshed for longitudinal Axis of Symmetry,
                arranged to represent their relative placement on the grillage model.
        """
        total_rows = self._grillage.N_longitudinal - 1    # Total number of plating zone rows on the grillage model
        plating_zone_array = self.hc_plate_zone_reference_ID_array()
        middle_row_index = int(np.ceil(total_rows / 2))
        self.plating_zones_ref_array = plating_zone_array[0:middle_row_index, :]

    def transverse_symm_plate_ref_array(self):
        """
        :return: 2D array of plating zone IDs to be meshed for both longitudinal and transverse Axis of Symmetry,
                arranged to represent their relative placement on the grillage model.
        """
        total_columns = self._grillage.N_transverse - 1   # Total number of plating zone columns on the grillage model
        plating_zone_array = self.hc_plate_zone_reference_ID_array()
        middle_column_index = int(np.ceil(total_columns / 2))
        self.plating_zones_ref_array = plating_zone_array[:, 0:middle_column_index]

    def both_symm_plate_ref_array(self):
        """
        :return: 2D array of plating zone IDs to be meshed for both longitudinal and transverse Axis of Symmetry,
                arranged to represent their relative placement on the grillage model.
        """
        total_rows = self._grillage.N_longitudinal - 1    # Total number of plating zone rows on the grillage model
        total_columns = self._grillage.N_transverse - 1   # Total number of plating zone columns on the grillage model
        plating_zone_array = self.hc_plate_zone_reference_ID_array()
        middle_row_index = int(np.ceil(total_rows / 2))
        middle_column_index = int(np.ceil(total_columns / 2))
        self.plating_zones_ref_array = plating_zone_array[0:middle_row_index, 0:middle_column_index]

    def grillage_plate_extent(self):
        """
        :return: Determines limits of plate mesh generation based on input Axis of Symmetry value.
                Calls specific methods for identifying which plating zones will be fully or partially meshed.
                If grillage has no axis of symmetry, all plating zones inside grillage.plating() will be meshed.
        """
        if self.axis_of_symm == AOS.LONGITUDINAL:
            self.identify_long_full_plate_zones()
            self.identify_long_half_plate_zones()
            self.longitudinal_symm_plate_ref_array()

        elif self.axis_of_symm == AOS.TRANSVERSE:
            self.identify_tran_full_plate_zones()
            self.identify_tran_half_plate_zones()
            self.transverse_symm_plate_ref_array()

        elif self.axis_of_symm == AOS.BOTH:
            self.identify_both_full_plate_zones()
            self.identify_both_half_plate_zones()
            self.identify_quarter_plate_zone()
            self.both_symm_plate_ref_array()

        else:
            self.full_plate_zones = self._grillage.plating()
            self.plating_zones = self._grillage.plating()
            self.plating_zones_ref_array = self.hc_plate_zone_reference_ID_array()

    def get_plate_dim(self, plate: Plate, plate_dim):
        # OVU METODU ELIMINIRATI, KORISTITI JEDNU OD DRUGIH
        """
        :param plate: Selected plating zone.
        :param plate_dim: Full plating zone dimension.
        :return: Method returns half of the full plating zone dimension if selected plate is split by any axis of symmetry.
        """
        if plate.id in self.full_plate_zones:
            return plate_dim
        elif plate.id in self.long_half_plate_zones or plate.id in self.tran_half_plate_zones or plate.id in self.quarter_plate_zone:
            return plate_dim / 2

    def get_long_plate_dim(self, plate: Plate):
        if plate.id in self.full_plate_zones or plate.id in self.long_half_plate_zones:
            return plate.plate_longitudinal_dim() * 1000
        elif plate.id in self.tran_half_plate_zones or plate.id in self.quarter_plate_zone:
            return (plate.plate_longitudinal_dim() / 2) * 1000

    def get_tran_plate_dim(self, plate: Plate):
        if plate.id in self.full_plate_zones or plate.id in self.tran_half_plate_zones:
            return plate.plate_transverse_dim() * 1000
        elif plate.id in self.long_half_plate_zones or plate.id in self.quarter_plate_zone:
            return (plate.plate_transverse_dim() / 2) * 1000

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
        :return: Quad element dimension perpendicular to stiffener direction based only on number of elements between stiffeners.
        """
        stiff_spacing = plate.get_stiffener_spacing() * 1000    # Stiffener spacing in [mm]
        n_elem_bs = self._min_num_ebs                           # Number of elements between stiffeners
        perp_dim = stiff_spacing / n_elem_bs                    # Quad element dimension between stiffeners
        return perp_dim

    def element_size_para_to_stiffeners(self, plate: Plate, plate_dim):
        """
        :param plate: Selected plating zone.
        :param plate_dim: Plating zone dimension parallel to stiffener direction, takes into account Axis of Symmetry.
        :return: Quad element dimension parlallel to stiffener direction, based only on quad element dimension
                 perpendicular to stiffener direction, desired and maximum aspect ratio. Method attempts to find
                 a divisor of plating zone dimension parallel to stiffeners that results in quad element size not
                 greater than the desired element dimension parallel to stiffener direction. If no divisors are found
                 or the divisor is the plating zone dimension itself, method refine_plate_element is used.
        """
        # Default bi trebao biti za varijante osim V1:
        y = self.element_size_perp_to_stiffeners(plate)             # Element dimension perpendicular to stiffener direction
        des_x_val = y * self.des_plate_aspect_ratio                 # Desired element dimension parallel to stiffener direction
        n_elem = self.find_largest_divisor(plate_dim, des_x_val)            # Number of elements with desired dimension des_x_val
        if n_elem is not None:
            x = plate_dim / n_elem                                          # Element dimension parallel to stiffener direction
            ar = self.element_aspect_ratio(x, y)
            if ar > self.plate_aspect_ratio:                        # Check aspect ratio in case L is a prime number
                x = self.refine_plate_element(plate_dim, des_x_val)
        else:                                                       # There is no divisor in case length L is not a integer
            x = self.refine_plate_element(plate_dim, des_x_val)
        return x

    def get_web_el_height(self, segment: Segment):
        """
        :param segment: Segment of a primary supporting member.
        :return: Quad element dimension along the height of a primary supporting member.
        """
        hw = segment.beam_prop.hw
        dim = hw / self._min_num_eweb
        return dim

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

    def element_size_plating_zone_perp(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Quad element dimension perpendicular to stiffener direction based on two criteria:
                1.) Number of elements between stiffeners: method element_size_perp_to_stiffeners
                2.) Maximum element dimension: method get_min_flange_el_length
        """
        y = self.element_size_perp_to_stiffeners(plate)     # Element dimension perpendicular to stiffeners

        if plate.stiff_dir == BeamDirection.LONGITUDINAL:
            max_dim = self.get_min_flange_el_length(plate.trans_seg1, plate.trans_seg2)  # Maximum allowed element dimension
        else:
            max_dim = self.get_min_flange_el_length(plate.long_seg1, plate.long_seg2)    # Maximum allowed element dimension

        if y > max_dim != 0:
            y = self.refine_plate_element(plate.get_stiffener_spacing() * 1000, max_dim)
        return y

    def element_size_plating_zone_para(self, plate: Plate, plate_dim):
        """
        :param plate: Selected plating zone.
        :param plate_dim: Plating zone dimension parallel to stiffener direction.
        :return: Quad element dimension parallel to stiffener direction based on two criteria:
                1.) Element dimensiom: method element_size_para_to_stiffeners
                2.) Maximum element dimension: method get_min_flange_el_length
        """
        x = self.element_size_para_to_stiffeners(plate, plate_dim)

        if plate.stiff_dir == BeamDirection.LONGITUDINAL:
            max_dim = self.get_min_flange_el_length(plate.long_seg1, plate.long_seg2)  # Maximum allowed element dimension
        else:
            max_dim = self.get_min_flange_el_length(plate.trans_seg1, plate.trans_seg2)    # Maximum allowed element dimension

        if x > max_dim != 0:
            x = self.refine_plate_element(plate_dim, max_dim)
        return x

    def element_size_plating_zone(self, plate: Plate, plate_dim):
        """
        Method for local consideration of base mesh dimensions dim_x and dim_y, for each plating zone individually.

        :param plate: Selected plating zone.
        :param plate_dim: Plating zone dimension parallel to stiffener direction. For MeshV1 this value should be calculated using
                method get_reduced_plate_dim and for others plate_dim_parallel_to_stiffeners.
        :return: Base mesh dimensions x and y that are below the maximum values based on plate aspect ratio and flange dimensions,
                considering each plating zone individually.
        """
        if plate.stiff_dir == BeamDirection.LONGITUDINAL:
            dim_x = self.element_size_plating_zone_para(plate, plate_dim)
            dim_y = self.element_size_plating_zone_perp(plate)
        else:
            dim_x = self.element_size_plating_zone_perp(plate)
            dim_y = self.element_size_para_to_stiffeners(plate, plate_dim)

        dim_xf = self.get_min_flange_el_length(plate.long_seg1, plate.long_seg2)        # Maximum base dim_x
        dim_yf = self.get_min_flange_el_length(plate.trans_seg1, plate.trans_seg2)      # Maximum base dim_y

        if self.element_aspect_ratio(dim_x, dim_y) > self._plate_aspect_ratio:  # Check plating element aspect ratio
            dim_x_limit = np.minimum(dim_xf, dim_y * self.plate_aspect_ratio)  # Maximum allowed dim_x dimension
            dim_y_limit = np.minimum(dim_yf, dim_x * self.plate_aspect_ratio)  # Maximum allowed dim_y dimension

            if plate.stiff_dir == BeamDirection.LONGITUDINAL:
                if dim_y > dim_x:
                    dim_y = self.refine_plate_element(plate.get_stiffener_spacing() * 1000, dim_y_limit)
                else:
                    dim_x = self.refine_plate_element(plate_dim, dim_x_limit)

            elif plate.stiff_dir == BeamDirection.TRANSVERSE:
                if dim_y > dim_x:
                    dim_y = self.refine_plate_element(plate_dim, dim_y_limit)
                else:
                    dim_x = self.refine_plate_element(plate.get_stiffener_spacing() * 1000, dim_x_limit)

        return dim_x, dim_y

    def calc_element_base_size(self):
        """
        :return: Calculates the quad element size based on stiffener spacing and maximum allowed aspect ratio for all plating zones.
            Returns dictionaries of x and y base dimensions for all plating zones.
        """
        mesh_dim_x = {}  # Dimension x for all plating zones, based on element_size_plating_zone
        mesh_dim_y = {}  # Dimension y for all plating zones, based on element_size_plating_zone

        # Calculate the quad element size based on stiffener spacing and maximum allowed aspect ratio for all plating zones
        for plate in self.grillage.plating().values():
            plate_dim = self.get_plate_dim(plate, plate.plate_dim_parallel_to_stiffeners() * 1000)
            dim_x, dim_y = self.element_size_plating_zone(plate, plate_dim)
            mesh_dim_x[plate.id] = dim_x
            mesh_dim_y[plate.id] = dim_y

        return mesh_dim_x, mesh_dim_y

    def assign_base_dim_x(self, mesh_dim_x):
        """
        Method for global consideration of base mesh dimensions dim_x, for each column of plating zones on the grillage model.

        :param mesh_dim_x: Input dictionary of quad element x dimensions based on stiffener spacing
                and maximum allowed aspect ratio for all plating zones.
        :return: Assigns dimension x for each column of plating zones between transverse primary supporting members,
                identified using plating zones reference array based on Axis of Symmetry input.

        """
        ref_array = self.plating_zones_ref_array
        n_rows, n_columns = np.shape(ref_array)
        final_mesh_dim_x = {}

        for column in range(1, n_columns + 1):              # Column of plating zones between transverse primary supporting members
            plating_zone_IDs = ref_array[:, column - 1]     # List of plating zone IDs in the selected column
            plating_zones = [self.grillage.plating()[plate_id] for plate_id in plating_zone_IDs]  # List of plate objects in the selected column
            tran1 = self._grillage.transverse_members()[column]
            tran2 = self._grillage.transverse_members()[column + 1]

            if any(plate.stiff_dir == BeamDirection.TRANSVERSE for plate in plating_zones):
                for plate in plating_zones:
                    if plate.stiff_dir == BeamDirection.TRANSVERSE:  # If dimension x is limited by transverse stiffener spacing
                        max_x = self.get_min_flange_el_length_between_psm(tran1, tran2)
                        dim_x = mesh_dim_x[plate.id]  # Quad element size in the x direction for plate in list plating_zones
                        if dim_x > max_x != 0:  # If dimension x exceeds the maximum allowed
                            dim_x = self.refine_plate_element(plate.get_stiffener_spacing() * 1000, max_x)
                            final_mesh_dim_x[column] = dim_x  # Save value of dim_x
                        else:
                            final_mesh_dim_x[column] = dim_x
                        break  # Stop after finding the first zone with transverse stiffeners
            else:
                dim_x_list = [mesh_dim_x[plate.id] for plate in plating_zones]
                dim_x = np.amin(dim_x_list)  # Use minimum value of all saved dim_x for plating zones between transverse psm
                final_mesh_dim_x[column] = dim_x

        return final_mesh_dim_x

    def assign_base_dim_y(self, mesh_dim_y):
        """
        Method for global consideration of base mesh dimensions dim_y, for each row of plating zones on the grillage model.

        :param mesh_dim_y: Input dictionary of quad element y dimensions based on stiffener spacing
                and maximum allowed aspect ratio for all plating zones.
        :return: Assigns dimension y for each row of plating zones between longitudinal primary supporting members,
                identified using plating zones reference array based on Axis of Symmetry input.
        """
        ref_array = self.plating_zones_ref_array
        n_rows, n_columns = np.shape(ref_array)
        final_mesh_dim_y = {}

        for row in range(1, n_rows + 1):                # Row of plating zones between longitudinal primary supporting members
            plating_zone_IDs = ref_array[row - 1, :]    # List of plating zone IDs in the selected row
            plating_zones = [self.grillage.plating()[plate_id] for plate_id in plating_zone_IDs]  # List of plate objects in the selected row
            long1 = self._grillage.longitudinal_members()[row]
            long2 = self._grillage.longitudinal_members()[row + 1]

            if any(plate.stiff_dir == BeamDirection.LONGITUDINAL for plate in plating_zones):
                for plate in plating_zones:
                    if plate.stiff_dir == BeamDirection.LONGITUDINAL:  # If dimension y is limited by longitudinal stiffener spacing
                        max_y = self.get_min_flange_el_length_between_psm(long1, long2)
                        dim_y = mesh_dim_y[plate.id]  # Base element size in the y direction for plate in list plating_zones
                        if dim_y > max_y != 0:  # If dimension y exceeds the maximum allowed
                            dim_y = self.refine_plate_element(plate.get_stiffener_spacing() * 1000, max_y)
                            final_mesh_dim_y[row] = dim_y  # Save value of dim_y
                        else:  # If dim_y does not exceed the maximum max_y
                            final_mesh_dim_y[row] = dim_y
                        break  # Stop after finding the first zone with longitudinal stiffeners
            else:
                dim_y_list = [mesh_dim_y[plate.id] for plate in plating_zones]
                dim_y = np.amin(dim_y_list)  # Use minimum value of all saved dim_y for plating zones between longitudinal psm
                final_mesh_dim_y[row] = dim_y

        return final_mesh_dim_y

    def get_mesh_dim_x(self, plate: Plate):
        pass

    def get_mesh_dim_y(self, plate: Plate):
        pass


class MeshV1(MeshSize):
    def __init__(self, grillage: Grillage, axis_of_symm=AOS.NONE):
        """
        Class for calculating mesh dimensions specific to meshing variant V1.

        Mesh variant V1 reflects primary supporting member flange elements onto plating and has the following limitations:
            1.) Flange width of a primary supporting member has to be the same on all of its segments.
            2.) All primary supporting members need to have the same web height.
            3.) Flange element overlap has to be in the same plane.
            4.) Grillage plating can not be defined with any camber.

        :param grillage: Input grillage model.
        :param axis_of_symm: Optional argument: global Axis of Symmetry of the grillage model.

        """
        super().__init__(grillage, axis_of_symm)
        self._mesh_dim_x = {}   # Dictionary of final base mesh x dimensions (dim_x) for one column of plating zones
        self._mesh_dim_y = {}   # Dictionary of final base mesh y dimensions (dim_y) for one row of plating zones
        self._tr_el_dim_x = []  # List of transition element x dimensions
        self._tr_el_dim_y = []  # List of transition element y dimensions

    @property
    def mesh_dim_x(self):
        return self._mesh_dim_x

    @mesh_dim_x.setter
    def mesh_dim_x(self, value):
        self._mesh_dim_x = value

    @property
    def mesh_dim_y(self):
        return self._mesh_dim_y

    @mesh_dim_y.setter
    def mesh_dim_y(self, value):
        self._mesh_dim_y = value

    def get_reduced_plate_dim(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Reduced plate dimensions based on plate stiffener orientation and flange width for MeshVariantV1.
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

    def calc_element_base_size(self):
        """
        Difference in this method specific to mesh variant V1 is in plate dimension plate_dim used for calculating base mesh dimensions.
        If flange elements are reflected onto plating mesh, a reduced plate dimension should be used.

        :return: Calculates the quad element size based on stiffener spacing and maximum allowed aspect ratio for all plating zones.
            Returns dictionaries of x and y base dimensions for all plating zones.
        """
        mesh_dim_x = {}  # Dimension x for all plating zones, based on element_size_plating_zone
        mesh_dim_y = {}  # Dimension y for all plating zones, based on element_size_plating_zone

        # Calculate the quad element size based on stiffener spacing and maximum allowed aspect ratio for all plating zones
        for plate in self.plating_zones.values():
            plate_dim = self.get_plate_dim(plate, self.get_reduced_plate_dim(plate))
            dim_x, dim_y = self.element_size_plating_zone(plate, plate_dim)
            mesh_dim_x[plate.id] = dim_x
            mesh_dim_y[plate.id] = dim_y

        return mesh_dim_x, mesh_dim_y

    def element_base_size_mesh(self):
        """
        Method for global consideration of base mesh dimensions dim_x and dim_y.

        :return: Assigns base dim_x and dim_y to each row and column of plating zones that will be meshed.
                Dimensions are saved for each row and column in dictionary mesh_dim_x and mesh_dim_y.
        """
        # Global consideration of base mesh dimensions dim_x and dim_y
        if self._grillage.hc_variant_check() is True:         # Calculate element size only if grillage model passes hc_variant_check
            base_dim_x, base_dim_y = self.calc_element_base_size()
            self.mesh_dim_x = self.assign_base_dim_x(base_dim_x)
            self.mesh_dim_y = self.assign_base_dim_y(base_dim_y)

    def get_base_dim_x(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Base quad element x dimension for any plating zone. Returns the value based on longitudinal segment ID
                from the list of all x dimensions: mesh_dim_x.
        """
        if self._mesh_dim_x:
            segment_id = plate.long_seg1.id
            return self._mesh_dim_x[segment_id]
        else:
            print("ERROR: Mesh x dimensions dictionary is blank! Calculate mesh quad element size first.")

    def get_base_dim_y(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Base quad element y dimension for any plating zone. Returns the value based on transverse segment ID
                from the list of all y dimensions: mesh_dim_y.
        """
        if self._mesh_dim_y:
            segment_id = plate.trans_seg1.id
            return self._mesh_dim_y[segment_id]
        else:
            print("ERROR: Mesh y dimensions dictionary is blank! Calculate mesh quad element size first.")

    def transition_element_size_plating_zone(self, plate: Plate, segment_id):
        # Method for local consideration of transition mesh dimensions dim_tr_x and dim_tr_y, for each plating zone individually.
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

    def assign_transition_dim_x(self):
        """
        Method for global consideration of transition mesh dimension x, for each column of plating zones.

        :return: Assigns transition elemenet x dimension for each column of plating zones between transverse primary supporting members,
                identified using plating zones reference array based on Axis of Symmetry input.
        """
        ref_array = self.plating_zones_ref_array
        n_rows, n_columns = np.shape(ref_array)
        self._tr_el_dim_x = np.zeros((2, n_columns))

        for column in range(1, n_columns + 1):  # Column of plating zones between transverse primary supporting members
            plating_zone_IDs = ref_array[:, column - 1]  # List of plating zone IDs in the selected column
            plating_zones = [self.grillage.plating()[plate_id] for plate_id in plating_zone_IDs]

            for plate in plating_zones:
                tr_dim_x1 = self.transition_element_size_plating_zone(plate, 1)[0]
                tr_dim_x2 = self.transition_element_size_plating_zone(plate, 2)[0]

                if plate.stiff_dir == BeamDirection.TRANSVERSE:
                    self._tr_el_dim_x[0, column - 1] = tr_dim_x1
                    self._tr_el_dim_x[1, column - 1] = tr_dim_x2
                    break
                else:
                    self._tr_el_dim_x[0, column - 1] = tr_dim_x1
                    self._tr_el_dim_x[1, column - 1] = tr_dim_x2

    def assign_transition_dim_y(self):
        """
        Method for global consideration of transition mesh dimension y, for each row of plating zones.

        :return: Assigns transition elemenet y dimension for each row of plating zones between longitudinal primary supporting members,
                identified using plating zones reference array based on Axis of Symmetry input.
        """
        ref_array = self.plating_zones_ref_array
        n_rows, n_columns = np.shape(ref_array)
        self._tr_el_dim_y = np.zeros((2, n_rows))

        for row in range(1, n_rows + 1):  # Row of plating zones between longitudinal primary supporting members
            plating_zone_IDs = ref_array[row - 1, :]  # List of plating zone IDs in the selected row
            plating_zones = [self.grillage.plating()[plate_id] for plate_id in plating_zone_IDs]

            for plate in plating_zones:
                tr_dim_y1 = self.transition_element_size_plating_zone(plate, 1)[1]
                tr_dim_y2 = self.transition_element_size_plating_zone(plate, 2)[1]

                if plate.stiff_dir == BeamDirection.LONGITUDINAL:
                    self._tr_el_dim_y[0, row - 1] = tr_dim_y1
                    self._tr_el_dim_y[1, row - 1] = tr_dim_y2
                    break
                else:
                    self._tr_el_dim_y[0, row - 1] = tr_dim_y1
                    self._tr_el_dim_y[1, row - 1] = tr_dim_y2

    def element_transition_size_mesh(self):
        # Global consideration of transition element mesh dimensions tr_dim_x and tr_dim_y
        self.assign_transition_dim_x()
        self.assign_transition_dim_y()

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
        tr_el_longitudinal = 2      # Number of transitional elements on the plating zone in the longitudinal direction
        tr_el_transverse = 2        # Number of transitional elements on the plating zone in the transverse direction

        if plate.id in self.long_half_plate_zones:
            tr_el_transverse -= 1

        elif plate.id in self.tran_half_plate_zones:
            tr_el_longitudinal -= 1

        elif plate.id in self.quarter_plate_zone:
            tr_el_transverse -= 1
            tr_el_longitudinal -= 1

        tr_x1, tr_x2 = self.get_tr_dim_x(plate)
        if tr_x1 == 0 and tr_el_longitudinal > 0:
            tr_el_longitudinal -= 1
        if tr_x2 == 0 and tr_el_longitudinal > 0:
            tr_el_longitudinal -= 1

        tr_y1, tr_y2 = self.get_tr_dim_y(plate)
        if tr_y1 == 0 and tr_el_transverse > 0:
            tr_el_transverse -= 1
        if tr_y2 == 0 and tr_el_transverse > 0:
            tr_el_transverse -= 1

        return tr_el_longitudinal, tr_el_transverse

    def get_flange_element_num(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Number of flange elements along x and y dimension reflected onto the plating zone.
        """
        flange_element_number_x = self.num_eaf * 2
        flange_element_number_y = self.num_eaf * 2

        if plate.id in self.long_half_plate_zones:
            flange_element_number_y -= self.num_eaf

        elif plate.id in self.tran_half_plate_zones:
            flange_element_number_x -= self.num_eaf

        elif plate.id in self.quarter_plate_zone:
            flange_element_number_x -= self.num_eaf
            flange_element_number_y -= self.num_eaf

        if self.get_flange_el_width(plate.trans_seg1) == 0 and flange_element_number_x > 0:
            flange_element_number_x -= self.num_eaf
        if self.get_flange_el_width(plate.trans_seg2) == 0 and flange_element_number_x > 0:
            flange_element_number_x -= self.num_eaf

        if self.get_flange_el_width(plate.long_seg1) == 0 and flange_element_number_y > 0:
            flange_element_number_y -= self.num_eaf
        if self.get_flange_el_width(plate.long_seg2) == 0 and flange_element_number_y > 0:
            flange_element_number_y -= self.num_eaf

        return flange_element_number_x, flange_element_number_y

    def get_element_number(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Number of all elements along x and y dimension of the plating zone.
            x - number of elements in the longitudinal direction
            y - number of elements in the transverse direction
        """

        L = self.get_long_plate_dim(plate)  # Longitudinal plating zone dimension based on Axis of Symmetry input, [mm]
        B = self.get_tran_plate_dim(plate)  # Transverse plating zone dimension, based on Axis of Symmetry input, [mm]

        dim_x = self.get_base_dim_x(plate)
        dim_y = self.get_base_dim_y(plate)

        tr_el_dim_x1, tr_el_dim_x2 = self.get_tr_dim_x(plate)
        tr_el_dim_y1, tr_el_dim_y2 = self.get_tr_dim_y(plate)

        fl_dim_x1 = self.get_flange_el_width(plate.trans_seg1)
        fl_dim_x2 = self.get_flange_el_width(plate.trans_seg2)
        fl_dim_y1 = self.get_flange_el_width(plate.long_seg1)
        fl_dim_y2 = self.get_flange_el_width(plate.long_seg2)

        if plate.id in self.long_half_plate_zones:
            fl_dim_y2 = 0

        elif plate.id in self.tran_half_plate_zones:
            fl_dim_x2 = 0

        elif plate.id in self.quarter_plate_zone:
            fl_dim_x2 = 0
            fl_dim_y2 = 0

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
        :return: Distance between nodes along both longitudinal edges (along x axis), in order, for the selected plating zone.
        """

        tr_dim_x1, tr_dim_x2 = self.get_tr_dim_x(plate)
        fl_dim_x1 = self.get_flange_el_width(plate.trans_seg1)
        fl_dim_x2 = self.get_flange_el_width(plate.trans_seg2)
        base_dim_x = self.get_base_dim_x(plate)

        n_elem_tr = int(self.get_tr_element_num(plate)[0])
        n_elem_flange = int(self.get_flange_element_num(plate)[0])
        n_elem_x = int(self.get_element_number(plate)[0]) - n_elem_flange - n_elem_tr
        # print("Zona", plate.id, n_elem_x)

        element_dim = {}    # Dictionary of all quad element dimensions along the x axis of the plating zone
        element_id = 1

        # Flange elements of transverse segment 1
        if n_elem_flange == 1:
            start = element_id
            end = int(start + n_elem_flange)
            for element in range(start, end):
                element_dim[element_id] = fl_dim_x1
                element_id += 1
        else:
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
        :return: Distance between nodes along both longitudinal edges (along x axis), in order, for the selected plating zone.
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
        if n_elem_flange == 1:
            end = int(start + n_elem_flange)        # If there is no flange, n_elem_flange = 0 and end = start, so the for loop skips
            for element in range(start, end):
                element_dim[element_id] = fl_dim_y1
                element_id += 1
        else:
            end = int(start + n_elem_flange / 2)  # If there is no flange, n_elem_flange = 0 and end = start, so the for loop skips
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

    def calculate_mesh_dimensions(self):
        """
        :return: Calculates all element dimensions and stores spacing between edge nodes, along both x and y axis,
            for each plating zone to be meshed into dictionaries plate_edge_node_x and plate_edge_node_y.
            These values need to be calculated only once and will be used for all node and element generation.
        """

        # ************ Razmisliti o zaokruživanju zbog numeričke greške ***********

        self.grillage_plate_extent()            # Calculate mesh limits
        self.element_base_size_mesh()           # Calculate base mesh size
        self.element_transition_size_mesh()     # Calculate transition mesh size

        for plate in self.plating_zones.values():
            self._plate_edge_node_x[plate.id] = self.get_mesh_dim_x(plate)  # Spacing between all edge nodes along x axis
            self._plate_edge_node_y[plate.id] = self.get_mesh_dim_y(plate)  # Spacing between all edge nodes along y axis


class MeshV2(MeshSize):
    def __init__(self, grillage: Grillage, axis_of_symm=AOS.NONE):
        """
        Class for calculating mesh dimensions specific to meshing variant V2.

        WIP

        :param grillage: Input grillage model.
        :param axis_of_symm: Optional argument: global Axis of Symmetry of the grillage model.
        """
        super().__init__(grillage, axis_of_symm)
        self._mesh_dim_x = []  # List of base mesh x dimensions (dim_x) in the longitudinal direction
        self._mesh_dim_y = []  # List of base mesh y dimensions (dim_y) in the transverse direction


class PlatingZoneMesh:
    def __init__(self, mesh_size: MeshSize, plate: Plate, start_n_id, start_e_id, split_along=AOS.NONE):
        """
        Class for generating FE mesh on a selected plating zone.

        :param mesh_size: Selected input mesh size for plate zone mesh generation.
        :param plate: Selected plate for calculating node locations, generating Node and Element objects.
        :param start_n_id: Starting node ID which allows continued numeration after other methods.
        :param start_e_id: Starting element ID which allows continued numeration after other methods.
        :param split_along: Optional argument, determines the meshing limits of the selected plating zone based on Axis Of Symmetry.
        :return: Determines node coordinates and generates finite element Node and Element objects on the selected plating zone.
                Returns last node and element ID, to continue node and element numbering on the next plating zone or segment.
        """
        self._mesh_size = mesh_size
        self._plate = plate
        self._start_node_id = start_n_id
        self._start_element_id = start_e_id
        self._split_along = split_along

        self._edge_nodes_x = self._mesh_size.plate_edge_node_x[plate.id]    # Input: Distance between edge nodes, in order along x axis
        self._edge_nodes_y = self._mesh_size.plate_edge_node_y[plate.id]    # Input: Distance between edge nodes, in order along y axis

    def get_mesh_limits(self):
        """
        :return: Row and column limit values for generating nodes and elements.
        """
        row_limit = len(self._edge_nodes_y) + 1  # Number of node rows on the entire plating zone
        column_limit = len(self._edge_nodes_x) + 1  # Number of node columns on the entire plating zone

        if self._split_along == AOS.LONGITUDINAL:  # Longitudinal axis of symmetry splits the plating zone
            row_limit = int(np.floor(len(self._edge_nodes_y) / 2) + 1)

        elif self._split_along == AOS.TRANSVERSE:  # Transverse axis of symmetry splits the plating zone
            column_limit = int(np.floor(len(self._edge_nodes_x) / 2) + 1)

        elif self._split_along == AOS.BOTH:  # Both longitudinal and transverse axis of symmetry splits the plating zone
            row_limit = int(np.floor(len(self._edge_nodes_y) / 2) + 1)
            column_limit = int(np.floor(len(self._edge_nodes_x) / 2) + 1)

        return row_limit, column_limit

    def reference_node_ID_array(self, row_limit, column_limit):
        """
        :return: 2D array of node IDs arranged to represent relative placement of nodes on the plating zone.
                Used as a reference for quad element generation.
        """
        node_id_array = np.zeros((row_limit, column_limit))
        node_id = self._start_node_id
        for row in range(0, row_limit):
            for column in range(0, column_limit):
                node_id_array[row, column] = node_id
                node_id += 1
        return node_id_array

    def generate_nodes(self):
        ref_node1 = Segment.get_segment_node1(self._plate.long_seg1)        # Reference node 1 coordinates in [mm]
        ref_node2 = Segment.get_segment_node2(self._plate.long_seg1)        # Reference node 2 coordinates in [mm]
        ref_vector = np.subtract(ref_node2, ref_node1)                      # Reference vector in the direction of the reference segment
        unit_ref_vector = ref_vector / np.linalg.norm(ref_vector)           # Unit reference vector
        normal_vector = np.array((0, 0, 1))                                 # Vector normal to the plating surface
        perpendicular_vector = np.cross(normal_vector, unit_ref_vector)     # Unit perpendicular vector

        spacing_vector_x = np.zeros(3)
        spacing_vector_y = np.zeros(3)
        dim_y_index = 1
        node_id = self._start_node_id
        row_limit, column_limit = self.get_mesh_limits()
        for row in range(0, row_limit):                     # Row of nodes along x axis
            dim_x_index = 1
            if row > 0:
                spacing_vector_y += self._edge_nodes_y[dim_y_index] * perpendicular_vector
                dim_y_index += 1
            else:
                spacing_vector_y = np.zeros(3)

            for column in range(0, column_limit):           # Column of nodes along y axis
                if column > 0:
                    spacing_vector_x += self._edge_nodes_x[dim_x_index] * unit_ref_vector
                    dim_x_index += 1
                else:
                    spacing_vector_x = np.zeros(3)

                spacing_vector = spacing_vector_x + spacing_vector_y  # Node position vector in the local coordinate system
                node = spacing_vector + ref_node1  # Node coordinats in the global coordinate system, ref_node1 = position vector

                # Ovdje instanciraj objekt čvora

                print("Node ID:", node_id, ", koordinate:", node)
                node_id += 1

        return node_id

    def generate_elements(self):
        row_limit, column_limit = self.get_mesh_limits()
        node_id_array = self.reference_node_ID_array(row_limit, column_limit)
        element_id = self._start_element_id
        for row in range(0, row_limit - 1):
            for column in range(0, column_limit - 1):
                node1_id = node_id_array[row, column]
                node2_id = node_id_array[row, column + 1]
                node3_id = node_id_array[row + 1, column + 1]
                node4_id = node_id_array[row + 1, column]

                # Ovdje instanciraj objekt quad elementa, dodaj čvorove, svojstvo elementu, itd.
                # Quad elementu dodijeliti svojstvo identificirano u self._unique_properties

                print("Quad element ID", element_id, ", Node IDs:", node1_id, node2_id, node3_id, node4_id)
                element_id += 1
        return element_id

    def generate_mesh(self):
        """
        :return: Generates all nodes and elements. Returns last node and element ID to continue numeration.
        """
        nodes = self.generate_nodes()
        elements = self.generate_elements()
        return nodes, elements


class PlateMesh:
    def __init__(self, mesh_size: MeshSize):
        """
        Class for generating FE mesh on all plating zones.

        :param mesh_size: Calculated mesh dimensions.
        """
        self._mesh_size = mesh_size
        self._axis_of_symm = mesh_size.axis_of_symm
        self._grillage = mesh_size.grillage

    def generate_mesh(self):
        node_id = 1
        element_id = 1

        for plate in self._mesh_size.full_plate_zones.values():
            node, element = PlatingZoneMesh(self._mesh_size, plate, node_id, element_id).generate_mesh()
            node_id = node
            element_id = element

        for plate in self._mesh_size.long_half_plate_zones.values():
            node, element = PlatingZoneMesh(self._mesh_size, plate, node_id, element_id, AOS.LONGITUDINAL).generate_mesh()
            node_id = node
            element_id = element

        for plate in self._mesh_size.tran_half_plate_zones.values():
            node, element = PlatingZoneMesh(self._mesh_size, plate, node_id, element_id, AOS.TRANSVERSE).generate_mesh()
            node_id = node
            element_id = element

        for plate in self._mesh_size.quarter_plate_zone.values():
            node, element = PlatingZoneMesh(self._mesh_size, plate, node_id, element_id, AOS.BOTH).generate_mesh()
            node_id = node
            element_id = element


#  ************ WIP **********


class SegmentMesh:
    def __init__(self, mesh_size: MeshSize, segment: Segment, start_n_id, start_e_id, split=False):
        self._mesh_size = mesh_size
        self._segment = segment
        self._start_node_id = start_n_id        # Starting node ID
        self._start_element_id = start_e_id     # Starting element ID
        self._split = split                     # Optional argument: set as True if any axis of symmetry splits the segment in half
        self._edge_nodes_x = {}                 # Input: Distance between edge nodes, in order along x axis
        self._edge_nodes_y = {}                 # Input: Distance between edge nodes, in order along y axis
        self._edge_nodes_z = {}                 # Input: Distance between edge nodes, in order along z axis

    def generate_nodes(self):
        pass

    def generate_elements(self):
        pass

    def generate_mesh(self):
        """
        :return: Generates all nodes and elements. Returns last node and element ID to continue numeration.
        """
        # nodes = self.generate_nodes()
        # elements = self.generate_elements()
        # return nodes, elements
        pass

    # def generate_segment_web_elements(self, segment: Segment, start_n_id, start_e_id, split=False):
    #     """
    #     :param segment: Selected segment for calculating node locations, generating Noode and Element objects.
    #     :param start_n_id: Starting node ID which allows continued numeration after other methods.
    #     :param start_e_id: Starting element ID which allows continued numeration after other methods.
    #     :param split: Variable determines the meshing limits of the selected segment based on Axis Of Symmetry.
    #     :return: Determines node coordinates and generates finite element Node and Element objects on the selected segment.
    #             Returns last node and element ID, to continue node and element numbering on the next segment.
    #     """
    #
    #     print("Čvorovi struka segmenta", segment.id, "jakog nosača", segment.primary_supp_mem.id)
    #
    #     direction = segment.primary_supp_mem.direction
    #     dim_z = self.get_web_el_height(segment)  # Vertical dimension of every segment web element
    #
    #     mesh_long_dim = None
    #     for plate in self._grillage.plating().values():  # Identify which plating zone the segment belongs to
    #         segment_defines_plate = plate.test_plate_segment(segment)
    #         if segment_defines_plate:
    #             if direction == BeamDirection.LONGITUDINAL:
    #                 mesh_long_dim = self.get_mesh_dim_x(plate)  # Mesh dimensions along the entire length of the segment
    #             elif direction == BeamDirection.TRANSVERSE:
    #                 mesh_long_dim = self.get_mesh_dim_y(plate)
    #             break  # Stop after finding the first plating zone segment belongs to
    #
    #     # **** DODATAK **** WIP
    #     # Ako segment presjeca os simetrije:
    #
    #     # Modify column limits based on split input
    #     column_limit = len(mesh_long_dim) + 1  # Number of node columns on the entire segment
    #     if split:
    #         pass
    #
    #     # Node coordinates calculation
    #     node_id = start_n_id
    #     ref_node1 = Segment.get_segment_node1(segment)  # Reference node 1 coordinates in [mm], origin of the local csy
    #     ref_node2 = Segment.get_segment_node2(segment)  # Reference node 2 coordinates in [mm]
    #     ref_vector = np.subtract(ref_node2, ref_node1)  # Reference vector in the direction of the reference segment
    #     unit_ref_vector = ref_vector / np.linalg.norm(ref_vector)  # Unit reference vector
    #     perpendicular_vector = np.array((0, 0, -1))  # Vector in the direction of psm flange, opposite of global z axis direction
    #     long_spacing_vector = np.zeros(3)  # Lonngitudinal spacing vector in the direction of segment psm
    #
    #     for row in range(0, self._min_num_eweb + 1):  # Total number of rows of web element nodes is equal to min_num_ewb + 1
    #         vertical_spacing_vector = perpendicular_vector * dim_z * row
    #         long_dim_index = 1
    #         for column in range(0, column_limit):
    #             if column > 0:
    #                 long_spacing_vector += mesh_long_dim[long_dim_index] * unit_ref_vector
    #                 long_dim_index += 1
    #             else:
    #                 long_spacing_vector = np.zeros(3)
    #
    #             position_vector = long_spacing_vector + vertical_spacing_vector  # Node position vector in the local coordinate system
    #             node_coords = position_vector + ref_node1
    #             print("Node ID:", node_id, node_coords)
    #             node_id += 1
    #
    #     # Node ID array - reference for quad element generation
    #     plate_id_array = np.zeros((self._min_num_eweb + 1, column_limit))
    #     node_id = start_n_id
    #     for row in range(0, self._min_num_eweb + 1):
    #         for column in range(0, column_limit):
    #             plate_id_array[row, column] = node_id
    #             node_id += 1
    #     # print(plate_id_array)
    #
    #     print(" \n Elementi segmenta", segment.id, "jakog nosača", segment.primary_supp_mem.id)
    #
    #     # Generate elements
    #     element_id = start_e_id
    #     for row in range(0, self._min_num_eweb):
    #         for column in range(0, column_limit - 1):
    #             node1_id = plate_id_array[row, column]
    #             node2_id = plate_id_array[row, column + 1]
    #             node3_id = plate_id_array[row + 1, column + 1]
    #             node4_id = plate_id_array[row + 1, column]
    #             # Ovdje instanciraj objekt quad elementa, dodaj čvorove, svojstvo elementu, itd.
    #             # Quad elementu dodijeliti svojstvo identificirano u self._unique_properties
    #             print("Quad element ID", element_id, ", Node IDs:", node1_id, node2_id, node3_id, node4_id)
    #             element_id += 1
    #
    #     return node_id, element_id


class TBeamMesh(SegmentMesh):
    def __init__(self, mesh_size: MeshSize, segment: Segment, start_n_id, start_e_id):
        super().__init__(mesh_size, segment, start_n_id, start_e_id)


class LBeamMesh(TBeamMesh):
    def __init__(self, mesh_size: MeshSize, segment: Segment, start_n_id, start_e_id):
        super().__init__(mesh_size, segment, start_n_id, start_e_id)


class FBBeamMesh(TBeamMesh):
    def __init__(self, mesh_size: MeshSize, segment: Segment, start_n_id, start_e_id):
        super().__init__(mesh_size, segment, start_n_id, start_e_id)


class TMeshV1(TBeamMesh):
    def __init__(self, mesh_size: MeshSize, segment: Segment, start_n_id, start_e_id):
        super().__init__(mesh_size, segment, start_n_id, start_e_id)


class LMeshV1(LBeamMesh):
    def __init__(self, mesh_size: MeshSize, segment: Segment, start_n_id, start_e_id):
        super().__init__(mesh_size, segment, start_n_id, start_e_id)


class FBMeshV1(FBBeamMesh):
    def __init__(self, mesh_size: MeshSize, segment: Segment, start_n_id, start_e_id):
        super().__init__(mesh_size, segment, start_n_id, start_e_id)


class PrimarySuppMemMesh:
    def __init__(self, psm: PrimarySuppMem):
        self._psm = psm
        self._axis_of_symm = MeshSize.axis_of_symm


class GrillageMesh:
    def __int__(self,  grillage: Grillage):
        self._grillage = grillage
        self._start_node_id = 1
        self._start_element_id = 1
        # Koncept:
        """
    def generate_mesh_V1(self):
        PlateMesh(self._grillage)
        for member in self._grillage.longitudinal_members().values():
        
        # Ova petlja bi trebala biti u klasi koja generira cijeli jaki nosač:
            for segment in member.segments:
                beam_type = segment.beam_prop.beam_type

                if beam_type == "T":
                    return TMeshV1(segment, self._start_node_id, self._start_element_id)
                elif beam_type == "L":
                    return LMeshV1(segment, self._start_node_id, self._start_element_id)
                elif beam_type == "FB":
                    return FBMeshV1(segment, self._start_node_id, self._start_element_id)
        """
