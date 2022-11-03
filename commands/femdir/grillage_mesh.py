"""
Tool for Grillage Structure Analysis
University of Zagreb, Faculty of Mechanical Engineering and Naval Architecture
Department of Naval Architecture and Ocean Engineering

Master's thesis project

    Gordan Kos, univ.bacc.ing.nav.arch.
    Dr.sc. Pero Prebeg, dipl.ing.


MODULE FOR GRILLAGE FINITE ELEMENT MESH DEFINITION

"""
import itertools
from grillage.grillage_model import *
np.set_printoptions(linewidth=400)

# import femdir.geofementity as gfe
# import femdir.geofem


class ModelCheck:
    def __init__(self, grillage: Grillage):
        """
        Class for checking the grillage model for symmetry and mesh generation feasibility.

        Symmetry checks include:
            Relative distances of symmetric primary supporting members: longitudinal_psm_symmetry, transverse_psm_symmetry.
            Central primary supporting member position: central_longitudinal, central_transversal.
            Plate property, stiffener layout and stiffener direction of plating zones: longitudinal_plate_symmetry, transverse_plate_symmetry.
            Beam property of segments: longitudinal_segment_symmetry, transverse_segment_symmetry.

        :param grillage: Input grillage model.
        """
        self._grillage = grillage

    def longitudinal_psm_symmetry(self):
        """
        :return: True if relative distances of longitudinal primary supporting members are symmetric.
        """
        test = False
        for member in self._grillage.longitudinal_members().values():
            if member.rel_dist < 0.5:
                y1 = member.rel_dist
                y2 = member.symmetric_member.rel_dist
                test = np.isclose(y1, 1 - y2)
                if not test:
                    break
        return test

    def transverse_psm_symmetry(self):
        """
        :return: True if relative distances of transverse primary supporting members are symmetric.
        """
        test = False
        for member in self._grillage.transverse_members().values():
            if member.rel_dist < 0.5:
                x1 = member.rel_dist
                x2 = member.symmetric_member.rel_dist
                test = np.isclose(x1, 1 - x2)
                if not test:
                    break
        return test

    def central_longitudinal(self):
        """
        :return: True if central longitudinal primary supporting member exists and has relative distance coordinate 0.5
        """
        n_long = self._grillage.N_longitudinal
        if np.mod(n_long, 2) != 0:
            central_member_id = int(np.ceil(n_long / 2))
            central_member = self._grillage.longitudinal_members()[central_member_id]
            if np.isclose(central_member.rel_dist, 0.5):
                return True
            else:
                return False
        else:
            return True

    def central_transversal(self):
        """
        :return: True if central transverse primary supporting member exists and has relative distance coordinate 0.5
        """
        n_tran = self._grillage.N_transverse
        if np.mod(n_tran, 2) != 0:
            central_member_id = int(np.ceil(n_tran / 2))
            central_member = self._grillage.transverse_members()[central_member_id]
            if np.isclose(central_member.rel_dist, 0.5):
                return True
            else:
                return False
        else:
            return True

    def plate_zone_ID_array(self):
        """
        :return: 2D array of all plating zone IDs arranged to represent relative placement of zones on the entire grillage model.
        """
        total_rows = self._grillage.N_longitudinal - 1        # Total number of plating zone rows on the grillage model
        total_columns = self._grillage.N_transverse - 1       # Total number of plating zone columns on the grillage model
        total_zones = total_rows * total_columns

        id_list = np.arange(1, total_zones + 1, 1)
        plating_zone_array = np.reshape(id_list, [total_rows, total_columns])
        return plating_zone_array

    def longitudinal_plate_symmetry(self):
        """
        :return: True if longitudinally symmetric plating zones have the same plate property, stiffener layout and stiffener direction.
        """
        same_plate_prop = True
        same_stiff_layout = True
        same_stiff_dir = True

        plate_ref_array = self.plate_zone_ID_array()
        long_symm_ref_array = np.flip(plate_ref_array, axis=0)

        row_limit = int(np.floor((self._grillage.N_longitudinal - 1) / 2))
        column_limit = self._grillage.N_transverse - 1

        for row in range(0, row_limit):
            for column in range(0, column_limit):
                plate_id = plate_ref_array[row, column]
                symm_plate_id = long_symm_ref_array[row, column]
                plate = self._grillage.plating()[plate_id]
                symm_plate = self._grillage.plating()[symm_plate_id]

                if plate.plate_prop is not symm_plate.plate_prop:
                    same_plate_prop = False
                    break

                if plate.stiff_layout is not symm_plate.stiff_layout:
                    same_stiff_layout = False
                    break

                if plate.stiff_dir is not symm_plate.stiff_dir:
                    same_stiff_dir = False
                    break

        tests = [same_plate_prop, same_stiff_layout, same_stiff_dir]
        if all(tests):
            return True
        else:
            return False

    def transverse_plate_symmetry(self):
        """
        :return: True if transversely symmetric plating zones have the same plate property, stiffener layout and stiffener direction.
        """
        same_plate_prop = True
        same_stiff_layout = True
        same_stiff_dir = True

        plate_ref_array = self.plate_zone_ID_array()
        tran_symm_ref_array = np.flip(plate_ref_array, axis=1)

        row_limit = self._grillage.N_longitudinal - 1
        column_limit = int(np.floor((self._grillage.N_transverse - 1) / 2))

        for row in range(0, row_limit):
            for column in range(0, column_limit):
                plate_id = plate_ref_array[row, column]
                symm_plate_id = tran_symm_ref_array[row, column]
                plate = self._grillage.plating()[plate_id]
                symm_plate = self._grillage.plating()[symm_plate_id]

                if plate.plate_prop is not symm_plate.plate_prop:
                    same_plate_prop = False
                    break

                if plate.stiff_layout is not symm_plate.stiff_layout:
                    same_stiff_layout = False
                    break

                if plate.stiff_dir is not symm_plate.stiff_dir:
                    same_stiff_dir = False
                    break

        tests = [same_plate_prop, same_stiff_layout, same_stiff_dir]
        if all(tests):
            return True
        else:
            return False

    def longitudinal_segment_symmetry(self):
        """
        :return: True if symmetric longitudinal segments have the same beam property.
        """
        same_beam_prop = True

        for member in self._grillage.longitudinal_members().values():
            if member.rel_dist < 0.5:
                for segment_id in range(1, self._grillage.N_transverse):
                    segment = member.segments[segment_id - 1]
                    symm_segment = member.symmetric_member.segments[segment_id - 1]

                    if segment.beam_prop is not symm_segment.beam_prop:
                        same_beam_prop = False
                        break
        return same_beam_prop

    def transverse_segment_symmetry(self):
        """
        :return: True if symmetric transverse segments have the same beam property.
        """
        same_beam_prop = True

        for member in self._grillage.transverse_members().values():
            if member.rel_dist < 0.5:
                for segment_id in range(1, self._grillage.N_longitudinal):
                    segment = member.segments[segment_id - 1]
                    symm_segment = member.symmetric_member.segments[segment_id - 1]

                    if segment.beam_prop is not symm_segment.beam_prop:
                        same_beam_prop = False
                        break
        return same_beam_prop

    def longitudinal_symmetry_tests(self):
        """
        :return: True if model passes all longitudinal symmetry tests.
        """
        tests = [self.longitudinal_psm_symmetry(),
                 self.central_longitudinal(),
                 self.longitudinal_segment_symmetry(),
                 self.longitudinal_plate_symmetry()]
        if all(tests):
            return True
        else:
            return False

    def transverse_symmetry_tests(self):
        """
        :return: True if model passes all transverse symmetry tests.
        """
        tests = [self.transverse_psm_symmetry(),
                 self.central_transversal(),
                 self.transverse_segment_symmetry(),
                 self.transverse_plate_symmetry()]
        if all(tests):
            return True
        else:
            return False

    def assign_symmetry(self):
        self._grillage.assign_symmetric_members()

        if self.longitudinal_symmetry_tests() and self.transverse_symmetry_tests():
            return AOS.BOTH
        elif self.longitudinal_symmetry_tests():
            return AOS.LONGITUDINAL
        elif self.transverse_symmetry_tests():
            return AOS.TRANSVERSE
        else:
            return AOS.NONE

    def mesh_feasibility(self):
        """
        :return: Returns True if input grillage model can be meshed.

        Method stops base mesh dimension calculations if grillage model does not meet the following criteria:
            1.) Plating zones between two adjacent Primary Supporting Members may not have the same stiffener
                orientation and different stiffener spacing.
            2.) Insert another impossibility here
        """
        hc_check = True
        # 1.)
        # Between longitudinal primary supporting members:
        for psm_id in range(1, self._grillage.N_longitudinal):
            psm_1 = self._grillage.longitudinal_members()[psm_id]
            psm_2 = self._grillage.longitudinal_members()[psm_id + 1]
            plate_list = self._grillage.plating_zones_between_psm(psm_1, psm_2)
            plate_combinations = itertools.combinations(plate_list, 2)

            for plate in list(plate_combinations):
                plate1 = plate[0]  # First plating zone of plate combinations
                plate2 = plate[1]  # Second plating zone of plate combinations
                spacing1 = Plate.get_stiffener_spacing(plate1)
                spacing2 = Plate.get_stiffener_spacing(plate2)

                if plate1.stiff_dir == BeamDirection.LONGITUDINAL and plate2.stiff_dir == BeamDirection.LONGITUDINAL and spacing1 != spacing2:
                    print("Stiffener spacing of longitudinal stiffeners on plating zones between adjacent longitudinal primary supporting members",
                          psm_1.id, "and", psm_2.id, "do not match! \n", "     Plating zone", plate1.id, "has stiffener spacing of",
                          spacing1, "m", ", while plating zone", plate2.id, "has spacing of", spacing2, "m")

                    hc_check = False

        # Between transverse primary supporting members:
        for psm_id in range(1, self._grillage.N_transverse):
            psm_1 = self._grillage.transverse_members()[psm_id]
            psm_2 = self._grillage.transverse_members()[psm_id + 1]
            plate_list = self._grillage.plating_zones_between_psm(psm_1, psm_2)
            plate_combinations = itertools.combinations(plate_list, 2)

            for plate in list(plate_combinations):
                plate1 = plate[0]  # First plating zone of plate combinations
                plate2 = plate[1]  # Second plating zone of plate combinations
                spacing1 = Plate.get_stiffener_spacing(plate1)
                spacing2 = Plate.get_stiffener_spacing(plate2)

                if plate1.stiff_dir == BeamDirection.TRANSVERSE and plate2.stiff_dir == BeamDirection.TRANSVERSE and spacing1 != spacing2:
                    print("Stiffener spacing of transverse stiffeners on plating zones between adjacent transverse primary supporting members",
                          psm_1.id, "and", psm_2.id, "do not match! \n", "     Plating zone", plate1.id, "has stiffener spacing of",
                          spacing1, "m", ", while plating zone", plate2.id, "has spacing of", spacing2, "m")
                    hc_check = False

        return hc_check


class MeshExtent:
    def __init__(self, grillage: Grillage,  axis_of_symm_override: AOS = None):
        """
        Class for calculating FE mesh extents for the selected grillage model and Axis of Symmetry.
        Contains dictionaries of all plating zones and segments to be fully or partially meshed.

        :param grillage: Input grillage model.
        :param axis_of_symm_override: Optional argument: overrides automatic Axis of Symmetry discovery.
        """
        self._grillage = grillage
        self._axis_of_symm_override = axis_of_symm_override
        self._aos_input = ModelCheck(self._grillage).assign_symmetry()

        self.plating_zones_ref_array = []   # 2D array of plating zone IDs to be meshed, based on Axis of Symmetry selection
        self.plating_zones = {}             # All plating zones included in mesh generation
        self.full_plate_zones = {}          # Plating zones to be fully meshed
        self.long_half_plate_zones = {}     # Plating zones to be split with a longitudinal axis of symmetry
        self.tran_half_plate_zones = {}     # Plating zones to be split with a transverse axis of symmetry
        self.quarter_plate_zone = {}        # Plating zone to be split with both axis of symmetry
        self.long_e_split_zone = {}         # Plating zones with longitudinal axis of symmetry passing between stiffeners
        self.tran_e_split_zone = {}         # Plating zones with transverse axis of symmetry passing between stiffeners

        self.full_segments = {}             # Segments to be fully meshed
        self.half_segments = {}             # Segments split in half by some axis of symmetry

    @property
    def grillage(self):
        return self._grillage

    def feasibility_test(self):
        return ModelCheck(self._grillage).mesh_feasibility()

    @property
    def axis_of_symm_override(self):
        return self._axis_of_symm_override

    @property
    def aos_input(self):
        return self._aos_input

    @property
    def axis_of_symm(self):
        if self._axis_of_symm_override:
            return self._axis_of_symm_override
        else:
            return self._aos_input

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

    def longitudinal_psm_extent(self):
        """
        :return: Dictionary of longitudinal primary supporting members to be considered for mesh generation,
                based on input Axis of Symmetry.
        """
        longitudinals = {}

        if self.axis_of_symm == AOS.LONGITUDINAL or self.axis_of_symm == AOS.BOTH:
            n_long = int(np.ceil(self._grillage.N_longitudinal / 2))
            for i in range(1, n_long + 1):
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
            n_tran = int(np.ceil(self._grillage.N_transverse / 2))
            for i in range(1, n_tran + 1):
                transversals[i] = self._grillage.transverse_members()[i]
            return transversals
        else:
            return self._grillage.transverse_members()

    def identify_long_full_segments(self):
        """
        :return: Identifies Segment objects for full mesh generation; grillage with a longitudinal axis of symmetry.
                Stores identified segments in full_segments dictionary.
        """
        n_long = self._grillage.N_longitudinal              # Number of longitudinal primary supporting members
        n_tran = self._grillage.N_transverse                # Number of transverse primary supporting members
        n = 1
        for member in self.longitudinal_psm_extent().values():
            for segment_id in range(0, n_tran - 1):
                self.full_segments[n] = member.segments[segment_id]
                n += 1

        n_tran_segments = int(np.floor((n_long - 1) / 2))   # Number of transverse segmenets to be fully meshed
        for member in self._grillage.transverse_members().values():
            for segment_id in range(0, n_tran_segments):
                self.full_segments[n] = member.segments[segment_id]
                n += 1

    def identify_tran_full_segments(self):
        """
        :return: Identifies Segment objects for full mesh generation; grillage with a transverse axis of symmetry.
                Stores identified segments in full_segments dictionary.
        """
        n_long = self._grillage.N_longitudinal              # Number of longitudinal primary supporting members
        n_tran = self._grillage.N_transverse                # Number of transverse primary supporting members
        n = 1
        n_long_segments = int(np.floor((n_tran - 1) / 2))           # Number of longitudinal segmenets to be fully meshed
        for member in self._grillage.longitudinal_members().values():
            for segment_id in range(0, n_long_segments):
                self.full_segments[n] = member.segments[segment_id]
                n += 1

        for member in self.transverse_psm_extent().values():
            for segment_id in range(0, n_long - 1):
                self.full_segments[n] = member.segments[segment_id]
                n += 1

    def identify_both_full_segments(self):
        """
        :return: Identifies Segment objects for full mesh generation; grillage with both axis of symmetry.
                Stores identified segments in full_segments dictionary.
        """
        n_long = self._grillage.N_longitudinal  # Number of longitudinal primary supporting members
        n_tran = self._grillage.N_transverse    # Number of transverse primary supporting members

        n = 1
        n_long_segments = int(np.floor((n_tran - 1) / 2))           # Number of longitudinal segmenets to be fully meshed
        for member in self.longitudinal_psm_extent().values():
            for segment_id in range(0, n_long_segments):
                self.full_segments[n] = member.segments[segment_id]
                n += 1

        n_tran_segments = int(np.floor((n_long - 1) / 2))           # Number of transverse segmenets to be fully meshed
        for member in self.transverse_psm_extent().values():
            for segment_id in range(0, n_tran_segments):
                self.full_segments[n] = member.segments[segment_id]
                n += 1

    def identify_none_full_segments(self):
        """
        :return: Identifies Segment objects for full mesh generation; grillage with no axis of symmetry.
                Stores identified segments in full_segments dictionary.
        """
        n_long = self._grillage.N_longitudinal  # Number of longitudinal primary supporting members
        n_tran = self._grillage.N_transverse    # Number of transverse primary supporting members

        n = 1
        for member in self._grillage.longitudinal_members().values():
            for segment_id in range(0, n_tran - 1):
                self.full_segments[n] = member.segments[segment_id]
                n += 1

        for member in self._grillage.transverse_members().values():
            for segment_id in range(0, n_long - 1):
                self.full_segments[n] = member.segments[segment_id]
                n += 1

    def identify_long_half_segments(self):
        """
        :return: Identifies Segment objects for half mesh generation; grillage with a longitudinal axis of symmetry.
                Stores identified segments in half_segments dictionary.
        """
        n_long = self._grillage.N_longitudinal                  # Number of longitudinal primary supporting members
        middle_segment_ID = int(np.ceil((n_long - 1) / 2))      # Middle transverse segment ID
        n = 1

        if np.mod(n_long, 2) == 0:  # Center transverse segment exists for even number of longitudinal members
            for member in self._grillage.transverse_members().values():
                self.half_segments[n] = member.segments[middle_segment_ID - 1]
                n += 1

    def identify_tran_half_segments(self):
        """
        :return: Identifies Segment objects for half mesh generation; grillage with a transverse axis of symmetry.
                Stores identified segments in half_segments dictionary.
        """
        n_tran = self._grillage.N_transverse                    # Number of transverse primary supporting members
        middle_segment_ID = int(np.ceil((n_tran - 1) / 2))      # Middle longitudinal segment ID
        n = 1

        if np.mod(n_tran, 2) == 0:    # Center longitudinal segment exists for even number of transverse members
            for member in self._grillage.longitudinal_members().values():
                self.half_segments[n] = member.segments[middle_segment_ID - 1]
                n += 1

    def identify_both_half_segments(self):
        """
        :return: Identifies Segment objects for half mesh generation; grillage with both axis of symmetry.
                Stores identified segments in half_segments dictionary.
        """
        n_long = self._grillage.N_longitudinal              # Number of longitudinal primary supporting members
        n_tran = self._grillage.N_transverse                # Number of transverse primary supporting members
        long_segment_ID = int(np.ceil((n_tran - 1) / 2))    # Middle longitudinal segment ID
        tran_segment_ID = int(np.ceil((n_long - 1) / 2))    # Middle transverse segment ID

        n = 1
        if np.mod(n_tran, 2) == 0:      # Center longitudinal segment exists for even number of transverse members
            for member in self.longitudinal_psm_extent().values():
                self.half_segments[n] = member.segments[long_segment_ID - 1]
                n += 1

        if np.mod(n_long, 2) == 0:      # Center transverse segment exists for even number of longitudinal members
            for member in self.transverse_psm_extent().values():
                self.half_segments[n] = member.segments[tran_segment_ID - 1]
                n += 1

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

    def get_plate_dim(self, plate: Plate, plate_dim):
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
        """
        :param plate: Selected plating zone.
        :return: Method returns longitudinal dimension of the selected plating zone.
                Returns half of the original value if the zone is split by transverse axis of symmetry.
        """
        if plate.id in self.full_plate_zones or plate.id in self.long_half_plate_zones:
            return plate.plate_longitudinal_dim() * 1000
        elif plate.id in self.tran_half_plate_zones or plate.id in self.quarter_plate_zone:
            return (plate.plate_longitudinal_dim() / 2) * 1000

    def get_tran_plate_dim(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Method returns transverse dimension of the selected plating zone.
                Returns half of the original value if the zone is split by longitudinal axis of symmetry.
        """
        if plate.id in self.full_plate_zones or plate.id in self.tran_half_plate_zones:
            return plate.plate_transverse_dim() * 1000
        elif plate.id in self.long_half_plate_zones or plate.id in self.quarter_plate_zone:
            return (plate.plate_transverse_dim() / 2) * 1000

    def aos_between_stiffeners(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: True if some axis of symmetry passes between stiffeners on the selected plating zone.
                False return does not indicate stiffener is located on some axis of symmetry!
                If a new stiffener layout definition type would be created with nonsymmetric stiffener placement, this method needs to be modified.
        """
        stiff_num = plate.get_stiffener_number()
        if (plate.id in self.long_half_plate_zones or plate.id in self.quarter_plate_zone) and\
                np.mod(stiff_num, 2) == 0 and plate.stiff_dir == BeamDirection.LONGITUDINAL:
            return True

        elif (plate.id in self.tran_half_plate_zones or plate.id in self.quarter_plate_zone) and\
                np.mod(stiff_num, 2) == 0 and plate.stiff_dir == BeamDirection.TRANSVERSE:
            return True
        else:
            return False

    def aos_on_stiffener(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: True if some axis of symmetry is located on and parallel to a stiffener on the selected plating zone.
                False return does not indicate axis of symmetry passes between stiffeners!
                If a new stiffener layout definition type would be created with nonsymmetric stiffener placement, this method needs to be modified.
        """
        stiff_num = plate.get_stiffener_number()
        if (plate.id in self.long_half_plate_zones or plate.id in self.quarter_plate_zone) and\
                np.mod(stiff_num, 2) != 0 and plate.stiff_dir == BeamDirection.LONGITUDINAL:
            return True

        elif (plate.id in self.tran_half_plate_zones or plate.id in self.quarter_plate_zone) and\
                np.mod(stiff_num, 2) != 0 and plate.stiff_dir == BeamDirection.TRANSVERSE:
            return True
        else:
            return False

    def aos_on_segment(self, segment: Segment):
        """
        :param segment: Segment of a primary supporting member.
        :return: True if some axis of symmetry is located on and parallel to a segment. Used to identify which segments
                should be assigned half of their original web thickness and be modelled with half of their flange.
        """
        rel_dist = segment.primary_supp_mem.rel_dist
        direction = segment.primary_supp_mem.direction
        if (self.axis_of_symm == AOS.LONGITUDINAL or self.axis_of_symm == AOS.BOTH)\
                and direction == BeamDirection.LONGITUDINAL and np.isclose(rel_dist, 0.5):
            return True
        elif (self.axis_of_symm == AOS.TRANSVERSE or self.axis_of_symm == AOS.BOTH)\
                and direction == BeamDirection.TRANSVERSE and np.isclose(rel_dist, 0.5):
            return True
        else:
            return False

    def identify_tran_split_elements(self):
        """
        :return: Identifies Plate objects where transverse axis of symmetry passes between stiffeners.
                Stores identified zones in tran_e_split_zone dictionary.
        """
        plating_zone_array = self.hc_plate_zone_reference_ID_array()
        for plate in self.tran_half_plate_zones.values():
            if self.aos_between_stiffeners(plate):
                column = np.where(plating_zone_array == plate.id)[1]    # Column index of the identified plating zone with element split
                split_element_zones = plating_zone_array[:, column]
                for i in split_element_zones:
                    for j in i:
                        self.tran_e_split_zone[j] = self._grillage.plating()[j]

    def identify_long_split_elements(self):
        """
        :return: Identifies Plate objects where longitudinal axis of symmetry passes between stiffeners.
                Stores identified zones in long_e_split_zone dictionary.
        """
        plating_zone_array = self.hc_plate_zone_reference_ID_array()
        for plate in self.long_half_plate_zones.values():
            if self.aos_between_stiffeners(plate):
                row = np.where(plating_zone_array == plate.id)[0]       # Row index of the identified plating zone with element split
                split_element_zones = plating_zone_array[row, :]
                for i in split_element_zones:
                    for j in i:
                        self.long_e_split_zone[j] = self._grillage.plating()[j]

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
        :return: 2D array of plating zone IDs to be meshed for transverse Axis of Symmetry,
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
                Saves a reference ID array of all plating zones to be meshed into plating_zones_ref_array.
        """
        if self.axis_of_symm == AOS.LONGITUDINAL:
            self.identify_long_full_plate_zones()
            self.identify_long_half_plate_zones()
            self.longitudinal_symm_plate_ref_array()
            self.identify_long_split_elements()

        elif self.axis_of_symm == AOS.TRANSVERSE:
            self.identify_tran_full_plate_zones()
            self.identify_tran_half_plate_zones()
            self.transverse_symm_plate_ref_array()
            self.identify_tran_split_elements()

        elif self.axis_of_symm == AOS.BOTH:
            self.identify_both_full_plate_zones()
            self.identify_both_half_plate_zones()
            self.identify_quarter_plate_zone()
            self.both_symm_plate_ref_array()
            self.identify_long_split_elements()
            self.identify_tran_split_elements()

        else:
            self.full_plate_zones = self._grillage.plating()
            self.plating_zones = self._grillage.plating()
            self.plating_zones_ref_array = self.hc_plate_zone_reference_ID_array()

    def grillage_segment_extent(self):
        """
        :return: Determines limits of segment mesh generation based on selected Axis of Symmetry.
                Calls specific methods for identifying which segments will be fully or partially meshed.
                If grillage has no axis of symmetry, all segments on the grillage model will be meshed.
        """
        if self.axis_of_symm == AOS.LONGITUDINAL:
            self.identify_long_full_segments()
            self.identify_long_half_segments()

        elif self.axis_of_symm == AOS.TRANSVERSE:
            self.identify_tran_full_segments()
            self.identify_tran_half_segments()

        elif self.axis_of_symm == AOS.BOTH:
            self.identify_both_full_segments()
            self.identify_both_half_segments()

        else:
            self.identify_none_full_segments()
        """
        # Provjera punih:
        print("Izrada pune mreže na segmentima jakih nosača:")
        for segment in self.full_segments.values():
            psm_id = segment.primary_supp_mem.id
            direction = segment.primary_supp_mem.direction
            print("Direction:", direction, ", PSM ID:", psm_id, ", segment ID:", segment.id)

        # Provjera polovicnih:
        print("Izrada polovične mreže na segmentima jakih nosača:")
        for segment in self.half_segments.values():
            psm_id = segment.primary_supp_mem.id
            direction = segment.primary_supp_mem.direction
            print("Direction:", direction, ", PSM ID:", psm_id, ", segment ID:", segment.id)
        """

    def grillage_mesh_extent(self):
        """
        :return: Determines limits of grillage mesh generation for plating zones and segments based on Axis of Symmetry value.
                Method stops mesh generation if the model does not pass the feasibility test.
        """
        if self.feasibility_test() is True:
            self.grillage_plate_extent()
            self.grillage_segment_extent()
        else:
            raise Exception("ERROR: Grillage model does not pass mesh feasibility test!")


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
    def __init__(self, mesh_extent: MeshExtent):
        """
        Class for calculating mesh dimensions on the selected grillage model.

        Final result of the following methods are distances between edge nodes of all structural elements,
        along both x and y axis, which will be used for all node and element generation. Dimensions for plate
        elementes are stored in plate_edge_node_x, plate_edge_node_y. Dimensions for segment flange elements
        are stored in flange_edge_node.

        :param mesh_extent: FE mesh extents for the selected grillage model and Axis of Symmetry.
        """
        # super().__init__()
        self._mesh_extent = mesh_extent

        self._grillage = self._mesh_extent.grillage
        self._axis_of_symm = self._mesh_extent.axis_of_symm

        self._min_num_ebs = 1               # Minimum number of elements between stiffeners; default = 1
        self._min_num_eweb = 3              # Minimum number of elements representing the web of a psm along its height; default = 3
        self._num_eaf = 1                   # Number of elements across primary supporting member flange; default = 1
        self._flange_aspect_ratio = 8       # Maximum aspect ratio value for primary supporting member flange quad elements; default = 8
        self._plate_aspect_ratio = 3        # Maximum aspect ratio value for plating and primary supporting member web quad elements; default = 3
        self._des_plate_aspect_ratio = 2    # Desirable plating aspect ratio value less than the maximum; default = 2
        self._unique_properties = {}        # Dictionary of unique plate thickness and material combinations used in the grillage model

        self._mesh_dim_x = {}               # Final base mesh x dimensions (dim_x) for each column of plating zones
        self._mesh_dim_y = {}               # Final base mesh y dimensions (dim_y) for each row of plating zones
        self._transition_dim_x = []         # 2D array of transition element x dimensions
        self._transition_dim_y = []         # 2D array of transition element y dimensions

        self._plate_edge_node_x = {}        # Distance between plating nodes along x axis, in order, for all meshed plating zones
        self._plate_edge_node_y = {}        # Distance between plating nodes along y axis, in order, for all meshed plating zones
        self._flange_edge_node = {}         # Distance between flange nodes, in order, for all meshed segments

    @property
    def mesh_extent(self):
        return self._mesh_extent

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

    @property
    def transition_dim_x(self):
        return self._transition_dim_x

    @transition_dim_x.setter
    def transition_dim_x(self, value):
        self._transition_dim_x = value

    @property
    def transition_dim_y(self):
        return self._transition_dim_y

    @transition_dim_y.setter
    def transition_dim_y(self, value):
        self._transition_dim_y = value

    @property
    def plate_edge_node_x(self):
        return self._plate_edge_node_x

    @plate_edge_node_x.setter
    def plate_edge_node_x(self, value):
        self._plate_edge_node_x = value

    @property
    def plate_edge_node_y(self):
        return self._plate_edge_node_y

    @plate_edge_node_y.setter
    def plate_edge_node_y(self, value):
        self._plate_edge_node_y = value

    @property
    def flange_edge_node(self):
        return self._flange_edge_node

    @flange_edge_node.setter
    def flange_edge_node(self, value):
        self._flange_edge_node = value

    @des_plate_aspect_ratio.setter
    def des_plate_aspect_ratio(self, value):
        self._des_plate_aspect_ratio = value
        if value > self.plate_aspect_ratio:
            raise Exception("Desired plate element aspect ratio can not be greater than the maximum plate aspect ratio!")

    @property
    def unique_properties(self):
        return self._unique_properties

    def create_unique_plate_property(self, plate_property: PlateProperty, thickness, material: MaterialProperty):
        """
        :param plate_property: Plate property used in the grillage model.
        :param thickness: Quad element plate property thickness to be added to unique properties.
        :param material: Quad element plate property material to be added to unique properties.
        :return: Checks for a duplicate of given thickness and material combination in unique_properties.
                If no duplicate exists, a new UniquePlateProperty is created, with a reference to PlateProperty used.
        """
        duplicate = False
        if self._unique_properties:                         # If the dictionary is not empty
            for unique_property in self._unique_properties.values():   # Check if values in dict have the same properties
                if unique_property.tp == thickness and unique_property.mat == material:       # Duplicates have the same plate thickness and material
                    duplicate = True
                    unique_property.plate_prop.append(plate_property)           # Save duplicate plate object to the list
                    break
                else:
                    duplicate = False

        if duplicate is False:                       # Create unique property if current plate is not a duplicate
            upp_id = len(self.unique_properties) + 1
            curr_upp = UniquePlateProperty(upp_id, thickness, material)
            curr_upp.plate_prop.append(plate_property)
            self._unique_properties[upp_id] = curr_upp

    def create_unique_beam_property(self, beam_property: BeamProperty, thickness, material: MaterialProperty):
        """
        :param beam_property: Beam property used in the grillage model.
        :param thickness: Quad element plate property thickness to be added to unique properties.
        :param material: Quad element plate property material to be added to unique properties.
        :return: Checks for a duplicate of given thickness and material combination in unique_properties.
                If no duplicate exists, a new UniquePlateProperty is created, with a reference to BeamProperty used.
        """
        duplicate = False
        for unique_property in self._unique_properties.values():           # Check if values in dict have the same properties as psm web
            if unique_property.tp == thickness and unique_property.mat == material:
                duplicate = True
                unique_property.beam_prop.append(beam_property)            # Save duplicate beam property object to the list
                break
            else:
                duplicate = False

        if duplicate is False:                                  # Create unique property if beam property is not a duplicate
            upp_id = len(self.unique_properties) + 1
            curr_upp = UniquePlateProperty(upp_id, thickness, material)
            curr_upp.beam_prop.append(beam_property)
            self._unique_properties[upp_id] = curr_upp

    def identify_unique_property(self):
        """
        :return: Identifies unique plate thickness and material combinations used in the grillage model for
                plating zones and segments. Separately checks net thickness and material of segmenet web and flange.
        """
        for plate_property in self._grillage.plate_props().values():
            tp = plate_property.tp_net(self._grillage.corrosion_addition()[1], plate_property.tp)
            mat = plate_property.plate_mat
            self.create_unique_plate_property(plate_property, tp, mat)

        for beam_property in self._grillage.beam_props().values():
            if beam_property.beam_type == "T" or beam_property.beam_type == "L":
                tw = beam_property.tw_net(self._grillage.corrosion_addition()[1])  # Net web thickness
                tf = beam_property.tf_net(self._grillage.corrosion_addition()[1])  # Net flange thickness
                mat = beam_property.mat
                self.create_unique_beam_property(beam_property, tw, mat)
                self.create_unique_beam_property(beam_property, tf, mat)

            if beam_property.beam_type == "FB":
                tw = beam_property.tw_net(self._grillage.corrosion_addition()[1])  # Net web thickness
                mat = beam_property.mat
                self.create_unique_beam_property(beam_property, tw, mat)

        for segment in self._mesh_extent.full_segments.values():
            if self._mesh_extent.aos_on_segment(segment):
                mat = segment.beam_prop.mat
                tw = segment.beam_prop.tw_net(self._grillage.corrosion_addition()[1]) / 2
                self.create_unique_beam_property(segment.beam_prop, tw, mat)

        for segment in self._mesh_extent.half_segments.values():
            if self._mesh_extent.aos_on_segment(segment):
                mat = segment.beam_prop.mat
                tw = segment.beam_prop.tw_net(self._grillage.corrosion_addition()[1]) / 2
                self.create_unique_beam_property(segment.beam_prop, tw, mat)

    @staticmethod
    def save_node_spacing_to_dict(dictionary: dict, n_dims: int, dimension):
        """
        Method used for saving edge node spacing dimensions (quad element dimensions) to a dictionary,
        for any number of dimensions of equal value.

        :param dictionary: Dictionary for all node spacing (element) dimensions on a plating zone.
        :param n_dims: Number of dimensions to be stored in the Dictionary.
        :param dimension: Element dimension to be stored.
        """
        if dictionary:                          # If the dictionary is not empty
            dimension_id = len(dictionary) + 1
        else:
            dimension_id = 1

        end = int(dimension_id + n_dims)
        for element in range(dimension_id, end):
            dictionary[dimension_id] = dimension
            dimension_id += 1

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

        for plate in self._mesh_extent.plating_zones.values():
            plate_dim = self._mesh_extent.get_plate_dim(plate, plate.plate_dim_parallel_to_stiffeners() * 1000)
            dim_x, dim_y = self.element_size_plating_zone(plate, plate_dim)
            mesh_dim_x[plate.id] = dim_x
            mesh_dim_y[plate.id] = dim_y

        return mesh_dim_x, mesh_dim_y

    def calc_element_transition_size_mesh(self):
        pass

    def assign_base_dim_x(self, mesh_dim_x):
        """
        Method for global consideration of base mesh dimensions dim_x, for each column of plating zones on the grillage model.

        :param mesh_dim_x: Input dictionary of quad element x dimensions based on stiffener spacing
                and maximum allowed aspect ratio for all plating zones.
        :return: Assigns dimension x for each column of plating zones between transverse primary supporting members,
                identified using plating zones reference array based on Axis of Symmetry input.

        """
        ref_array = self._mesh_extent.plating_zones_ref_array
        n_rows, n_columns = np.shape(ref_array)
        final_mesh_dim_x = {}

        for column in range(1, n_columns + 1):              # Column of plating zones between transverse primary supporting members
            plating_zone_IDs = ref_array[:, column - 1]     # List of plating zone IDs in the selected column
            plating_zones = [self._grillage.plating()[plate_id] for plate_id in plating_zone_IDs]  # List of plate objects in the selected column

            for plate in plating_zones:
                if plate.stiff_dir == BeamDirection.TRANSVERSE:  # If dimension x is limited by transverse stiffener spacing
                    tran1 = self._grillage.transverse_members()[column]
                    tran2 = self._grillage.transverse_members()[column + 1]
                    max_x = self.get_min_flange_el_length_between_psm(tran1, tran2)
                    dim_x = mesh_dim_x[plate.id]  # Quad element size in the x direction for plate in list plating_zones
                    if dim_x > max_x != 0:  # If dimension x exceeds the maximum allowed
                        dim_x = self.refine_plate_element(plate.get_stiffener_spacing() * 1000, max_x)
                        final_mesh_dim_x[column] = dim_x  # Save value of dim_x
                    else:
                        final_mesh_dim_x[column] = dim_x
                    break  # Stop after finding the first zone with transverse

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
        ref_array = self._mesh_extent.plating_zones_ref_array
        n_rows, n_columns = np.shape(ref_array)
        final_mesh_dim_y = {}

        for row in range(1, n_rows + 1):                # Row of plating zones between longitudinal primary supporting members
            plating_zone_IDs = ref_array[row - 1, :]    # List of plating zone IDs in the selected row
            plating_zones = [self._grillage.plating()[plate_id] for plate_id in plating_zone_IDs]  # List of plate objects in the selected row

            for plate in plating_zones:
                if plate.stiff_dir == BeamDirection.LONGITUDINAL:  # If dimension y is limited by longitudinal stiffener spacing
                    long1 = self._grillage.longitudinal_members()[row]
                    long2 = self._grillage.longitudinal_members()[row + 1]
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

    def calc_element_base_size_mesh(self):
        """
        Method for global consideration of base mesh dimensions dim_x and dim_y.

        :return: Assigns base dim_x and dim_y to each row and column of plating zones that will be meshed.
                Dimensions are saved for each row and column in dictionary mesh_dim_x and mesh_dim_y.
        """
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

    def get_tr_dim_x(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Transition quad element x dimensions for any plating zone. Returns the value based on longitudinal segment ID.
                First value (dim_1) represents the element closest to the first transverse segment (trans_seg_1) and the second
                (dim_2) represents the element closest to the second transverse segment (trans_seg_2) that define the plating zone.
        """
        segment_id = plate.long_seg1.id
        dim_id = segment_id - 1
        dim = self.transition_dim_x
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
        dim = self.transition_dim_y
        dim_1 = dim[0][dim_id]
        dim_2 = dim[1][dim_id]
        return dim_1, dim_2

    def get_base_element_number(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Number of base size elements along x and y dimension of the plating zone.
            x - number of elements in the longitudinal direction
            y - number of elements in the transverse direction
        """
        L = self._mesh_extent.get_long_plate_dim(plate)  # Longitudinal plating zone dimension based on Axis of Symmetry input, [mm]
        B = self._mesh_extent.get_tran_plate_dim(plate)  # Transverse plating zone dimension, based on Axis of Symmetry input, [mm]

        dim_x = self.get_base_dim_x(plate)      # Base mesh dimension x
        dim_y = self.get_base_dim_y(plate)      # Base mesh dimension y

        tr_el_dim_x1, tr_el_dim_x2 = self.get_tr_dim_x(plate)   # Transition element x dimensions
        tr_el_dim_y1, tr_el_dim_y2 = self.get_tr_dim_y(plate)   # Transition element y dimensions

        fl_dim_x1 = self.get_flange_el_width(plate.trans_seg1)
        fl_dim_x2 = self.get_flange_el_width(plate.trans_seg2)
        fl_dim_y1 = self.get_flange_el_width(plate.long_seg1)
        fl_dim_y2 = self.get_flange_el_width(plate.long_seg2)

        if plate.id in self._mesh_extent.long_half_plate_zones:
            fl_dim_y2 = 0
            tr_el_dim_y2 = 0

        elif plate.id in self._mesh_extent.tran_half_plate_zones:
            fl_dim_x2 = 0
            tr_el_dim_x2 = 0

        elif plate.id in self._mesh_extent.quarter_plate_zone:
            fl_dim_x2 = 0
            tr_el_dim_x2 = 0
            fl_dim_y2 = 0
            tr_el_dim_y2 = 0

        n_dim_x = np.floor((L - tr_el_dim_x1 - tr_el_dim_x2 - fl_dim_x1 - fl_dim_x2) / dim_x)    # Number of elements with dim_x along x axis
        n_dim_y = np.floor((B - tr_el_dim_y1 - tr_el_dim_y2 - fl_dim_y1 - fl_dim_y2) / dim_y)    # Number of elements with dim_y along y axis
        return n_dim_x, n_dim_y

    def get_flange_element_num(self, segment: Segment):
        """
        :param segment: Selected segment of a primary supporting member.
        :return: Number of elements across the flange. Method returns 0 for FB beam type.
        """
        flange_element_width = self.get_flange_el_width(segment)
        if flange_element_width == 0:
            return 0
        else:
            return self.num_eaf

    def get_long_split_element_num(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Number of base size elements being split in half by longitudinal Axis of Symmetry.
        """
        if plate.id in self._mesh_extent.long_e_split_zone:
            B = self._mesh_extent.get_tran_plate_dim(plate)
            fl_dim_y1 = self.get_flange_el_width(plate.long_seg1)
            tr_el_dim_y1 = self.get_tr_dim_y(plate)[0]
            dim_y = self.get_base_dim_y(plate)
            n_dim_y = self.get_base_element_number(plate)[1]
            remaining_dist = B - fl_dim_y1 - tr_el_dim_y1 - dim_y * n_dim_y
            if remaining_dist == 0:
                return 0
            elif remaining_dist < dim_y:
                return 1
            else:
                return 0
        else:
            return 0

    def get_tran_split_element_num(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Number of base size elements being split in half by transverse Axis of Symmetry.
        """
        if plate.id in self._mesh_extent.tran_e_split_zone:
            L = self._mesh_extent.get_long_plate_dim(plate)
            fl_dim_x1 = self.get_flange_el_width(plate.trans_seg1)
            tr_el_dim_x1 = self.get_tr_dim_x(plate)[0]
            dim_x = self.get_base_dim_x(plate)
            n_dim_x = self.get_base_element_number(plate)[0]
            remaining_dist = L - fl_dim_x1 - tr_el_dim_x1 - dim_x * n_dim_x
            if remaining_dist == 0:
                return 0
            elif remaining_dist < dim_x:
                return 1
            else:
                return 0
        else:
            return 0

    def get_longitudinal_transition_element_num(self, plate: Plate, segment: Segment):
        """
        :param plate: Selected plating zone.
        :param segment: Selected segment of a primary supporting member.
        :return: Number of transition elements between longitudinal segment flange and base elements.
        """
        transition_element_dim_1, transition_element_dim_2 = self.get_tr_dim_y(plate)
        if segment is plate.long_seg1 and transition_element_dim_1 != 0:
            return 1
        elif segment is plate.long_seg2 and transition_element_dim_2 != 0:
            return 1
        else:
            return 0

    def get_transverse_transition_element_num(self, plate: Plate, segment: Segment):
        """
        :param plate: Selected plating zone.
        :param segment: Selected segment of a primary supporting member.
        :return: Number of transition elements between transverse segment flange and base elements.
        """
        transition_element_dim_1, transition_element_dim_2 = self.get_tr_dim_x(plate)
        if segment is plate.trans_seg1 and transition_element_dim_1 != 0:
            return 1
        elif segment is plate.trans_seg2 and transition_element_dim_2 != 0:
            return 1
        else:
            return 0

    def edge_node_spacing_x(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Saves distance between plate nodes along longitudinal edges (along x axis), in order, for the selected plating zone.
        """
        fl_dim_x1 = self.get_flange_el_width(plate.trans_seg1)
        fl_dim_x2 = self.get_flange_el_width(plate.trans_seg2)
        base_dim_x = self.get_base_dim_x(plate)
        transition_element_dim_x1, transition_element_dim_x2 = self.get_tr_dim_x(plate)

        transition_element_num_seg_1 = int(self.get_transverse_transition_element_num(plate, plate.trans_seg1))
        transition_element_num_seg_2 = int(self.get_transverse_transition_element_num(plate, plate.trans_seg2))
        flange_element_num_tran_seg1 = int(self.get_flange_element_num(plate.trans_seg1))
        flange_element_num_tran_seg2 = int(self.get_flange_element_num(plate.trans_seg2))
        base_mesh_element_num = int(self.get_base_element_number(plate)[0])
        split_element_number = int(self.get_tran_split_element_num(plate))

        if plate.id in self._mesh_extent.tran_half_plate_zones or plate.id in self._mesh_extent.quarter_plate_zone:
            flange_element_num_tran_seg2 = 0
            transition_element_num_seg_2 = 0

        x_spacing = {}    # Dictionary of all quad element dimensions along the x axis of the selected plating zone

        self.save_node_spacing_to_dict(x_spacing, flange_element_num_tran_seg1, fl_dim_x1)
        self.save_node_spacing_to_dict(x_spacing, transition_element_num_seg_1, transition_element_dim_x1)
        self.save_node_spacing_to_dict(x_spacing, base_mesh_element_num, base_dim_x)
        self.save_node_spacing_to_dict(x_spacing, split_element_number, base_dim_x / 2)
        self.save_node_spacing_to_dict(x_spacing, transition_element_num_seg_2, transition_element_dim_x2)
        self.save_node_spacing_to_dict(x_spacing, flange_element_num_tran_seg2, fl_dim_x2)

        self.plate_edge_node_x[plate.id] = x_spacing  # Save to edge node dimension x of the selected plating zone

    def edge_node_spacing_y(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Saves distance between plate nodes along transverse edges (along y axis), in order, for the selected plating zone.
        """
        fl_dim_y1 = self.get_flange_el_width(plate.long_seg1)
        fl_dim_y2 = self.get_flange_el_width(plate.long_seg2)
        base_dim_y = self.get_base_dim_y(plate)
        transition_element_dim_y1, transition_element_dim_y2 = self.get_tr_dim_y(plate)

        transition_element_num_seg_1 = int(self.get_longitudinal_transition_element_num(plate, plate.long_seg1))
        transition_element_num_seg_2 = int(self.get_longitudinal_transition_element_num(plate, plate.long_seg2))
        flange_element_num_long_seg1 = int(self.get_flange_element_num(plate.long_seg1))
        flange_element_num_long_seg2 = int(self.get_flange_element_num(plate.long_seg2))
        base_mesh_element_num = int(self.get_base_element_number(plate)[1])
        split_element_number = int(self.get_long_split_element_num(plate))

        if plate.id in self._mesh_extent.long_half_plate_zones or plate.id in self._mesh_extent.quarter_plate_zone:
            flange_element_num_long_seg2 = 0
            transition_element_num_seg_2 = 0

        y_spacing = {}    # Dictionary of all quad element dimensions along the y axis of the plating zone

        self.save_node_spacing_to_dict(y_spacing, flange_element_num_long_seg1, fl_dim_y1)
        self.save_node_spacing_to_dict(y_spacing, transition_element_num_seg_1, transition_element_dim_y1)
        self.save_node_spacing_to_dict(y_spacing, base_mesh_element_num, base_dim_y)
        self.save_node_spacing_to_dict(y_spacing, split_element_number, base_dim_y / 2)
        self.save_node_spacing_to_dict(y_spacing, transition_element_num_seg_2, transition_element_dim_y2)
        self.save_node_spacing_to_dict(y_spacing, flange_element_num_long_seg2, fl_dim_y2)

        self.plate_edge_node_y[plate.id] = y_spacing

    def calculate_mesh_dimensions(self):
        """
        :return: Calculates all element dimensions and stores spacing between edge nodes, along both x and y axis,
            for each plating zone to be meshed into dictionaries plate_edge_node_x and plate_edge_node_y.
            These values need to be calculated only once and will be used for all node and element generation.
        """

        if self._mesh_extent.axis_of_symm_override:
            print("Selected Axis of Symmetry override:", self._mesh_extent.axis_of_symm_override)
            print("Automatic symmetry discovery would have selected:", self._mesh_extent.aos_input)

        else:
            print("Automatically discovered grillage model symmetry:", self._mesh_extent.aos_input)

        self._mesh_extent.grillage_mesh_extent()                # Calculate mesh limits
        self.calc_element_base_size_mesh()                      # Calculate base mesh size
        self.calc_element_transition_size_mesh()                # Calculate transition mesh size

        for plate in self._mesh_extent.plating_zones.values():
            self.edge_node_spacing_x(plate)                     # Spacing between all edge nodes along x axis
            self.edge_node_spacing_y(plate)                     # Spacing between all edge nodes along y axis


class MeshV1(MeshSize):
    def __init__(self, mesh_extent: MeshExtent):
        """
        Class for calculating mesh dimensions specific to meshing variant V1.

        Mesh variant V1 reflects primary supporting member flange elements onto plating and has the following limitations:
            1.) Flange width of a primary supporting member has to be the same on all of its segments.
            2.) All primary supporting members need to have the same web height.
            3.) Flange element overlap has to be in the same plane.
            4.) Grillage plating can not be defined with any camber.
        """
        super().__init__(mesh_extent)

    def get_reduced_plate_dim(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Reduced plate dimensions based on plate stiffener orientation and flange width for mesh variant V1.
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
        for plate in self._mesh_extent.plating_zones.values():
            plate_dim = self._mesh_extent.get_plate_dim(plate, self.get_reduced_plate_dim(plate))
            dim_x, dim_y = self.element_size_plating_zone(plate, plate_dim)
            mesh_dim_x[plate.id] = dim_x
            mesh_dim_y[plate.id] = dim_y

        return mesh_dim_x, mesh_dim_y

    def transition_element_size_plating_zone(self, plate: Plate, segment_id):
        """
        Method for local consideration of transition mesh dimensions dim_tr_x and dim_tr_y, for each plating zone individually.
        Specific to meshing variant V1

        :param plate: Selected plating zone.
        :param segment_id: 1 selects first logitudinal segment, 2 selects second longitudinal segment
        """
        dim_x = self.get_base_dim_x(plate)
        dim_y = self.get_base_dim_y(plate)
        stiff_offset = plate.get_equal_stiffener_offset() * 1000    # Stiffener offset of stiffners on the plating zone

        # If stiffener direction is longitudinal, transition elements next to transverse segments do not exist
        if plate.stiff_dir == BeamDirection.LONGITUDINAL:
            plate_segments = {1: plate.long_seg1, 2: plate.long_seg2}
            flange_width = self.get_flange_el_width(plate_segments[segment_id]) * self.num_eaf
            remaining_dist = stiff_offset - flange_width    # Between stiffener offset and flange width
            n_elem = np.floor(remaining_dist / dim_y)       # Number of elements with dimension dim_y that fit inside the remaining distance

            if n_elem == 0:                                 # If no base dimension y elements fit in the remaining distance
                dim_tr_y = stiff_offset - flange_width      # Transition element still exists
            else:
                dim_tr_y = stiff_offset - n_elem * dim_y - flange_width  # Transition element dimension x next to longitudinal segment 1

            if remaining_dist < 0:
                raise Exception("Transition element dimension y on plating zone", plate.id, "is negative!",
                                "Primary supporting member flange overlaps a stiffener, check flange width value.")
            if dim_tr_y != 0:
                ar = self.element_aspect_ratio(dim_tr_y, dim_x)
                if ar > self._plate_aspect_ratio and remaining_dist > dim_y:
                    dim_tr_y += dim_y
                    return 0, dim_tr_y
                else:
                    return 0, dim_tr_y

        # If stiffener direction is transverse, transition elements next to longitudinal segments do not exist
        else:
            plate_segments = {1: plate.trans_seg1, 2: plate.trans_seg2}
            flange_width = self.get_flange_el_width(plate_segments[segment_id]) * self.num_eaf
            remaining_dist = stiff_offset - flange_width
            n_elem = np.floor(remaining_dist / dim_x)

            if n_elem == 0:                                     # If no base dimension x elements fit in the remaining distance
                dim_tr_x = stiff_offset - flange_width          # Transition element still exists
            else:
                dim_tr_x = stiff_offset - n_elem * dim_x - flange_width

            if remaining_dist < 0:
                raise Exception("Transition element dimension x on plating zone", plate.id, "is negative!",
                                "Primary supporting member flange overlaps a stiffener, check flange width value.")

            if dim_tr_x != 0:
                ar = self.element_aspect_ratio(dim_tr_x, dim_y)
                if ar > self._plate_aspect_ratio and remaining_dist > dim_x:
                    dim_tr_x += dim_x
                    return dim_tr_x, 0
                else:
                    return dim_tr_x, 0

    def assign_transition_dim_x(self):
        """
        Method for global consideration of transition mesh dimension x, for each column of plating zones.

        :return: Assigns transition elemenet x dimension for each column of plating zones between transverse primary supporting members,
                identified using plating zones reference array based on Axis of Symmetry input.
        """
        ref_array = self._mesh_extent.plating_zones_ref_array
        n_rows, n_columns = np.shape(ref_array)
        tr_el_dim_x = np.zeros((2, n_columns))

        for column in range(1, n_columns + 1):  # Column of plating zones between transverse primary supporting members
            plating_zone_IDs = ref_array[:, column - 1]  # List of plating zone IDs in the selected column
            plating_zones = [self._grillage.plating()[plate_id] for plate_id in plating_zone_IDs]

            for plate in plating_zones:
                tr_dim_x1 = self.transition_element_size_plating_zone(plate, 1)[0]
                tr_dim_x2 = self.transition_element_size_plating_zone(plate, 2)[0]

                if plate.stiff_dir == BeamDirection.TRANSVERSE:
                    tr_el_dim_x[0, column - 1] = tr_dim_x1
                    tr_el_dim_x[1, column - 1] = tr_dim_x2
                    break
                else:
                    tr_el_dim_x[0, column - 1] = tr_dim_x1
                    tr_el_dim_x[1, column - 1] = tr_dim_x2

        self.transition_dim_x = tr_el_dim_x

    def assign_transition_dim_y(self):
        """
        Method for global consideration of transition mesh dimension y, for each row of plating zones.

        :return: Assigns transition elemenet y dimension for each row of plating zones between longitudinal primary supporting members,
                identified using plating zones reference array based on Axis of Symmetry input.
        """
        ref_array = self._mesh_extent.plating_zones_ref_array
        n_rows, n_columns = np.shape(ref_array)
        tr_el_dim_y = np.zeros((2, n_rows))

        for row in range(1, n_rows + 1):  # Row of plating zones between longitudinal primary supporting members
            plating_zone_IDs = ref_array[row - 1, :]  # List of plating zone IDs in the selected row
            plating_zones = [self._grillage.plating()[plate_id] for plate_id in plating_zone_IDs]

            for plate in plating_zones:
                tr_dim_y1 = self.transition_element_size_plating_zone(plate, 1)[1]
                tr_dim_y2 = self.transition_element_size_plating_zone(plate, 2)[1]

                if plate.stiff_dir == BeamDirection.LONGITUDINAL:
                    tr_el_dim_y[0, row - 1] = tr_dim_y1
                    tr_el_dim_y[1, row - 1] = tr_dim_y2
                    break
                else:
                    tr_el_dim_y[0, row - 1] = tr_dim_y1
                    tr_el_dim_y[1, row - 1] = tr_dim_y2

        self.transition_dim_y = tr_el_dim_y

    def calc_element_transition_size_mesh(self):
        """
        Method for global consideration of transition element mesh dimensions specific to mesh variant V1
        :return: Saves calculated transition element dimension to transition_dim_x, transition_dim_y.
        """
        self.assign_transition_dim_x()
        self.assign_transition_dim_y()


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

        self._edge_nodes_x = self._mesh_size.plate_edge_node_x[plate.id]    # Distance between edge nodes, in order along x axis
        self._edge_nodes_y = self._mesh_size.plate_edge_node_y[plate.id]    # Distance between edge nodes, in order along y axis

    def get_element_property(self):
        mat = self._plate.plate_prop.plate_mat
        tp = self._plate.plate_prop.tp
        for unique_property in self._mesh_size.unique_properties.values():
            for plate_property in unique_property.plate_prop:
                if plate_property.plate_mat is mat and plate_property.tp is tp:
                    return unique_property.id

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
        total_nodes = row_limit * column_limit
        id_list = np.arange(self._start_node_id, self._start_node_id + total_nodes, 1)
        node_id_array = np.reshape(id_list, [row_limit, column_limit])
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

        print("Čvorovi zone oplate", self._plate.id)

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
                node_coords = spacing_vector + ref_node1  # Node coordinats in the global coordinate system, ref_node1 = position vector

                # Ovdje instanciraj objekt čvora

                print(" Node ID:", node_id, ", koordinate:", node_coords)
                node_id += 1

        return node_id

    def generate_elements(self):
        row_limit, column_limit = self.get_mesh_limits()
        node_id_array = self.reference_node_ID_array(row_limit, column_limit)
        element_id = self._start_element_id

        print("Elementi zone oplate", self._plate.id)

        for row in range(0, row_limit - 1):
            for column in range(0, column_limit - 1):
                node1_id = node_id_array[row, column]
                node2_id = node_id_array[row, column + 1]
                node3_id = node_id_array[row + 1, column + 1]
                node4_id = node_id_array[row + 1, column]

                # Ovdje instanciraj objekt quad elementa, dodaj čvorove, svojstvo elementu, itd.
                # Quad elementu dodijeliti svojstvo identificirano u self._unique_properties

                print(" Quad element ID", element_id, ", Node IDs:", node1_id, node2_id, node3_id, node4_id)
                element_id += 1
        return element_id

    def generate_mesh(self):
        """
        :return: Generates all nodes and elements on the selected plating zone. Returns last node and element ID to continue numeration.
        """
        nodes = self.generate_nodes()
        elements = self.generate_elements()
        return nodes, elements


class SegmentMesh:
    def __init__(self, mesh_size: MeshSize, segment: Segment, start_n_id, start_e_id, split=False):
        """
        Class for generating FE mesh on a selected segment.

        :param mesh_size: Selected input mesh size for segment mesh generation.
        :param segment: Selected segment for calculating node locations, generating Noode and Element objects.
        :param start_n_id: Starting node ID which allows continued numeration after other methods.
        :param start_e_id: Starting element ID which allows continued numeration after other methods.
        :param split: Variable determines the meshing limits of the selected segment based on Axis of Symmetry.
        :return: Determines node coordinates and generates finite element Node and Element objects on the selected segment.
                Returns last node and element ID, to continue node and element numbering on the next segment.
        """

        self._mesh_size = mesh_size
        self._mesh_extent = self._mesh_size.mesh_extent
        self._segment = segment
        self._start_node_id = start_n_id        # Starting node ID
        self._start_element_id = start_e_id     # Starting element ID
        self._split = split                     # Optional argument: set as True if any axis of symmetry splits the segment in half

        self._edge_plate_nodes = {}     # Distances between nodes, at the web connection with plate
        self._edge_flange_nodes = {}    # Distances between nodes, at the web connection with flange
        self._edge_nodes_z = {}         # Distances between nodes, in order along z axis

        self._web_node_ref_array = np.zeros((1, 1))
        self._flange_node_ref_array = np.zeros((1, 1))

    def get_plate_edge_nodes(self):
        """
        :return: Identifies a plating zone the segment defines and gets distances between edge nodes in the appropriate direction,
                at the connection of segment web and plating zone.
        """
        direction = self._segment.primary_supp_mem.direction
        for plate in self._mesh_extent.plating_zones.values():            # Identify which plating zone the segment belongs to
            segment_defines_plate = plate.test_plate_segment(self._segment)
            if segment_defines_plate:
                if direction == BeamDirection.LONGITUDINAL:
                    self._edge_plate_nodes = self._mesh_size.plate_edge_node_x[plate.id]
                elif direction == BeamDirection.TRANSVERSE:
                    self._edge_plate_nodes = self._mesh_size.plate_edge_node_y[plate.id]
                break                                                   # Stop after finding the first plating zone the segment defines

    def node_column_limit(self):
        """
        :return: Column limit value along the length of the segment for generating nodes and elements.
        """
        if self._split is True:                 # Limit for half mesh generation - axis of symmetry splits the segment in half
            column_limit = int(np.floor(len(self._edge_plate_nodes) / 2) + 1)
        else:                                   # Limit for full mesh generation
            column_limit = int(np.floor(len(self._edge_plate_nodes)) + 1)
        return column_limit

    def get_web_element_property(self):
        """
        :return: Quad element plate property ID used for primary supporting member segment web elements.
        """
        grillage = self._segment.primary_supp_mem.grillage
        mat = self._segment.beam_prop.mat

        if self._mesh_extent.aos_on_segment(self._segment):
            tw = self._segment.beam_prop.tw_net(grillage.corrosion_addition()[1]) / 2
            for unique_property in self._mesh_size.unique_properties.values():
                if unique_property.mat is mat and unique_property.tp == tw:
                    return unique_property.id

        else:
            tw = self._segment.beam_prop.tw_net(grillage.corrosion_addition()[1])
            for unique_property in self._mesh_size.unique_properties.values():
                if unique_property.mat is mat and unique_property.tp == tw:
                    return unique_property.id

    def get_flange_element_property(self):
        """
        :return: Quad element plate property ID used for primary supporting member segment flange elements.
        """
        pass

    def reference_web_node_ID_array(self):
        pass

    def generate_web_nodes(self):
        pass

    def generate_web_elements(self):
        pass

    def reference_T_flange_node_ID_array(self):
        """
        :return: 2D array of node IDs arranged to represent relative placement of nodes on the primary supporting member segment flange
                of T beam type. Used as a reference for quad element generation. Method uses common nodes at the connection of web and flange.
        """
        web_node_id_array = self._web_node_ref_array
        last_web_node_id = web_node_id_array[-1, -1]        # Last node ID on the segment web
        last_web_node_row = web_node_id_array[-1, :]        # List of common node IDs at the connection of segment web and flange

        column_limit = self.node_column_limit()              # Number of nodes along the local longitudinal axis of the segment
        row_limit = int(self._mesh_size.num_eaf * 2) + 1    # Number of nodes in the direction of flange width
        total_nodes = column_limit * (row_limit - 1)        # Total number of flange nodes, excluding the middle row

        id_list = np.arange(last_web_node_id + 1, last_web_node_id + total_nodes + 1, 1)
        id_array = np.reshape(id_list, [row_limit - 1, column_limit])
        flange_node_id_array = np.insert(id_array, self._mesh_size.num_eaf, last_web_node_row, axis=0)  # Insert last web node row in the middle
        """
        # PROVJERA
        print("Čvorovi struka:")
        print(web_node_id_array)

        print("Zadnji red čvorova struka:")
        print(last_web_node_row)

        print("Čvorovi prirubnice prije ubacivanja zadnjeg reda:")
        print(id_array)

        print("Čvorovi prirubnice nakon ubacivanja:")
        print(flange_node_id_array)
        """
        self._flange_node_ref_array = flange_node_id_array

    def reference_L_flange_node_ID_array(self):
        """
        :return: 2D array of node IDs arranged to represent relative placement of nodes on the primary supporting member segment flange
                of L beam type. Used as a reference for quad element generation. Method uses common nodes at the connection of web and flange.
        """
        web_node_id_array = self._web_node_ref_array
        last_web_node_id = web_node_id_array[-1, -1]        # Last node ID on the segment web
        last_web_node_row = web_node_id_array[-1, :]        # List of common node IDs at the connection of segment web and flange

        column_limit = self.node_column_limit()              # Number of nodes along the local longitudinal axis of the segment
        row_limit = int(self._mesh_size.num_eaf) + 1        # Number of nodes in the direction of flange width
        total_nodes = column_limit * (row_limit - 1)        # Total number of flange nodes, excluding the middle row

        id_list = np.arange(last_web_node_id + 1, last_web_node_id + total_nodes + 1, 1)
        id_array = np.reshape(id_list, [row_limit - 1, column_limit])
        flange_node_id_array = np.insert(id_array, 0, last_web_node_row, axis=0)  # Insert last web node row at the start

        """
        # PROVJERA
        print("Čvorovi struka:")
        print(web_node_id_array)

        print("Zajednički zadnji red čvorova struka:")
        print(last_web_node_row)

        print("Čvorovi prirubnice L profila prije ubacivanja zadnjeg reda:")
        print(id_array)

        print("Čvorovi prirubnice L profila nakon ubacivanja:")
        print(flange_node_id_array)
        """
        self._flange_node_ref_array = flange_node_id_array

    def generate_flange_nodes(self):
        pass

    def generate_flange_elements(self, start_element_id):
        pass

    def generate_mesh(self):
        """
        :return: Generates all nodes and elements. Returns last node and element ID to continue numeration.
        """
        nodes = self._start_node_id
        elements = self._start_element_id

        self.get_plate_edge_nodes()
        beam_type = self._segment.beam_prop.beam_type
        self.reference_web_node_ID_array()

        if beam_type == "T":
            self.generate_web_nodes()
            last_web_element = self.generate_web_elements()
            self.reference_T_flange_node_ID_array()
            nodes = self.generate_flange_nodes()
            elements = self.generate_flange_elements(last_web_element)

        elif beam_type == "L":
            self.generate_web_nodes()
            last_web_element = self.generate_web_elements()
            self.reference_L_flange_node_ID_array()
            nodes = self.generate_flange_nodes()
            elements = self.generate_flange_elements(last_web_element)

        elif beam_type == "FB":
            nodes = self.generate_web_nodes()
            elements = self.generate_web_elements()

        return nodes, elements


class SegmentV1(SegmentMesh):
    def __init__(self, mesh_size: MeshSize, segment: Segment, start_n_id, start_e_id, split=False):
        """
        CLass for segment mesh generation specific to meshing variant V1.
        """
        super().__init__(mesh_size, segment, start_n_id, start_e_id, split)
        self._mesh_size = mesh_size
        self._segment = segment
        self._start_node_id = start_n_id        # Starting node ID
        self._start_element_id = start_e_id     # Starting element ID
        self._split = split

    def reference_web_node_ID_array(self):
        """
        :return: 2D array of node IDs arranged to represent relative placement of nodes on the primary supporting member segment web.
                Used as a reference for quad element generation.
                Version 1 assumes equal number of nodes in each row and quad elements with edges parallel to the global coordinate axis.
        """
        column_limit = self.node_column_limit()                      # Number of nodes along the local longitudinal axis of the segment
        row_limit = self._mesh_size.min_num_eweb + 1                # Number of nodes along the z axis
        total_nodes = row_limit * column_limit
        id_list = np.arange(self._start_node_id, self._start_node_id + total_nodes, 1)
        web_node_id_array = np.reshape(id_list, [row_limit, column_limit])
        self._web_node_ref_array = web_node_id_array

    def generate_web_nodes(self):
        column_limit = self.node_column_limit()                          # Number of nodes along the local longitudinal axis of the segment
        dim_z = self._mesh_size.get_web_el_height(self._segment)        # Vertical dimension of every segment web element
        mesh_dim = self._edge_plate_nodes

        node_id = self._start_node_id
        ref_node1 = Segment.get_segment_node1(self._segment)            # Reference node 1 coordinates in [mm], origin of the local csy
        ref_node2 = Segment.get_segment_node2(self._segment)            # Reference node 2 coordinates in [mm]
        ref_vector = np.subtract(ref_node2, ref_node1)                  # Reference vector in the direction of the reference segment
        unit_ref_vector = ref_vector / np.linalg.norm(ref_vector)       # Unit reference vector
        perpendicular_vector = np.array((0, 0, -1))                     # Vector along segment web, opposite of global z axis direction
        long_spacing_vector = np.zeros(3)                               # Lonngitudinal spacing vector in the direction of segment psm

        print("Čvorovi struka segmenta", self._segment.id, "jakog", self._segment.beam_prop.beam_type, "nosača",
              self._segment.primary_supp_mem.id, self._segment.primary_supp_mem.direction)

        for row in range(0, self._mesh_size.min_num_eweb + 1):          # Total number of rows of web element nodes is equal to min_num_ewb + 1
            vertical_spacing_vector = perpendicular_vector * dim_z * row
            long_dim_index = 1
            for column in range(0, column_limit):
                if column > 0:
                    long_spacing_vector += mesh_dim[long_dim_index] * unit_ref_vector
                    long_dim_index += 1
                else:
                    long_spacing_vector = np.zeros(3)

                position_vector = long_spacing_vector + vertical_spacing_vector  # Node position vector in the local coordinate system
                node_coords = position_vector + ref_node1
                print("Node ID:", node_id, node_coords)
                node_id += 1
        return node_id

    def generate_web_elements(self):
        column_limit = self.node_column_limit()                          # Number of nodes along the local longitudinal axis of the segment
        element_id = self._start_element_id
        node_id_array = self._web_node_ref_array

        print("Elementi struka segmenta", self._segment.id, "jakog", self._segment.beam_prop.beam_type, "nosača",
              self._segment.primary_supp_mem.id, self._segment.primary_supp_mem.direction)

        for row in range(0, self._mesh_size.min_num_eweb):
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

    def generate_flange_nodes(self):
        row_limit, column_limit = np.shape(self._flange_node_ref_array)
        mesh_dim = self._edge_plate_nodes                                   # Distances between flange nodes are equal to plate node distances on V1

        ref_node1 = Segment.get_segment_node1(self._segment)                # Reference node 1 coordinates in [mm], origin of the local csy
        ref_node2 = Segment.get_segment_node2(self._segment)                # Reference node 2 coordinates in [mm]
        ref_vector = np.subtract(ref_node2, ref_node1)                      # Reference vector in the direction of the reference segment
        unit_ref_vector = ref_vector / np.linalg.norm(ref_vector)           # Unit reference vector
        direction_unit_vector = self._segment.primary_supp_mem.flange_direction
        tran_spacing_vector = np.zeros(3)
        long_spacing_vector = np.zeros(3)

        flange_element_width = self._mesh_size.get_flange_el_width(self._segment)
        start_node = ref_node1 + direction_unit_vector * flange_element_width * self._mesh_size.num_eaf
        node_id = self._flange_node_ref_array[0, 0]

        print("Čvorovi prirubnice segmenta", self._segment.id, "jakog", self._segment.beam_prop.beam_type, "nosača",
              self._segment.primary_supp_mem.id, self._segment.primary_supp_mem.direction)

        for row in range(0, row_limit):          # Total number of rows of web element nodes is equal to min_num_ewb + 1
            if row > 0:
                tran_spacing_vector += -direction_unit_vector * flange_element_width
            else:
                tran_spacing_vector = np.zeros(3)
            long_dim_index = 1

            for column in range(0, column_limit):
                if column > 0:
                    long_spacing_vector += mesh_dim[long_dim_index] * unit_ref_vector
                    long_dim_index += 1
                else:
                    long_spacing_vector = np.zeros(3)

                position_vector = long_spacing_vector + tran_spacing_vector  # Node position vector in the local coordinate system
                node_coords = position_vector + start_node
                print("Node ID:", node_id, node_coords)
                node_id += 1
        return node_id

    def generate_flange_elements(self, start_element_id):
        # element_id = self.generate_web_elements()       # OVO NEKAKO DRUGAČIJE POSLOŽITI
        element_id = start_element_id
        row_limit, column_limit = np.shape(self._flange_node_ref_array)
        flange_id_array = self._flange_node_ref_array

        print("Elementi prirubnice segmenta", self._segment.id, "jakog", self._segment.beam_prop.beam_type, "nosača",
              self._segment.primary_supp_mem.id, self._segment.primary_supp_mem.direction)

        for row in range(0, row_limit - 1):
            for column in range(0, column_limit - 1):
                node1_id = flange_id_array[row, column]
                node2_id = flange_id_array[row, column + 1]
                node3_id = flange_id_array[row + 1, column + 1]
                node4_id = flange_id_array[row + 1, column]

                # Ovdje instanciraj objekt quad elementa, dodaj čvorove, svojstvo elementu, itd.
                # Quad elementu dodijeliti svojstvo identificirano u self._unique_properties

                print("Quad element ID", element_id, ", Node IDs:", node1_id, node2_id, node3_id, node4_id)
                element_id += 1
        return element_id


class GrillageMesh:
    def __init__(self, mesh_size: MeshSize):
        """
        Class for generating FE mesh on the entire grillage model.

        Contains methods for correcting node overlaps and flange element overlaps.

        :param mesh_size: Calculated mesh dimensions.
        """
        self._mesh_size = mesh_size
        self._mesh_extent = self._mesh_size.mesh_extent

        self._start_node_id = 1
        self._start_element_id = 1

    @property
    def start_node_id(self):
        return self._start_node_id

    @start_node_id.setter
    def start_node_id(self, value):
        self._start_node_id = value

    @property
    def start_element_id(self):
        return self._start_element_id

    @start_element_id.setter
    def start_element_id(self, value):
        self._start_element_id = value

    def generate_plate_mesh(self):
        node_id = self.start_node_id
        element_id = self.start_element_id

        for plate in self._mesh_extent.full_plate_zones.values():
            node_id, element_id = PlatingZoneMesh(self._mesh_size, plate, node_id, element_id).generate_mesh()

        for plate in self._mesh_extent.long_half_plate_zones.values():
            node_id, element_id = PlatingZoneMesh(self._mesh_size, plate, node_id, element_id, AOS.LONGITUDINAL).generate_mesh()

        for plate in self._mesh_extent.tran_half_plate_zones.values():
            node_id, element_id = PlatingZoneMesh(self._mesh_size, plate, node_id, element_id, AOS.TRANSVERSE).generate_mesh()

        for plate in self._mesh_extent.quarter_plate_zone.values():
            node_id, element_id = PlatingZoneMesh(self._mesh_size, plate, node_id, element_id, AOS.BOTH).generate_mesh()

        self.start_node_id = node_id
        self.start_element_id = element_id

    def generate_psm_mesh(self):
        node_id = self.start_node_id
        element_id = self.start_element_id

        if isinstance(self._mesh_size, MeshV1):
            for segment in self._mesh_extent.full_segments.values():
                node_id, element_id = SegmentV1(self._mesh_size, segment, node_id, element_id).generate_mesh()

            for segment in self._mesh_extent.half_segments.values():
                node_id, element_id = SegmentV1(self._mesh_size, segment, node_id, element_id).generate_mesh()

        """
        elif isinstance(self._mesh_size, MeshV2):
            for segment in self._mesh_extent.full_segments.values():
                node_id, element_id = SegmentV2(self._mesh_size, segment, node_id, element_id).generate_mesh()
    
            for segment in self._mesh_extent.half_segments.values():
                node_id, element_id = SegmentV2(self._mesh_size, segment, node_id, element_id).generate_mesh()
        """

    def generate_mesh(self):
        self.generate_plate_mesh()
        self.generate_psm_mesh()

    def check_node_overlap(self):
        # koristiti np.isclose()
        pass
