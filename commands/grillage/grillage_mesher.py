"""
Tool for Grillage Structure Analysis
University of Zagreb, Faculty of Mechanical Engineering and Naval Architecture
Department of Naval Architecture and Ocean Engineering

Master's thesis project

    Gordan Kos, univ.bacc.ing.nav.arch.
    Dr.sc. Pero Prebeg, dipl.ing.


MODULE FOR GRILLAGE FINITE ELEMENT MESH GENERATION

"""
import itertools
from grillage.grillage_model import *
from femdir.custom_exceptions import *
from grillage.grillage_fem import GeoGrillageFEM

# np.set_printoptions(linewidth=400)


# class MeshVariant(Enum):
#     V1 = 1
#     V2 = 2


class ModelCheck:
    def __init__(self, grillage: Grillage):
        """
        Class for checking the grillage model for symmetry
        and mesh generation feasibility.

        Symmetry checks include:
            Relative distances of symmetric primary supporting members.
            Central primary supporting member position.
            Plate property, stiffener layout and stiffener direction.
            Beam property of segments.

        :param grillage: Input grillage model.
        """
        self._grillage = grillage

    def longitudinal_psm_symmetry(self):
        """
        :return: True if relative distances of longitudinal primary supporting
            members are symmetric.
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
        :return: True if relative distances of transverse primary supporting
            members are symmetric.
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
        :return: True if central longitudinal primary supporting member exists
            and has relative distance coordinate 0.5
        """
        n_long = self._grillage.N_longitudinal
        if np.mod(n_long, 2) != 0:
            member_id = int(np.ceil(n_long / 2))
            central_member = self._grillage.longitudinal_members()[member_id]
            if np.isclose(central_member.rel_dist, 0.5):
                return True
            else:
                return False
        else:
            return True

    def central_transversal(self):
        """
        :return: True if central transverse primary supporting member exists
            and has relative distance coordinate 0.5
        """
        n_tran = self._grillage.N_transverse
        if np.mod(n_tran, 2) != 0:
            member_id = int(np.ceil(n_tran / 2))
            central_member = self._grillage.transverse_members()[member_id]
            if np.isclose(central_member.rel_dist, 0.5):
                return True
            else:
                return False
        else:
            return True

    def plate_zone_ID_array(self):
        """
        :return: 2D array of all plating zone IDs arranged to represent
            relative placement of zones on the entire grillage model.
        """
        total_rows = self._grillage.N_longitudinal - 1
        total_columns = self._grillage.N_transverse - 1
        total_zones = total_rows * total_columns

        id_list = np.arange(1, total_zones + 1, 1)
        plating_zone_array = np.reshape(id_list, [total_rows, total_columns])
        return plating_zone_array

    def longitudinal_plate_symmetry(self):
        """
        :return: True if longitudinally symmetric plating zones have the same
            plate property, stiffener layout and stiffener direction.
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
        :return: True if transversely symmetric plating zones have the same
            plate property, stiffener layout and stiffener direction.
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
        :return: True if symmetric longitudinal segments have the same
            beam property.
        """
        same_beam_prop = True

        for member in self._grillage.longitudinal_members().values():
            if member.rel_dist < 0.5:
                for segment_id in range(1, self._grillage.N_transverse):
                    segment = member.segments[segment_id - 1]
                    symm_seg = member.symmetric_member.segments[segment_id - 1]

                    if segment.beam_prop is not symm_seg.beam_prop:
                        same_beam_prop = False
                        break
        return same_beam_prop

    def transverse_segment_symmetry(self):
        """
        :return: True if symmetric transverse segments have the same
            beam property.
        """
        same_beam_prop = True

        for member in self._grillage.transverse_members().values():
            if member.rel_dist < 0.5:
                for segment_id in range(1, self._grillage.N_longitudinal):
                    segment = member.segments[segment_id - 1]
                    symm_seg = member.symmetric_member.segments[segment_id - 1]

                    if segment.beam_prop is not symm_seg.beam_prop:
                        same_beam_prop = False
                        break
        return same_beam_prop

    def long_symmetry_tests(self):
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

    def tran_symmetry_tests(self):
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

        if self.long_symmetry_tests() and self.tran_symmetry_tests():
            return AOS.BOTH
        elif self.long_symmetry_tests():
            return AOS.LONGITUDINAL
        elif self.tran_symmetry_tests():
            return AOS.TRANSVERSE
        else:
            return AOS.NONE

    def mesh_feasibility(self):
        """
        :return: Calls custom exception if input grillage model
            can not be meshed.

        Method stops calculation of limits of grillage mesh generation if
        grillage model does not meet the following criteria:
            1.) Plating zones between two adjacent Primary Supporting Members
                may not have the same stiffener orientation and different
                stiffener spacing.
            2.) Insert another impossibility here
        """
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

                if plate1.stiff_dir == BeamDirection.LONGITUDINAL and \
                        plate2.stiff_dir == BeamDirection.LONGITUDINAL and \
                        spacing1 != spacing2:
                    raise FeasibilityTestFailLong(psm_1.id, psm_2.id, plate1.id,
                                                  plate2.id, spacing1, spacing2)

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

                if plate1.stiff_dir == BeamDirection.TRANSVERSE and \
                        plate2.stiff_dir == BeamDirection.TRANSVERSE and \
                        spacing1 != spacing2:
                    raise FeasibilityTestFailTran(psm_1.id, psm_2.id, plate1.id,
                                                  plate2.id, spacing1, spacing2)

    @staticmethod
    def mesh_V2_feasibility(plate_edge_nodes, flange_edge_nodes):
        p1_end1 = plate_edge_nodes[1]
        p2_end1 = plate_edge_nodes[2]
        p_end1_total = p1_end1 + p2_end1

        p1_end2 = plate_edge_nodes[len(plate_edge_nodes)]
        p2_end2 = plate_edge_nodes[len(plate_edge_nodes) - 1]
        p_end2_total = p1_end2 + p2_end2

        f_end1 = flange_edge_nodes[1]
        f_end2 = flange_edge_nodes[len(flange_edge_nodes)]

        if f_end1 > p_end1_total:
            raise MeshV2FeasibilityFail(p_end1_total, f_end1)

        if f_end2 > p_end2_total:
            raise MeshV2FeasibilityFail(p_end2_total, f_end2)


class MeshExtent:
    def __init__(self, grillage: Grillage, axis_of_symm_override):
        """
        Class for calculating FE mesh extents for the selected grillage model
        and Axis of Symmetry. Contains dictionaries of all plating zones and
        segments to be fully or partially meshed. Generates GeoFEM materials,
        plate and beam properties based on identified mesh extents.

        :param grillage: Input grillage model.
        :param axis_of_symm_override: overrides automatic Axis of Symmetry discovery.

        Class contains the following data:

        plating_zones_ref_array - 2D array of plating zone IDs to be meshed,
                                  based on Axis of Symmetry selection (AOS)
        all_plating_zones - All plating zones included in mesh generation
        full_plate_zones - Plating zones to be fully meshed
        long_half_plate_zones - Plating zones to be split with a longitudinal
                                axis of symmetry
        tran_half_plate_zones - Plating zones to be split with a transverse AOS
        quarter_plate_zone - Plating zone to be split with both AOS
        long_e_split_zone - Plating zones with longitudinal axis of symmetry
                            passing between stiffeners
        tran_e_split_zone - Plating zones with transverse axis of symmetry
                            passing between stiffeners

        all_segments - All segments included in mesh generation.
        full_segments - Segments to be fully meshed
        half_segments - Segments split in half by some axis of symmetry
        """
        self._grillage = grillage
        self.model_check = ModelCheck(self._grillage)
        self._axis_of_symm_override = axis_of_symm_override
        self._aos_input = self.model_check.assign_symmetry()

        self.plating_zones_ref_array = []
        self.all_plating_zones = {}
        self.full_plate_zones = {}
        self.long_half_plate_zones = {}
        self.tran_half_plate_zones = {}
        self.quarter_plate_zone = {}
        self.long_e_split_zone = {}
        self.tran_e_split_zone = {}

        self._all_segment_count = 0
        self.all_segments = {}
        self.full_segments = {}
        self.half_segments = {}

    @property
    def grillage(self):
        return self._grillage

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

    def hc_plate_zone_ref_ID_array(self):
        """
        :return: 2D array of all plating zone IDs arranged to represent
            relative placement of zones on the entire grillage model.

        total_rows - Total number of plating zone rows on the grillage model
        total_columns - Total number of plating zone columns on the model
        """
        total_rows = self._grillage.N_longitudinal - 1
        total_columns = self._grillage.N_transverse - 1
        total_zones = total_rows * total_columns

        id_list = np.arange(1, total_zones + 1, 1)
        plating_zone_array = np.reshape(id_list, [total_rows, total_columns])
        return plating_zone_array

    def longitudinal_psm_extent(self):
        """
        :return: Dictionary of longitudinal primary supporting members to be
            considered for mesh generation, based on input Axis of Symmetry.
        """
        longitudinals = {}

        if self.axis_of_symm is AOS.LONGITUDINAL or \
                self.axis_of_symm is AOS.BOTH:
            n_long = int(np.ceil(self._grillage.N_longitudinal / 2))
            for i in range(1, n_long + 1):
                longitudinals[i] = self._grillage.longitudinal_members()[i]
            return longitudinals
        else:
            return self._grillage.longitudinal_members()

    def transverse_psm_extent(self):
        """
        :return: Dictionary of transverse primary supporting members to be
            considered for mesh generation, based on input Axis of Symmetry.
        """
        transversals = {}

        if self.axis_of_symm is AOS.TRANSVERSE or self.axis_of_symm is AOS.BOTH:
            n_tran = int(np.ceil(self._grillage.N_transverse / 2))
            for i in range(1, n_tran + 1):
                transversals[i] = self._grillage.transverse_members()[i]
            return transversals
        else:
            return self._grillage.transverse_members()

    def identify_long_full_segments(self):
        """
        :return: Identifies Segment objects for full mesh generation; grillage
            with a longitudinal axis of symmetry.
            Stores identified segments in full_segments dictionary.

        n_tran_segments - Number of transverse segmenets to be fully meshed
        """
        n_long = self._grillage.N_longitudinal
        n_tran = self._grillage.N_transverse
        n = 1
        for member in self.longitudinal_psm_extent().values():
            for segment_id in range(0, n_tran - 1):
                self.full_segments[n] = member.segments[segment_id]
                self.add_to_all_segments(member.segments[segment_id])
                n += 1

        n_tran_segments = int(np.floor((n_long - 1) / 2))
        for member in self._grillage.transverse_members().values():
            for segment_id in range(0, n_tran_segments):
                self.full_segments[n] = member.segments[segment_id]
                self.add_to_all_segments(member.segments[segment_id])
                n += 1

    def identify_tran_full_segments(self):
        """
        :return: Identifies Segment objects for full mesh generation; grillage
            with a transverse axis of symmetry.
            Stores identified segments in full_segments dictionary.

        n_long_segments - Number of longitudinal segments to be fully meshed
        """
        n_long = self._grillage.N_longitudinal
        n_tran = self._grillage.N_transverse
        n = 1
        n_long_segments = int(np.floor((n_tran - 1) / 2))
        for member in self._grillage.longitudinal_members().values():
            for segment_id in range(0, n_long_segments):
                self.full_segments[n] = member.segments[segment_id]
                self.add_to_all_segments(member.segments[segment_id])
                n += 1

        for member in self.transverse_psm_extent().values():
            for segment_id in range(0, n_long - 1):
                self.full_segments[n] = member.segments[segment_id]
                self.add_to_all_segments(member.segments[segment_id])
                n += 1

    def identify_both_full_segments(self):
        """
        :return: Identifies Segment objects for full mesh generation; grillage
            with both axis of symmetry.
            Stores identified segments in full_segments dictionary.

        n_long_segments - Number of longitudinal segments to be fully meshed
        n_tran_segments - Number of transverse segmenets to be fully meshed
        """
        n_long = self._grillage.N_longitudinal
        n_tran = self._grillage.N_transverse

        n = 1
        n_long_segments = int(np.floor((n_tran - 1) / 2))
        for member in self.longitudinal_psm_extent().values():
            for segment_id in range(0, n_long_segments):
                self.full_segments[n] = member.segments[segment_id]
                self.add_to_all_segments(member.segments[segment_id])
                n += 1

        n_tran_segments = int(np.floor((n_long - 1) / 2))
        for member in self.transverse_psm_extent().values():
            for segment_id in range(0, n_tran_segments):
                self.full_segments[n] = member.segments[segment_id]
                self.add_to_all_segments(member.segments[segment_id])
                n += 1

    def identify_none_full_segments(self):
        """
        :return: Identifies Segment objects for full mesh generation;
            grillage with no axis of symmetry.
            Stores identified segments in full_segments dictionary.
        """
        n_long = self._grillage.N_longitudinal
        n_tran = self._grillage.N_transverse

        n = 1
        for member in self._grillage.longitudinal_members().values():
            for segment_id in range(0, n_tran - 1):
                self.full_segments[n] = member.segments[segment_id]
                self.add_to_all_segments(member.segments[segment_id])
                n += 1

        for member in self._grillage.transverse_members().values():
            for segment_id in range(0, n_long - 1):
                self.full_segments[n] = member.segments[segment_id]
                self.add_to_all_segments(member.segments[segment_id])
                n += 1

    def identify_long_half_segments(self):
        """
        :return: Identifies Segment objects for half mesh generation;
            grillage with a longitudinal axis of symmetry. Center transverse
            segment exists for even number of longitudinal members.
            Stores identified segments in half_segments dictionary.
        """
        n_long = self._grillage.N_longitudinal
        middle_segment_ID = int(np.ceil((n_long - 1) / 2))
        n = 1

        if np.mod(n_long, 2) == 0:
            for member in self._grillage.transverse_members().values():
                self.half_segments[n] = member.segments[middle_segment_ID - 1]
                self.add_to_all_segments(member.segments[middle_segment_ID - 1])

                n += 1

    def identify_tran_half_segments(self):
        """
        :return: Identifies Segment objects for half mesh generation;
            grillage with a transverse axis of symmetry. Center longitudinal
            segment exists for even number of transverse members.
            Stores identified segments in half_segments dictionary.
        """
        n_tran = self._grillage.N_transverse
        middle_segment_ID = int(np.ceil((n_tran - 1) / 2))
        n = 1

        if np.mod(n_tran, 2) == 0:
            for member in self._grillage.longitudinal_members().values():
                self.half_segments[n] = member.segments[middle_segment_ID - 1]
                self.add_to_all_segments(member.segments[middle_segment_ID - 1])
                n += 1

    def identify_both_half_segments(self):
        """
        :return: Identifies Segment objects for half mesh generation;
            grillage with both axis of symmetry.
            Stores identified segments in half_segments dictionary.
        """
        n_long = self._grillage.N_longitudinal
        n_tran = self._grillage.N_transverse
        long_segment_ID = int(np.ceil((n_tran - 1) / 2))
        tran_segment_ID = int(np.ceil((n_long - 1) / 2))

        n = 1
        if np.mod(n_tran, 2) == 0:
            for member in self.longitudinal_psm_extent().values():
                self.half_segments[n] = member.segments[long_segment_ID - 1]
                self.add_to_all_segments(member.segments[long_segment_ID - 1])
                n += 1

        if np.mod(n_long, 2) == 0:
            for member in self.transverse_psm_extent().values():
                self.half_segments[n] = member.segments[tran_segment_ID - 1]
                self.add_to_all_segments(member.segments[tran_segment_ID - 1])
                n += 1

    def identify_long_full_plate_zones(self):
        """
        :return: Identifies Plate objects for full mesh generation on the entire
            plating zone; grillage with a longitudinal axis of symm.
            Stores identified zones in full_plate_zones dictionary.
        """
        total_rows = self._grillage.N_longitudinal - 1
        plating_zone_array = self.hc_plate_zone_ref_ID_array()

        n_full_rows = int(np.floor(total_rows / 2))
        zone_arr_split = plating_zone_array[0:n_full_rows, :]
        for i in zone_arr_split:
            for j in i:
                self.full_plate_zones[j] = self._grillage.plating()[j]
                self.all_plating_zones[j] = self._grillage.plating()[j]

    def identify_tran_full_plate_zones(self):
        """
        :return: Identifies Plate objects for full mesh generation on the entire
            plating zone; grillage with a transverse axis of symm.
            Stores identified zones in full_plate_zones dictionary.
        """
        total_columns = self._grillage.N_transverse - 1
        plating_zone_array = self.hc_plate_zone_ref_ID_array()

        n_full_columns = int(np.floor(total_columns / 2))
        zone_arr_split = plating_zone_array[:, 0:n_full_columns]
        for i in zone_arr_split:
            for j in i:
                self.full_plate_zones[j] = self._grillage.plating()[j]
                self.all_plating_zones[j] = self._grillage.plating()[j]

    def identify_both_full_plate_zones(self):
        """
        :return: Identifies Plate objects for full mesh generation on the entire
            plating zone; grillage with both axis of symm.
            Stores identified zones in full_plate_zones dictionary.
        """
        total_rows = self._grillage.N_longitudinal - 1
        total_columns = self._grillage.N_transverse - 1
        plating_zone_array = self.hc_plate_zone_ref_ID_array()

        mid_row_id = int(np.floor(total_rows / 2))
        mid_col_id = int(np.floor(total_columns / 2))
        zone_arr_split = plating_zone_array[0:mid_row_id, 0:mid_col_id]
        for i in zone_arr_split:
            for j in i:
                self.full_plate_zones[j] = self._grillage.plating()[j]
                self.all_plating_zones[j] = self._grillage.plating()[j]

    def identify_long_half_plate_zones(self):
        """
        :return: Identifies Plate objects for half mesh generation; plating
            zones split with a longitudinal axis of symmetry.
            Stores identified zones in long_half_plate_zones dictionary.
        """
        total_rows = self._grillage.N_longitudinal - 1
        plating_zone_array = self.hc_plate_zone_ref_ID_array()

        if np.mod(self._grillage.N_longitudinal, 2) == 0:
            mid_row_id = int(np.floor(total_rows / 2))
            middle_row = plating_zone_array[mid_row_id, :]
            for i in middle_row:
                self.long_half_plate_zones[i] = self._grillage.plating()[i]
                self.all_plating_zones[i] = self._grillage.plating()[i]

    def identify_tran_half_plate_zones(self):
        """
        :return: Identifies Plate objects for half mesh generation; plating
            zones split with a transverse axis of symmetry.
            Stores identified zones in tran_half_plate_zones dictionary.
        """
        total_columns = self._grillage.N_transverse - 1
        plating_zone_array = self.hc_plate_zone_ref_ID_array()

        if np.mod(self._grillage.N_transverse, 2) == 0:
            mid_col_id = int(np.floor(total_columns / 2))
            middle_column = plating_zone_array[:, mid_col_id]
            for i in middle_column:
                self.tran_half_plate_zones[i] = self._grillage.plating()[i]
                self.all_plating_zones[i] = self._grillage.plating()[i]

    def identify_both_half_plate_zones(self):
        """
        :return: Identifies Plate objects for half mesh generation; grillage
            with both axis of symm. Stores identified zones in
            long_half_plate_zones and tran_half_plate_zones dictionary.
        """
        total_rows = self._grillage.N_longitudinal - 1
        total_columns = self._grillage.N_transverse - 1
        plating_zone_array = self.hc_plate_zone_ref_ID_array()
        mid_row_id = int(np.floor(total_rows / 2))
        mid_col_id = int(np.floor(total_columns / 2))

        if np.mod(self._grillage.N_longitudinal, 2) == 0:
            middle_row = plating_zone_array[mid_row_id, 0:mid_col_id]
            for i in middle_row:
                self.long_half_plate_zones[i] = self._grillage.plating()[i]
                self.all_plating_zones[i] = self._grillage.plating()[i]

        if np.mod(self._grillage.N_transverse, 2) == 0:
            middle_column = plating_zone_array[0:mid_row_id, mid_col_id]
            for i in middle_column:
                self.tran_half_plate_zones[i] = self._grillage.plating()[i]
                self.all_plating_zones[i] = self._grillage.plating()[i]

    def identify_quarter_plate_zone(self):
        """
        :return: Identifies Plate object for quarter mesh generation;
            plating zone split with both axis of symmetry.
            Only possible if grillage has both axis of symmetry, with even
            number of longitudinal and transverse primary supporting members.
            There can be only one plating zone split with both axis of symmetry.
            Identified zone is saved in quarter_plate_zone dictionary
            for consistency with other methods.
        """
        if np.mod(self._grillage.N_longitudinal, 2) == 0 and \
                np.mod(self._grillage.N_transverse, 2) == 0:
            total_rows = self._grillage.N_longitudinal - 1
            total_columns = self._grillage.N_transverse - 1
            plating_zone_array = self.hc_plate_zone_ref_ID_array()

            mid_row_id = int(np.floor(total_rows / 2))
            mid_col_id = int(np.floor(total_columns / 2))
            plateid = plating_zone_array[mid_row_id, mid_col_id]

            self.quarter_plate_zone[plateid] = self._grillage.plating()[plateid]
            self.all_plating_zones[plateid] = self._grillage.plating()[plateid]

    def get_plate_dim(self, plate: Plate, plate_dim):
        """
        :param plate: Selected plating zone.
        :param plate_dim: Full plating zone dimension.
        :return: Method returns half of the full plating zone dimension if
            selected plate is split by any axis of symmetry.
        """
        if plate.id in self.full_plate_zones:
            return plate_dim

        elif plate.id in self.long_half_plate_zones or \
                plate.id in self.tran_half_plate_zones or \
                plate.id in self.quarter_plate_zone:
            return plate_dim / 2

    def get_long_plate_dim(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Method returns longitudinal dimension of the selected plating
            zone.Returns half of the original value if the zone is split by
            transverse axis of symmetry.
        """
        if plate.id in self.full_plate_zones or \
                plate.id in self.long_half_plate_zones:
            return plate.plate_longitudinal_dim() * 1000

        elif plate.id in self.tran_half_plate_zones or \
                plate.id in self.quarter_plate_zone:
            return (plate.plate_longitudinal_dim() / 2) * 1000

    def get_tran_plate_dim(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Method returns transverse dimension of the selected plating
            zone.Returns half of the original value if the zone is split by
            longitudinal axis of symmetry.
        """
        if plate.id in self.full_plate_zones or \
                plate.id in self.tran_half_plate_zones:
            return plate.plate_transverse_dim() * 1000

        elif plate.id in self.long_half_plate_zones or \
                plate.id in self.quarter_plate_zone:
            return (plate.plate_transverse_dim() / 2) * 1000

    def aos_between_stiffeners(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: True if some axis of symmetry passes between stiffeners on the
            selected plating zone. False return does not indicate stiffener is
            located on some axis of symmetry! If a new stiffener layout
            definition type would be created with nonsymmetric stiffener
            placement, this method needs to be modified.
        """
        stiff_num = plate.get_stiffener_number()
        if (plate.id in self.long_half_plate_zones or
            plate.id in self.quarter_plate_zone) and \
                np.mod(stiff_num, 2) == 0 and \
                plate.stiff_dir == BeamDirection.LONGITUDINAL:
            return True

        elif (plate.id in self.tran_half_plate_zones or
              plate.id in self.quarter_plate_zone) and \
                np.mod(stiff_num, 2) == 0 and \
                plate.stiff_dir == BeamDirection.TRANSVERSE:
            return True
        else:
            return False

    def aos_on_stiffener(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: True if some axis of symmetry is located on and parallel to a
            stiffener on the selected plating zone. False return does not
            indicate axis of symmetry passes between stiffeners! If a new
            stiffener layout definition type would be created with nonsymmetric
            stiffener placement, this method needs to be modified.
        """
        stiff_num = plate.get_stiffener_number()
        if (plate.id in self.long_half_plate_zones or
            plate.id in self.quarter_plate_zone) and \
                np.mod(stiff_num, 2) != 0 and \
                plate.stiff_dir == BeamDirection.LONGITUDINAL:
            return True

        elif (plate.id in self.tran_half_plate_zones or
              plate.id in self.quarter_plate_zone) and \
                np.mod(stiff_num, 2) != 0 and \
                plate.stiff_dir == BeamDirection.TRANSVERSE:
            return True
        else:
            return False

    def aos_on_segment(self, segment: Segment):
        """
        :param segment: Segment of a primary supporting member.
        :return: True if some axis of symmetry is located on and parallel to a
            segment. Used to identify which segments should be assigned half of
            their original web thickness and half of their flange.
        """
        rel_dist = segment.primary_supp_mem.rel_dist
        direction = segment.primary_supp_mem.direction
        aos = self.axis_of_symm
        if (aos == AOS.LONGITUDINAL or aos == AOS.BOTH) and \
                direction == BeamDirection.LONGITUDINAL and \
                np.isclose(rel_dist, 0.5):
            return True

        elif (aos == AOS.TRANSVERSE or aos == AOS.BOTH) and \
                direction == BeamDirection.TRANSVERSE and \
                np.isclose(rel_dist, 0.5):
            return True
        else:
            return False

    def identify_tran_split_elements(self):
        """
        :return: Identifies Plate objects where transverse axis of symmetry
            passes between stiffeners.
            Stores identified zones in tran_e_split_zone dictionary.
        """
        plating_zone_array = self.hc_plate_zone_ref_ID_array()
        for plate in self.tran_half_plate_zones.values():
            if self.aos_between_stiffeners(plate):
                column_id = np.where(plating_zone_array == plate.id)[1]
                split_element_zones = plating_zone_array[:, column_id]
                for i in split_element_zones:
                    for j in i:
                        self.tran_e_split_zone[j] = self._grillage.plating()[j]

    def identify_long_split_elements(self):
        """
        :return: Identifies Plate objects where longitudinal axis of symmetry
            passes between stiffeners.
            Stores identified zones in long_e_split_zone dictionary.
        """
        plating_zone_array = self.hc_plate_zone_ref_ID_array()
        for plate in self.long_half_plate_zones.values():
            if self.aos_between_stiffeners(plate):
                row_id = np.where(plating_zone_array == plate.id)[0]
                split_element_zones = plating_zone_array[row_id, :]
                for i in split_element_zones:
                    for j in i:
                        self.long_e_split_zone[j] = self._grillage.plating()[j]

    def longitudinal_symm_plate_ref_array(self):
        """
        :return: 2D array of plating zone IDs to be meshed for longitudinal
            Axis of Symmetry, arranged to represent their relative placement
            on the grillage model.
        """
        total_rows = self._grillage.N_longitudinal - 1
        plating_zone_array = self.hc_plate_zone_ref_ID_array()
        mid_row_id = int(np.ceil(total_rows / 2))
        self.plating_zones_ref_array = plating_zone_array[0:mid_row_id, :]

    def transverse_symm_plate_ref_array(self):
        """
        :return: 2D array of plating zone IDs to be meshed for transverse
            Axis of Symmetry, arranged to represent their relative placement
            on the grillage model.
        """
        total_columns = self._grillage.N_transverse - 1
        plating_zone_array = self.hc_plate_zone_ref_ID_array()
        mid_col_id = int(np.ceil(total_columns / 2))
        self.plating_zones_ref_array = plating_zone_array[:, 0:mid_col_id]

    def both_symm_plate_ref_array(self):
        """
        :return: 2D array of plating zone IDs to be meshed for both longitudinal
            and transverse Axis of Symmetry, arranged to represent their
            relative placement on the grillage model.
        """
        total_rows = self._grillage.N_longitudinal - 1
        total_columns = self._grillage.N_transverse - 1
        ref_array = self.hc_plate_zone_ref_ID_array()
        mid_row_id = int(np.ceil(total_rows / 2))
        mid_col_id = int(np.ceil(total_columns / 2))
        self.plating_zones_ref_array = ref_array[0:mid_row_id, 0:mid_col_id]

    def grillage_plate_extent(self):
        """
        :return: Determines limits of plate mesh generation based on input
            Axis of Symmetry value. Calls specific methods for identifying which
            plating zones will be fully or partially meshed. If grillage has no
            axis of symmetry, all plating zones inside grillage.plating() will
            be meshed. Saves a reference ID array of all plating zones to be
            meshed into plating_zones_ref_array.
        """
        if self.axis_of_symm is AOS.LONGITUDINAL:
            self.identify_long_full_plate_zones()
            self.identify_long_half_plate_zones()
            self.longitudinal_symm_plate_ref_array()
            self.identify_long_split_elements()

        elif self.axis_of_symm is AOS.TRANSVERSE:
            self.identify_tran_full_plate_zones()
            self.identify_tran_half_plate_zones()
            self.transverse_symm_plate_ref_array()
            self.identify_tran_split_elements()

        elif self.axis_of_symm is AOS.BOTH:
            self.identify_both_full_plate_zones()
            self.identify_both_half_plate_zones()
            self.identify_quarter_plate_zone()
            self.both_symm_plate_ref_array()
            self.identify_long_split_elements()
            self.identify_tran_split_elements()

        else:
            self.full_plate_zones = self._grillage.plating()
            self.all_plating_zones = self._grillage.plating()
            self.plating_zones_ref_array = self.hc_plate_zone_ref_ID_array()

    def add_to_all_segments(self, segment):
        self._all_segment_count += 1
        self.all_segments[self._all_segment_count] = segment

    def grillage_segment_extent(self):
        """
        :return: Determines limits of segment mesh generation based on selected
            Axis of Symmetry. Calls specific methods for identifying which
            segments will be fully or partially meshed. If grillage has no axis
            of symmetry, all segments on the grillage model will be meshed.
        """
        if self.axis_of_symm is AOS.LONGITUDINAL:
            self.identify_long_full_segments()
            self.identify_long_half_segments()

        elif self.axis_of_symm is AOS.TRANSVERSE:
            self.identify_tran_full_segments()
            self.identify_tran_half_segments()

        elif self.axis_of_symm is AOS.BOTH:
            self.identify_both_full_segments()
            self.identify_both_half_segments()

        else:
            self.identify_none_full_segments()

    def add_unique_plate_prop(self, prop: PlateProperty, fem: GeoGrillageFEM):
        """
        Method adds plating zone plate finite element property.
        Checks for duplicate PlateProperty ID in plate_property_IDs.
        If no duplicate exists, a new GeoFEM plate property is created.
        :param prop: Grillage model PlateProperty object.
        :param fem:
        """
        corr_add = self._grillage.corrosion_addition()[1]
        gfe_prop_id = fem.id_prop_count
        if prop.id not in fem.plate_property_IDs.keys():
            tp = prop.tp_net(corr_add, prop.tp)
            mat_id = prop.plate_mat.id
            mat = fem.getMaterial(mat_id)
            fem.add_plate_property(gfe_prop_id, tp, mat)
            fem.plate_property_IDs[prop.id] = gfe_prop_id

    def add_unique_web_prop(self, prop: TBeamProperty, fem: GeoGrillageFEM):
        """
        Method adds segment web plate finite element property.
        Checks for duplicate BeamProperty ID in web_property_IDs.
        If no duplicate exists, a new GeoFEM plate property is created.
        :param prop: Grillage model TBeamProperty object.
        :param fem:
        """
        corr_add = self._grillage.corrosion_addition()[1]
        gfe_prop_id = fem.id_prop_count
        if prop.id not in fem.web_property_IDs.keys():
            tw = prop.tw_net(corr_add)
            mat_id = prop.mat.id
            mat = fem.getMaterial(mat_id)
            fem.add_plate_property(gfe_prop_id, tw, mat)
            fem.web_property_IDs[prop.id] = gfe_prop_id

    def add_unique_half_web_prop(self, prop: TBeamProperty, fem: GeoGrillageFEM):
        """
        Method adds segment web plate finite element property with half tw.
        Checks for duplicate BeamProperty ID in half_web_property_IDs.
        If no duplicate exists, a new GeoFEM plate property is created.
        :param prop: Grillage model TBeamProperty object.
        :param fem:
        """
        corr_add = self._grillage.corrosion_addition()[1]
        gfe_prop_id = fem.id_prop_count
        if prop.id not in fem.half_web_property_IDs.keys():
            tw = prop.tw_net(corr_add) / 2
            mat_id = prop.mat.id
            mat = fem.getMaterial(mat_id)
            fem.add_plate_property(gfe_prop_id, tw, mat)
            fem.half_web_property_IDs[prop.id] = gfe_prop_id

    def add_unique_flange_prop(self, prop: TBeamProperty, fem: GeoGrillageFEM):
        """
        Method adds segment flange plate finite element property.
        Checks for duplicate BeamProperty ID in flange_property_IDs.
        If no duplicate exists, a new GeoFEM plate property is created.
        :param prop: Grillage model TBeamProperty object.
        :param fem:
        """
        corr_add = self._grillage.corrosion_addition()[1]
        gfe_prop_id = fem.id_prop_count
        if prop.id not in fem.flange_property_IDs.keys():
            tf = prop.tf_net(corr_add)
            mat_id = prop.mat.id
            mat = fem.getMaterial(mat_id)
            fem.add_plate_property(gfe_prop_id, tf, mat)
            fem.flange_property_IDs[prop.id] = gfe_prop_id

    def add_unique_T_beam_prop(self, prop: TBeamProperty, fem: GeoGrillageFEM):
        """
        Method adds stiffener layout beam element property.
        Checks for duplicate BeamProperty ID in stiff_beam_prop_IDs.
        If no duplicate exists, a new GeoFEM Beam property is created.
        :param prop: Grillage model TBeamProperty object.
        :param fem:
        """
        corr_add = self._grillage.corrosion_addition()[1]
        gfe_prop_id = fem.id_prop_count
        if prop.id not in fem.stiff_beam_prop_IDs.keys():
            name_str = "T_" + str(prop.hw) + "x" + str(prop.tw)
            name_str += "/" + str(prop.bf) + "x" + str(prop.tf)
            hw = prop.hw_net(corr_add, 0.0)
            tw = prop.tw_net(corr_add)
            bf = prop.bf_net(corr_add)
            tf = prop.tf_net(corr_add)
            mat_id = prop.mat.id
            mat = fem.getMaterial(mat_id)
            fem.add_T_beam_property(name_str, hw, tw, bf, tf, mat)
            fem.stiff_beam_prop_IDs[prop.id] = gfe_prop_id

    def add_unique_half_T_beam_prop(self, prop: TBeamProperty, fem: GeoGrillageFEM):
        """
        Method adds stiffener layout beam element property.
        Checks for duplicate BeamProperty ID in stiff_beam_prop_IDs.
        If no duplicate exists, a new GeoFEM Beam property is created.
        :param prop: Grillage model TBeamProperty object.
        :param fem:
        """
        corr_add = self._grillage.corrosion_addition()[1]
        gfe_prop_id = fem.id_prop_count
        if prop.id not in fem.half_stiff_beam_prop_IDs.keys():
            name_str = "T_" + str(prop.hw) + "x" + str(prop.tw)
            name_str += "/" + str(prop.bf) + "x" + str(prop.tf)
            hw = prop.hw_net(corr_add, 0.0)
            tw = prop.tw_net(corr_add)
            bf = prop.bf_net(corr_add)
            tf = prop.tf_net(corr_add)
            mat_id = prop.mat.id
            mat = fem.getMaterial(mat_id)
            fem.add_half_T_beam_property(name_str, hw, tw, bf, tf, mat)
            fem.half_stiff_beam_prop_IDs[prop.id] = gfe_prop_id

    def add_unique_L_beam_prop(self, prop: LBeamProperty, fem: GeoGrillageFEM):
        """
        Method adds stiffener layout beam element property.
        Checks for duplicate BeamProperty ID in stiff_beam_prop_IDs.
        If no duplicate exists, a new GeoFEM Beam property is created.
        :param prop: Grillage model LBeamProperty object.
        :param fem:
        """
        corr_add = self._grillage.corrosion_addition()[1]
        gfe_prop_id = fem.id_prop_count
        if prop.id not in fem.stiff_beam_prop_IDs.keys():
            name_str = "L_" + str(prop.hw) + "x" + str(prop.tw)
            name_str += "/" + str(prop.bf) + "x" + str(prop.tf)
            hw = prop.hw_net(corr_add, 0.0)
            tw = prop.tw_net(corr_add)
            bf = prop.bf_net(corr_add)
            tf = prop.tf_net(corr_add)
            mat_id = prop.mat.id
            mat = fem.getMaterial(mat_id)
            fem.add_L_beam_property(name_str, hw, tw, bf, tf, mat)
            fem.stiff_beam_prop_IDs[prop.id] = gfe_prop_id

    def add_unique_half_L_beam_prop(self, prop: LBeamProperty, fem: GeoGrillageFEM):
        """
        Method adds stiffener layout beam element property.
        Checks for duplicate BeamProperty ID in stiff_beam_prop_IDs.
        If no duplicate exists, a new GeoFEM Beam property is created.
        :param prop: Grillage model LBeamProperty object.
        :param fem:
        """
        corr_add = self._grillage.corrosion_addition()[1]
        gfe_prop_id = fem.id_prop_count
        if prop.id not in fem.half_stiff_beam_prop_IDs.keys():
            name_str = "L_" + str(prop.hw) + "x" + str(prop.tw)
            name_str += "/" + str(prop.bf) + "x" + str(prop.tf)
            hw = prop.hw_net(corr_add, 0.0)
            tw = prop.tw_net(corr_add)
            bf = prop.bf_net(corr_add)
            tf = prop.tf_net(corr_add)
            mat_id = prop.mat.id
            mat = fem.getMaterial(mat_id)
            fem.add_half_L_beam_property(name_str, hw, tw, bf, tf, mat)
            fem.half_stiff_beam_prop_IDs[prop.id] = gfe_prop_id

    def add_unique_FB_beam_prop(self, prop: FBBeamProperty, fem: GeoGrillageFEM):
        """
        Method adds stiffener layout beam element property.
        Checks for duplicate BeamProperty ID in stiff_beam_prop_IDs.
        If no duplicate exists, a new GeoFEM Beam property is created.
        :param prop: Grillage model FBBeamProperty object.
        :param fem:
        """
        corr_add = self._grillage.corrosion_addition()[1]
        gfe_prop_id = fem.id_prop_count
        if prop.id not in fem.stiff_beam_prop_IDs.keys():
            name_str = "FB_" + str(prop.hw) + "x" + str(prop.tw)
            hw = prop.hw_net(corr_add, 0.0)
            tw = prop.tw_net(corr_add)
            mat_id = prop.mat.id
            mat = fem.getMaterial(mat_id)
            fem.add_FB_beam_property(name_str, hw, tw, mat)
            fem.stiff_beam_prop_IDs[prop.id] = gfe_prop_id

    def add_unique_half_FB_beam_prop(self, prop: FBBeamProperty, fem: GeoGrillageFEM):
        """
        Method adds stiffener layout beam element property.
        Checks for duplicate BeamProperty ID in stiff_beam_prop_IDs.
        If no duplicate exists, a new GeoFEM Beam property is created.
        :param prop: Grillage model FBBeamProperty object.
        :param fem:
        """
        corr_add = self._grillage.corrosion_addition()[1]
        gfe_prop_id = fem.id_prop_count
        if prop.id not in fem.half_stiff_beam_prop_IDs.keys():
            name_str = "FB_" + str(prop.hw) + "x" + str(prop.tw)
            hw = prop.hw_net(corr_add, 0.0)
            tw = prop.tw_net(corr_add)
            mat_id = prop.mat.id
            mat = fem.getMaterial(mat_id)
            fem.add_half_FB_beam_property(name_str, hw, tw, mat)
            fem.half_stiff_beam_prop_IDs[prop.id] = gfe_prop_id

    def add_unique_Bulb_beam_prop(self, prop: BulbBeamProperty, fem: GeoGrillageFEM):
        """
        Method adds stiffener layout beam element property.
        Checks for duplicate BeamProperty ID in stiff_beam_prop_IDs.
        If no duplicate exists, a new GeoFEM Beam property is created.
        :param prop: Grillage model BulbBeamProperty object.
        :param fem:
        """
        corr_add = self._grillage.corrosion_addition()[1]
        gfe_prop_id = fem.id_prop_count
        if prop.id not in fem.stiff_beam_prop_IDs.keys():
            name_str = "HP_" + str(prop.hw_HP) + "x" + str(prop.tw_HP)
            hw = prop.hw_ekv_net(corr_add)
            tw = prop.tw_ekv_net(corr_add)
            bf = prop.bf_ekv_net(corr_add)
            tf = prop.tf_ekv_net(corr_add)
            mat_id = prop.mat.id
            mat = fem.getMaterial(mat_id)
            fem.add_Bulb_beam_property(name_str, hw, tw, bf, tf, mat)
            fem.stiff_beam_prop_IDs[prop.id] = gfe_prop_id

    def add_unique_half_Bulb_beam_prop(self, prop: BulbBeamProperty, fem: GeoGrillageFEM):
        """
        Method adds stiffener layout beam element property.
        Checks for duplicate BeamProperty ID in stiff_beam_prop_IDs.
        If no duplicate exists, a new GeoFEM Beam property is created.
        :param prop: Grillage model BulbBeamProperty object.
        :param fem:
        """
        corr_add = self._grillage.corrosion_addition()[1]
        gfe_prop_id = fem.id_prop_count
        if prop.id not in fem.half_stiff_beam_prop_IDs.keys():
            name_str = "HP_" + str(prop.hw_HP) + "x" + str(prop.tw_HP)
            hw = prop.hw_ekv_net(corr_add)
            tw = prop.tw_ekv_net(corr_add)
            bf = prop.bf_ekv_net(corr_add)
            tf = prop.tf_ekv_net(corr_add)
            mat_id = prop.mat.id
            mat = fem.getMaterial(mat_id)
            fem.add_half_Bulb_beam_property(name_str, hw, tw, bf, tf, mat)
            fem.half_stiff_beam_prop_IDs[prop.id] = gfe_prop_id

    def add_unique_Hat_beam_prop(self, prop: HatBeamProperty, fem: GeoGrillageFEM):
        """
        Method adds stiffener layout beam element property.
        Checks for duplicate BeamProperty ID in stiff_beam_prop_IDs.
        If no duplicate exists, a new GeoFEM Beam property is created.
        :param prop: Grillage model HatBeamProperty object.
        :param fem:
        """
        corr_add = self._grillage.corrosion_addition()[1]
        gfe_prop_id = fem.id_prop_count
        if prop.id not in fem.stiff_beam_prop_IDs.keys():
            name_str = "Hat_" + str(prop.h) + "x" + str(prop.t)
            name_str += "x" + str(prop.bf) + "x" + str(prop.fi) + "°"
            h = prop.h_net(corr_add)
            t = prop.t_net(corr_add)
            bf = prop.bf_net(corr_add)
            fi = prop.fi
            mat_id = prop.mat.id
            mat = fem.getMaterial(mat_id)
            fem.add_Hat_beam_property(name_str, h, t, bf, fi, mat)
            fem.stiff_beam_prop_IDs[prop.id] = gfe_prop_id

    def add_unique_half_Hat_beam_prop(self, prop: HatBeamProperty,
                                      fem: GeoGrillageFEM):
        """
        Method adds stiffener layout beam element property.
        Checks for duplicate BeamProperty ID in stiff_beam_prop_IDs.
        If no duplicate exists, a new GeoFEM Beam property is created.
        :param prop: Grillage model HatBeamProperty object.
        :param fem:
        """
        corr_add = self._grillage.corrosion_addition()[1]
        gfe_prop_id = fem.id_prop_count
        if prop.id not in fem.half_stiff_beam_prop_IDs.keys():
            name_str = "Hat_" + str(prop.h) + "x" + str(prop.t)
            name_str += "x" + str(prop.bf) + "x" + str(prop.fi) + "°"
            h = prop.h_net(corr_add)
            t = prop.t_net(corr_add)
            bf = prop.bf_net(corr_add)
            fi = prop.fi
            mat_id = prop.mat.id
            mat = fem.getMaterial(mat_id)
            fem.add_half_Hat_beam_property(name_str, h, t, bf, fi, mat)
            fem.half_stiff_beam_prop_IDs[prop.id] = gfe_prop_id

    def generate_FEM_material(self, fem: GeoGrillageFEM):
        """
        :return: Generates GeoFEM materials from grillage model materials.
        """
        for mat in self._grillage.material_props().values():
            fem.add_material(mat)

    def generate_FEM_plate_property(self, fem: GeoGrillageFEM):
        """
        :return: Generates GeoFEM plate property from grillage BeamProperty.
        """
        for plate in self.all_plating_zones.values():
            self.add_unique_plate_prop(plate.plate_prop, fem)

        for segment in self.all_segments.values():
            beam_prop = segment.beam_prop
            beam_type = beam_prop.beam_type

            if self.aos_on_segment(segment):
                self.add_unique_half_web_prop(beam_prop, fem)
            else:
                self.add_unique_web_prop(beam_prop, fem)

            if beam_type == BeamType.T or beam_type == BeamType.L:
                self.add_unique_flange_prop(beam_prop, fem)

    def generate_FEM_beam_property(self, fem: GeoGrillageFEM):
        """
        :return: Generates GeoFEM beam property from grillage BeamProperty.
            Does not recognize stiffeners on Axis Of Symmetry. Stiffener on
            AOS identification is done at the moment of beam element generation.
        """
        for plate in self.all_plating_zones.values():
            stiff_prop = plate.stiff_layout.beam_prop

            if stiff_prop.beam_type is BeamType.T:
                self.add_unique_T_beam_prop(stiff_prop, fem)

            elif stiff_prop.beam_type is BeamType.L:
                self.add_unique_L_beam_prop(stiff_prop, fem)

            elif stiff_prop.beam_type is BeamType.FB:
                self.add_unique_FB_beam_prop(stiff_prop, fem)

            elif stiff_prop.beam_type is BeamType.Bulb:
                self.add_unique_Bulb_beam_prop(stiff_prop, fem)

            elif stiff_prop.beam_type is BeamType.Hat:
                self.add_unique_Hat_beam_prop(stiff_prop, fem)

    def generate_half_FEM_beam_property(self, fem: GeoGrillageFEM):
        for plate in self.all_plating_zones.values():
            stiff_prop = plate.stiff_layout.beam_prop
            if self.aos_on_stiffener(plate):
                if stiff_prop.beam_type is BeamType.T:
                    self.add_unique_half_T_beam_prop(stiff_prop, fem)

                elif stiff_prop.beam_type is BeamType.L:
                    self.add_unique_half_L_beam_prop(stiff_prop, fem)

                elif stiff_prop.beam_type is BeamType.FB:
                    self.add_unique_half_FB_beam_prop(stiff_prop, fem)

                elif stiff_prop.beam_type is BeamType.Bulb:
                    self.add_unique_half_Bulb_beam_prop(stiff_prop, fem)

                elif stiff_prop.beam_type is BeamType.Hat:
                    self.add_unique_half_Hat_beam_prop(stiff_prop, fem)

    def grillage_mesh_extent(self):
        """
        :return: Determines limits of grillage mesh generation for plating zones
            and segments based on Axis of Symmetry value. Method stops mesh
             generation if the model does not pass the feasibility test.
        """
        self.model_check.mesh_feasibility()
        self.grillage_plate_extent()
        self.grillage_segment_extent()

    def identify_boundary_coords(self):
        """
        :return: Identifies coordinates for pinned boundary conditions and
            returns list of end_coords.
        """
        symmetry = self.axis_of_symm
        end_coords = []

        for member in self.longitudinal_psm_extent().values():
            end1_coords, end2_coords = member.end_nodes
            if symmetry is AOS.TRANSVERSE or symmetry is AOS.BOTH:
                end_coords.append(end1_coords * [1, 1, 0])
            else:
                end_coords.append(end1_coords * [1, 1, 0])
                end_coords.append(end2_coords * [1, 1, 0])

        for member in self.transverse_psm_extent().values():
            end1_coords, end2_coords = member.end_nodes
            if symmetry is AOS.LONGITUDINAL or symmetry is AOS.BOTH:
                end_coords.append(end1_coords * [1, 1, 0])
            else:
                end_coords.append(end1_coords * [1, 1, 0])
                end_coords.append(end2_coords * [1, 1, 0])
        return end_coords


class MeshSize:
    def __init__(self, mesh_extent: MeshExtent,
                 min_num_ebs: int,
                 min_num_eweb: int,
                 num_eaf: int,
                 flange_aspect_ratio: float,
                 plate_aspect_ratio: float,
                 des_plate_aspect_ratio: float):
        """
        Class for calculating finite element dimensions on the selected
        grillage model.

        Final result of the following methods are distances between edge nodes
        of all structural elements, along both x and y axis, which will be used
        for all node and element generation.

        Mesh control parameters:
        min_num_ebs - Minimum number of elements between stiffeners; default = 1
        min_num_eweb - Minimum number of elements representing the web of a
            primary supporting member along its height; default = 3
        num_eaf - Number of elements across primary supporting member flange;
            default = 1
        flange_aspect_ratio - Maximum aspect ratio value for primary supporting
            member flange quad elements; default = 8
        plate_aspect_ratio - Maximum aspect ratio value for plating and primary
            supporting member web quad elements; default = 3
        des_plate_aspect_ratio - Desirable plating aspect ratio value less than
            the maximum; default = 2

        Mesh dimensions are saved into:
        mesh_dim_x - Final base mesh x dimensions (dim_x) for each column of
            plating zones on the grillage model
        mesh_dim_y - Final base mesh y dimensions (dim_y) for each row of
            plating zones on the grillage model

        :param mesh_extent: FE mesh extents for the selected grillage model
            and Axis of Symmetry.
        """
        self._mesh_extent = mesh_extent

        self._grillage = self._mesh_extent.grillage
        self._axis_of_symm = self._mesh_extent.axis_of_symm

        self._min_num_ebs = min_num_ebs
        self._min_num_eweb = min_num_eweb
        self._num_eaf = num_eaf
        self._flange_aspect_ratio = flange_aspect_ratio
        self._plate_aspect_ratio = plate_aspect_ratio
        self._des_plate_aspect_ratio = des_plate_aspect_ratio

        self._mesh_dim_x = {}
        self._mesh_dim_y = {}

        self.start_nodes = {}       # Dict of starting node ID for all meshed plating zones

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

    @des_plate_aspect_ratio.setter
    def des_plate_aspect_ratio(self, value):
        self._des_plate_aspect_ratio = value
        if value > self._plate_aspect_ratio:
            raise InvalidDesiredAspectRatio(value, self._plate_aspect_ratio)

    @staticmethod
    def save_node_spacing(dictionary: dict, n_dims, dimension):
        """
        Method used for saving edge node spacing dimensions to a dictionary,
        for any number of dimensions of equal value.

        :param dictionary: Dictionary for all node spacing (element) dimensions
            on a plating zone.
        :param n_dims: Number of dimensions to be stored in the Dictionary.
        :param dimension: Element dimension to be stored.
        """
        if dictionary:  # If the dictionary is not empty
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
        :param length: Length L which should be divided into n equal parts, each
            with length x.
        :param spacing: Value to which length x should be closes to.
        :return: Closest divisor of length L, which results in a length x
            closest to given spacing value. If a divisor does not exist,
            this method returns None.
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

            if not res:
                return None
            else:
                for i in range(0, len(res)):
                    if min_diff > abs((length / res[i]) - spacing):
                        min_diff = abs((length / res[i]) - spacing)
                        min_div_id = i
                n = res[min_div_id]
                return n

    @staticmethod
    def find_largest_divisor(length, max_val):
        """
        :param length: Length L which should be divided into n equal parts,
            each with length x.
        :param max_val: Maximum value which length x can not exceed.
        :return: Largest divisor of length L, which results in a length x
            smaller or equal to the given maximum value. If a divisor does not
            exist, such as when length L is decimal, this method returns None.
        """
        i = 1
        divisors = []
        while i <= length:
            if np.mod(length, i) == 0:
                divisors.append(i)
            i += 1

        curr_x_max = 0  # Current maximum value of x
        n = None

        if not divisors:
            return None
        else:  # If a divisor exists
            for i in divisors:
                curr_x = length / i  # Current dimension x value
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
        :param length: Dimension which is to be equally divided.
        :param dim_limit: Maximum dimension allowed for the element along the
            given length, defined by aspect ratio or other limits.
        :return: Element dimension value that is less than the maximum allowed
            and equally divides given length.
        """
        if dim_limit == 0:
            return length
        else:
            n_elements = np.ceil(length / dim_limit)
            dim = length / n_elements
            return dim

    def element_size_perp_to_stiffeners(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Quad element dimension perpendicular to stiffener direction
            based only on number of elements between stiffeners, in [mm].
        """
        stiff_spacing = plate.get_stiffener_spacing() * 1000
        perp_dim = stiff_spacing / self._min_num_ebs
        return perp_dim

    def element_size_para_to_stiffeners(self, plate: Plate, plate_dim):
        """
        :param plate: Selected plating zone.
        :param plate_dim: Plating zone dimension parallel to stiffener
            direction, takes into account Axis of Symmetry.
        :return: Quad element dimension parlallel to stiffener direction, based
            only on quad element dimension perpendicular to stiffener direction,
            desired and maximum aspect ratio. Method attempts to find a divisor
            of plating zone dimension parallel to stiffeners that results in
            quad element size not greater than the desired element dimension
            parallel to stiffener direction. If no divisors are found or the
            divisor is the plating zone dimension itself,
            method refine_plate_element is used.
        """
        y = self.element_size_perp_to_stiffeners(plate)
        des_x_val = y * self._des_plate_aspect_ratio  # Desired element dim
        x = self.refine_plate_element(plate_dim, des_x_val)
        return x

    def get_web_el_height(self, segment: Segment):
        """
        :param segment: Segment of a primary supporting member.
        :return: Quad element dimension along the height of a primary
            supporting member.
        """
        hw = segment.beam_prop.hw
        dim = hw / self._min_num_eweb
        return dim

    def get_flange_el_width(self, segment: Segment):
        """
        :param segment: Segment of a primary supporting member.
        :return: Flange quad element dimension across the width of the flange
            of the selected segment (perpendicular to the segment direction)
            based on the number of elements across the flange. For longitudinal
            segments this method returns dimension dim_yf and for transverse
            segments returns dimension dim_xf.
        """
        corr_add = self._grillage.corrosion_addition()[1]
        bf_net = TBeamProperty.bf_net(segment.beam_prop, corr_add)

        if segment.beam_prop.beam_type == BeamType.T:
            return bf_net / (self._num_eaf * 2)

        elif segment.beam_prop.beam_type == BeamType.L:
            return bf_net / self._num_eaf

        elif segment.beam_prop.beam_type == BeamType.FB:
            return 0

    def get_end1_max_flange_width(self, segment: Segment):
        beam_type = segment.beam_prop.beam_type
        corr_add = self._grillage.corrosion_addition()[1]
        eaf = self._num_eaf
        if beam_type is BeamType.T:
            self_bf_net = segment.beam_prop.bf_net(corr_add) / (eaf + 1)
        elif beam_type is BeamType.L:
            self_bf_net = segment.beam_prop.bf_net(corr_add) / eaf
        else:
            self_bf_net = 0

        seg_end1 = self._grillage.get_parallel_segments(segment)[0]

        if seg_end1:
            end1_type = seg_end1.beam_prop.beam_type
            if end1_type is BeamType.T:
                end1_bf_net = seg_end1.beam_prop.bf_net(corr_add) / (eaf + 1)
            elif end1_type is BeamType.L:
                end1_bf_net = seg_end1.beam_prop.bf_net(corr_add) / eaf
            else:
                end1_bf_net = 0
        else:
            end1_bf_net = 0

        bf_net = np.amax([self_bf_net, end1_bf_net])
        return bf_net

    def get_end2_max_flange_width(self, segment: Segment):
        beam_type = segment.beam_prop.beam_type
        corr_add = self._grillage.corrosion_addition()[1]
        eaf = self._num_eaf
        if beam_type is BeamType.T:
            self_bf_net = segment.beam_prop.bf_net(corr_add) / (eaf + 1)
        elif beam_type is BeamType.L:
            self_bf_net = segment.beam_prop.bf_net(corr_add) / eaf
        else:
            self_bf_net = 0

        seg_end2 = self._grillage.get_parallel_segments(segment)[1]

        if seg_end2:
            end2_type = seg_end2.beam_prop.beam_type
            if end2_type is BeamType.T:
                end2_bf_net = seg_end2.beam_prop.bf_net(corr_add) / (eaf + 1)
            elif end2_type is BeamType.L:
                end2_bf_net = seg_end2.beam_prop.bf_net(corr_add) / eaf
            else:
                end2_bf_net = 0
        else:
            end2_bf_net = 0

        bf_net = np.amax([self_bf_net, end2_bf_net])
        return bf_net

    def get_flange_el_length(self, segment: Segment):
        """
        :param segment: Segment of a primary supporting member.
        :return: Maximum flange quad element length based on flange element
            width and maximum flange aspect ratio (parallel to the segment
            direction). For longitudinal segments this method returns dimension
            dim_xf, and for transverse segments returns dimension dim_yf.
        """
        if self.get_flange_el_width(segment) != 0:
            return self._flange_aspect_ratio * self.get_flange_el_width(segment)
        else:
            return 0

    def get_min_fl_el_len(self, segment1: Segment, segment2: Segment):
        """
        :param segment1: First selected segment.
        :param segment2: Second selected segment.
        :return: Minimum value of two segment flange element lengths.
            If both segments do not have a flange, method returns 0.
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

    def get_min_fl_el_len_between_psm(self, member1: PrimarySuppMem,
                                      member2: PrimarySuppMem):
        """
        :param member1: First selected Primary supporting member.
        :param member2: Second selected Primary supporting member.
        :return: Method identifies all segments between two adjacent primary
            supporting members and returns the minimum flange element length of
            all identified segments. If all flange element length values are 0,
            the el_length_list will be empty and this method will return 0.
        """
        segments_list = self._grillage.segments_between_psm(member1, member2)
        el_length_list = [self.get_flange_el_length(segment)
                          for segment in segments_list
                          if self.get_flange_el_length(segment) != 0]

        if not el_length_list:
            return 0
        else:
            return np.amin(el_length_list)

    def element_size_plating_zone_perp(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Quad element dimension perpendicular to stiffener direction
            based on two criteria:
                1.) Number of elements between stiffeners
                2.) Maximum element dimension
        """
        y = self.element_size_perp_to_stiffeners(plate)

        if plate.stiff_dir == BeamDirection.LONGITUDINAL:
            max_dim = self.get_min_fl_el_len(plate.trans_seg1, plate.trans_seg2)
        else:
            max_dim = self.get_min_fl_el_len(plate.long_seg1, plate.long_seg2)

        if y > max_dim != 0:
            stiff_spacing = plate.get_stiffener_spacing() * 1000
            y = self.refine_plate_element(stiff_spacing, max_dim)
        return y

    def element_size_plating_zone_para(self, plate: Plate, plate_dim):
        """
        :param plate: Selected plating zone.
        :param plate_dim: Plating zone dimension parallel to stiffener direction.
        :return: Quad element dimension parallel to stiffener direction based
            on two criteria:
                1.) Element dimensiom: method element_size_para_to_stiffeners
                2.) Maximum element dimension: method get_min_flange_el_length
        """
        x = self.element_size_para_to_stiffeners(plate, plate_dim)
        if plate.stiff_dir == BeamDirection.LONGITUDINAL:
            max_dim = self.get_min_fl_el_len(plate.long_seg1, plate.long_seg2)
        else:
            max_dim = self.get_min_fl_el_len(plate.trans_seg1, plate.trans_seg2)

        if x > max_dim != 0:
            x = self.refine_plate_element(plate_dim, max_dim)
        return x

    def element_size_plating_zone(self, plate: Plate, plate_dim):
        """
        Method for local consideration of base mesh dimensions dim_x and dim_y,
            for each plating zone individually.

        :param plate: Selected plating zone.
        :param plate_dim: Plating zone dimension parallel to stiffener
            direction. For MeshV1 this value should be calculated using method
            get_reduced_plate_dim and for others using method
            plate_dim_parallel_to_stiffeners.
        :return: Base mesh dimensions x and y that are below the maximum values
            based on plate aspect ratio and flange dimensions,
            considering each plating zone individually.
        """
        if plate.stiff_dir == BeamDirection.LONGITUDINAL:
            dim_x = self.element_size_plating_zone_para(plate, plate_dim)
            dim_y = self.element_size_plating_zone_perp(plate)
        else:
            dim_x = self.element_size_plating_zone_perp(plate)
            dim_y = self.element_size_plating_zone_para(plate, plate_dim)

        if self.element_aspect_ratio(dim_x, dim_y) > self._plate_aspect_ratio:
            dim_xf = self.get_min_fl_el_len(plate.long_seg1, plate.long_seg2)
            dim_yf = self.get_min_fl_el_len(plate.trans_seg1, plate.trans_seg2)

            dim_x_limit = np.minimum(dim_xf, dim_y * self._plate_aspect_ratio)
            dim_y_limit = np.minimum(dim_yf, dim_x * self._plate_aspect_ratio)

            if plate.stiff_dir == BeamDirection.LONGITUDINAL:
                if dim_y > dim_x:
                    stiff_spacing = plate.get_stiffener_spacing() * 1000
                    dim_y = self.refine_plate_element(stiff_spacing, dim_y_limit)
                else:
                    dim_x = self.refine_plate_element(plate_dim, dim_x_limit)

            elif plate.stiff_dir == BeamDirection.TRANSVERSE:
                if dim_y > dim_x:
                    dim_y = self.refine_plate_element(plate_dim, dim_y_limit)
                else:
                    stiff_spacing = plate.get_stiffener_spacing() * 1000
                    dim_x = self.refine_plate_element(stiff_spacing, dim_x_limit)

        return dim_x, dim_y

    def calc_element_base_size(self):
        """
        :return: Calculates the quad element size based on stiffener spacing
            and maximum allowed aspect ratio for all plating zones.
            Returns dictionaries of x (dim_x) and y (dim_y) base dimensions
            for all plating zones.
        """
        mesh_dim_x = {}  # Dimension x for all plating zones
        mesh_dim_y = {}  # Dimension y for all plating zones

        for plate in self._mesh_extent.all_plating_zones.values():
            parallel_dim = plate.plate_dim_parallel_to_stiffeners() * 1000
            plate_dim = self._mesh_extent.get_plate_dim(plate, parallel_dim)
            dim_x, dim_y = self.element_size_plating_zone(plate, plate_dim)
            mesh_dim_x[plate.id] = dim_x
            mesh_dim_y[plate.id] = dim_y

        return mesh_dim_x, mesh_dim_y

    def assign_base_dim_x(self, mesh_dim_x):
        """
        Method for global consideration of base mesh dimensions dim_x,
            for each column of plating zones on the grillage model.

        :param mesh_dim_x: Input dictionary of quad element x dimensions based
            on stiffener spacing and maximum allowed aspect ratio.
        :return: Assigns dimension x for each column of plating zones between
            transverse primary supporting members, identified using plating
            zones reference array based on Axis of Symmetry input.
            If there are no transverse stiffeners on the column of plating
            zones, a minimum value of all saved dim_x for plating zones
            between transverse primary supporting members will be used.
        """
        ref_array = self._mesh_extent.plating_zones_ref_array
        n_rows, n_columns = np.shape(ref_array)
        final_mesh_dim_x = {}

        for column in range(1, n_columns + 1):
            plating_zone_IDs = ref_array[:, column - 1]
            plating_zones = [self._grillage.plating()[plate_id]
                             for plate_id in plating_zone_IDs]

            x_list = [mesh_dim_x[plate.id] for plate in plating_zones]
            dim_x = np.amin(x_list)
            final_mesh_dim_x[column] = dim_x

            for plate in plating_zones:
                if plate.stiff_dir == BeamDirection.TRANSVERSE:
                    tran1 = self._grillage.transverse_members()[column]
                    tran2 = self._grillage.transverse_members()[column + 1]
                    max_x = self.get_min_fl_el_len_between_psm(tran1, tran2)
                    dim_x = mesh_dim_x[plate.id]
                    if dim_x > max_x != 0:  # If dimension x exceeds the maximum
                        stiff_spacing = plate.get_stiffener_spacing() * 1000
                        dim_x = self.refine_plate_element(stiff_spacing, max_x)
                        final_mesh_dim_x[column] = dim_x
                    else:
                        final_mesh_dim_x[column] = dim_x
                    break  # Stop after finding transverse stiffeners

        return final_mesh_dim_x

    def assign_base_dim_y(self, mesh_dim_y):
        """
        Method for global consideration of base mesh dimensions dim_y,
            for each row of plating zones on the grillage model.

        :param mesh_dim_y: Input dictionary of quad element y dimensions based
            on stiffener spacing and maximum allowed aspect ratio.
        :return: Assigns dimension y for each row of plating zones between
            longitudinal primary supporting members, identified using plating
            zones reference array based on Axis of Symmetry input.
            If there are no longitudinal stiffeners on the column of plating
            zones, a minimum value of all saved dim_y for plating zones between
            longitudinal primary supporting members will be used.
        """
        ref_array = self._mesh_extent.plating_zones_ref_array
        n_rows, n_columns = np.shape(ref_array)
        final_mesh_dim_y = {}

        for row in range(1, n_rows + 1):
            plating_zone_IDs = ref_array[row - 1, :]
            plating_zones = [self._grillage.plating()[plate_id]
                             for plate_id in plating_zone_IDs]

            y_list = [mesh_dim_y[plate.id] for plate in plating_zones]
            dim_y = np.amin(y_list)
            final_mesh_dim_y[row] = dim_y

            for plate in plating_zones:
                if plate.stiff_dir == BeamDirection.LONGITUDINAL:
                    long1 = self._grillage.longitudinal_members()[row]
                    long2 = self._grillage.longitudinal_members()[row + 1]
                    max_y = self.get_min_fl_el_len_between_psm(long1, long2)
                    dim_y = mesh_dim_y[plate.id]
                    if dim_y > max_y != 0:  # If dimension y exceeds the maximum
                        stiff_spacing = plate.get_stiffener_spacing() * 1000
                        dim_y = self.refine_plate_element(stiff_spacing, max_y)
                        final_mesh_dim_y[row] = dim_y
                    else:  # If dim_y does not exceed the maximum max_y
                        final_mesh_dim_y[row] = dim_y
                    break  # Stop after finding longitudinal stiffeners

        return final_mesh_dim_y

    def get_base_dim_x(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Base quad element x dimension for any plating zone. Returns the
            value based on longitudinal segment ID from the list of all x
            dimensions: mesh_dim_x.
        """
        if self._mesh_dim_x:
            segment_id = plate.long_seg1.id
            return self._mesh_dim_x[segment_id]
        else:
            raise BaseMeshDimensionX

    def get_base_dim_y(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Base quad element y dimension for any plating zone. Returns the
            value based on transverse segment ID from the list of all y
            dimensions: mesh_dim_y.
        """
        if self._mesh_dim_y:
            segment_id = plate.trans_seg1.id
            return self._mesh_dim_y[segment_id]
        else:
            raise BaseMeshDimensionY

    def get_flange_element_num(self, segment: Segment):
        """
        :param segment: Selected segment of a primary supporting member.
        :return: Number of elements across the flange.
            Method returns 0 for FB beam type.
        """
        flange_element_width = self.get_flange_el_width(segment)
        if flange_element_width == 0:
            return 0
        else:
            return self._num_eaf

    def transition_dim_x(self, plate: Plate):
        """
        Method for local consideration of transition mesh dimensions dim_tr_x
        for each plating zone individually.

        If stiffener direction is transverse, there are transition elements
        next to transverse segments.

        :param plate: Selected plating zone.

        n_elem - Number of elements with dimension dim_x that fit inside
            the remaining distance
        """
        dim_tr_x = 0

        if plate.stiff_dir == BeamDirection.TRANSVERSE:
            dim_x = self.get_base_dim_x(plate)
            dim_y = self.get_base_dim_y(plate)
            stiff_offset = plate.get_equal_stiffener_offset() * 1000

            if stiff_offset < 0:
                raise NegativeRemainingDistance(plate.id)

            n_elem = np.floor(stiff_offset / dim_x)
            if n_elem == 0:
                dim_tr_x = stiff_offset
            else:
                dim_tr_x = stiff_offset - n_elem * dim_x

            if dim_tr_x != 0:
                ar = self.element_aspect_ratio(dim_tr_x, dim_y)
                if ar > self._plate_aspect_ratio and stiff_offset > dim_x:
                    dim_tr_x += dim_x

        return dim_tr_x, dim_tr_x

    def transition_dim_y(self, plate: Plate):
        """
        Method for local consideration of transition mesh dimensions dim_tr_y
        for each plating zone individually.

        If stiffener direction is longitudinal, there are transition elements
        next to longitudinal segments.

        :param plate: Selected plating zone.

        n_elem - Number of elements with dimension dim_y that fit inside
            the remaining distance
        """
        dim_tr_y = 0

        if plate.stiff_dir == BeamDirection.LONGITUDINAL:
            dim_x = self.get_base_dim_x(plate)
            dim_y = self.get_base_dim_y(plate)
            stiff_offset = plate.get_equal_stiffener_offset() * 1000

            if stiff_offset < 0:
                raise NegativeRemainingDistance(plate.id)

            n_elem = np.floor(stiff_offset / dim_y)
            if n_elem == 0:
                dim_tr_y = stiff_offset
            else:
                dim_tr_y = stiff_offset - n_elem * dim_y

            if dim_tr_y != 0:
                ar = self.element_aspect_ratio(dim_tr_y, dim_x)
                if ar > self._plate_aspect_ratio and stiff_offset > dim_y:
                    dim_tr_y += dim_y

        return dim_tr_y, dim_tr_y

    def assign_transition_dim_x(self):
        """
        Method for global consideration of transition mesh dimension x, for each
            column of plating zones.

        :return: Assigns transition elemenet x dimension for each column of
            plating zones between transverse primary supporting members,
            identified using plating zones reference array based on
            Axis of Symmetry.
        """
        ref_array = self._mesh_extent.plating_zones_ref_array
        n_rows, n_columns = np.shape(ref_array)
        tr_el_dim_x = np.zeros((2, n_columns))

        for column in range(1, n_columns + 1):
            plating_zone_IDs = ref_array[:, column - 1]
            plating_zones = [self._grillage.plating()[plate_id] for
                             plate_id in plating_zone_IDs]
            for plate in plating_zones:
                tr_dim_x1, tr_dim_x2 = self.transition_dim_x(plate)
                tr_el_dim_x[0, column - 1] = tr_dim_x1
                tr_el_dim_x[1, column - 1] = tr_dim_x2
                if plate.stiff_dir == BeamDirection.TRANSVERSE:
                    break

        return tr_el_dim_x

    def assign_transition_dim_y(self):
        """
        Method for global consideration of transition mesh dimension y, for
            each row of plating zones.

        :return: Assigns transition elemenet y dimension for each row of plating
            zones between longitudinal primary supporting members,
            identified using plating zones reference array based on
            Axis of Symmetry.
        """
        ref_array = self._mesh_extent.plating_zones_ref_array
        n_rows, n_columns = np.shape(ref_array)
        tr_el_dim_y = np.zeros((2, n_rows))

        for row in range(1, n_rows + 1):
            plating_zone_IDs = ref_array[row - 1, :]
            plating_zones = [self._grillage.plating()[plate_id] for
                             plate_id in plating_zone_IDs]
            for plate in plating_zones:
                tr_dim_y1, tr_dim_y2 = self.transition_dim_y(plate)
                tr_el_dim_y[0, row - 1] = tr_dim_y1
                tr_el_dim_y[1, row - 1] = tr_dim_y2
                if plate.stiff_dir == BeamDirection.LONGITUDINAL:
                    break

        return tr_el_dim_y

    def get_tr_dim_x(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Transition quad element x dimensions for any plating zone.
            Returns the value based on longitudinal segment ID. First value
            (dim_1) represents the element closest to the first transverse
            segment (trans_seg_1) and the second (dim_2) represents the element
            closest to the second transverse segment (trans_seg_2)
            that define the plating zone.
        """
        segment_id = plate.long_seg1.id
        dim_id = segment_id - 1
        dim = self.assign_transition_dim_x()
        dim_1 = dim[0][dim_id]
        dim_2 = dim[1][dim_id]
        return dim_1, dim_2

    def get_tr_dim_y(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Transition quad element x dimensions for any plating zone.
            Returns the value based on transverse segment ID. First value
            (dim_1) represents the element closest to the first longitudinal
            segment (long_seg_1) and the second (dim_2) represents the element
            closest to the second longitudinal segment (long_seg_2)
            that define the plating zone.
        """
        segment_id = plate.trans_seg1.id
        dim_id = segment_id - 1
        dim = self.assign_transition_dim_y()
        dim_1 = dim[0][dim_id]
        dim_2 = dim[1][dim_id]
        return dim_1, dim_2

    def get_base_element_number(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Number of base size elements along x and y dimension of
            the plating zone.

        Index reference:
            x - number of elements or dimension in the longitudinal direction
            y - number of elements or dimension in the transverse direction
        """
        L = self._mesh_extent.get_long_plate_dim(plate)
        B = self._mesh_extent.get_tran_plate_dim(plate)

        dim_x = self.get_base_dim_x(plate)
        dim_y = self.get_base_dim_y(plate)

        tr_el_dim_x1, tr_el_dim_x2 = self.get_tr_dim_x(plate)
        tr_el_dim_y1, tr_el_dim_y2 = self.get_tr_dim_y(plate)

        if plate.id in self._mesh_extent.long_half_plate_zones:
            tr_el_dim_y2 = 0

        elif plate.id in self._mesh_extent.tran_half_plate_zones:
            tr_el_dim_x2 = 0

        elif plate.id in self._mesh_extent.quarter_plate_zone:
            tr_el_dim_x2 = 0
            tr_el_dim_y2 = 0

        dim_x_range = L - tr_el_dim_x1 - tr_el_dim_x2
        dim_y_range = B - tr_el_dim_y1 - tr_el_dim_y2

        n_dim_x = np.round(dim_x_range / dim_x)  # Number of elements with dim_x
        n_dim_y = np.round(dim_y_range / dim_y)  # Number of elements with dim_y
        return n_dim_x, n_dim_y

    def get_long_split_element_num(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Number of base size elements being split in half by
            longitudinal Axis of Symmetry.
        """
        if plate.id in self._mesh_extent.long_e_split_zone:
            B = self._mesh_extent.get_tran_plate_dim(plate)
            tr_el_dim_y1 = self.get_tr_dim_y(plate)[0]
            dim_y = self.get_base_dim_y(plate)
            n_dim_y = self.get_base_element_number(plate)[1]
            remaining_dist = B - tr_el_dim_y1 - dim_y * n_dim_y
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
        :return: Number of base size elements being split in half by
            transverse Axis of Symmetry.
        """
        if plate.id in self._mesh_extent.tran_e_split_zone:
            L = self._mesh_extent.get_long_plate_dim(plate)
            tr_el_dim_x1 = self.get_tr_dim_x(plate)[0]
            dim_x = self.get_base_dim_x(plate)
            n_dim_x = self.get_base_element_number(plate)[0]
            remaining_dist = L - tr_el_dim_x1 - dim_x * n_dim_x
            if remaining_dist == 0:
                return 0
            elif remaining_dist < dim_x:
                return 1
            else:
                return 0
        else:
            return 0

    def get_long_tr_element_num(self, plate: Plate, segment: Segment):
        """
        :param plate: Selected plating zone.
        :param segment: Selected segment of a primary supporting member.
        :return: Number of transition elements between longitudinal segment
            and base elements.
        """
        tr_element_dim_1, tr_element_dim_2 = self.get_tr_dim_y(plate)
        if segment is plate.long_seg1 and tr_element_dim_1 != 0:
            return 1
        elif segment is plate.long_seg2 and tr_element_dim_2 != 0:
            return 1
        else:
            return 0

    def get_tran_tr_element_num(self, plate: Plate, segment: Segment):
        """
        :param plate: Selected plating zone.
        :param segment: Selected segment of a primary supporting member.
        :return: Number of transition elements between transverse segment
            and base elements.
        """
        tr_element_dim_1, tr_element_dim_2 = self.get_tr_dim_x(plate)
        if segment is plate.trans_seg1 and tr_element_dim_1 != 0:
            return 1
        elif segment is plate.trans_seg2 and tr_element_dim_2 != 0:
            return 1
        else:
            return 0

    def plate_edge_node_spacing_x(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Distance between plate nodes along longitudinal edges
         (along x axis), in order, for the selected plating zone.
        """
        base_dim_x = self.get_base_dim_x(plate)
        tr_element_dim_x1, tr_element_dim_x2 = self.get_tr_dim_x(plate)

        tr_el_num_seg_1 = self.get_tran_tr_element_num(plate, plate.trans_seg1)
        tr_el_num_seg_2 = self.get_tran_tr_element_num(plate, plate.trans_seg2)
        base_mesh_element_num = self.get_base_element_number(plate)[0]
        split_element_number = self.get_tran_split_element_num(plate)

        if plate.id in self._mesh_extent.tran_half_plate_zones or \
                plate.id in self._mesh_extent.quarter_plate_zone:
            tr_el_num_seg_2 = 0

        x_spacing = {}
        self.save_node_spacing(x_spacing, tr_el_num_seg_1, tr_element_dim_x1)
        self.save_node_spacing(x_spacing, base_mesh_element_num, base_dim_x)
        self.save_node_spacing(x_spacing, split_element_number, base_dim_x / 2)
        self.save_node_spacing(x_spacing, tr_el_num_seg_2, tr_element_dim_x2)

        return x_spacing

    def plate_edge_node_spacing_y(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Distance between plate nodes along transverse edges
            (along y axis), in order, for the selected plating zone.
        """
        base_dim_y = self.get_base_dim_y(plate)
        tr_element_dim_y1, tr_element_dim_y2 = self.get_tr_dim_y(plate)

        tr_el_num_seg_1 = self.get_long_tr_element_num(plate, plate.long_seg1)
        tr_el_num_seg_2 = self.get_long_tr_element_num(plate, plate.long_seg2)
        base_mesh_element_num = self.get_base_element_number(plate)[1]
        split_element_number = self.get_long_split_element_num(plate)

        if plate.id in self._mesh_extent.long_half_plate_zones or \
                plate.id in self._mesh_extent.quarter_plate_zone:
            tr_el_num_seg_2 = 0

        y_spacing = {}
        self.save_node_spacing(y_spacing, tr_el_num_seg_1, tr_element_dim_y1)
        self.save_node_spacing(y_spacing, base_mesh_element_num, base_dim_y)
        self.save_node_spacing(y_spacing, split_element_number, base_dim_y / 2)
        self.save_node_spacing(y_spacing, tr_el_num_seg_2, tr_element_dim_y2)

        return y_spacing

    def opposite_flange_width(self, segment: Segment):
        """
        :param segment: Selected segment of a primary supporting member.
        :return: Maximum net flange width of both perpendicular segments
            connected at the intersection of primary supporting members.
            Returns value bf_net at both ends of the selected segment.
            Includes checks for beam type and flange direction at both ends.
        """
        corr_add = self._grillage.corrosion_addition()[1]
        end1, end2 = self._grillage.get_perpendicular_segments(segment)
        segment_at_end1 = end1[0]
        segment_at_end2 = end2[0]

        bf_end1 = [segment.beam_prop.bf_net(corr_add) for segment in end1]
        bf_end2 = [segment.beam_prop.bf_net(corr_add) for segment in end2]
        bf_max1 = np.amax(bf_end1)
        bf_max2 = np.amax(bf_end2)

        t_beam_at_end1 = [segment.beam_prop.beam_type for segment in end1 if
                          segment.beam_prop.beam_type is BeamType.T]

        t_beam_at_end2 = [segment.beam_prop.beam_type for segment in end2 if
                          segment.beam_prop.beam_type is BeamType.T]

        if any(t_beam_at_end1):
            bf_max1 = bf_max1 / 2
        if any(t_beam_at_end2):
            bf_max2 = bf_max2 / 2

        l_beam_at_end1 = [segment.beam_prop.beam_type for segment in end1
                          if segment.beam_prop.beam_type is BeamType.L]
        l_beam_at_end2 = [segment.beam_prop.beam_type for segment in end2
                          if segment.beam_prop.beam_type is BeamType.L]

        flange_dir1 = segment_at_end1.primary_supp_mem.flange_direction
        flange_dir2 = segment_at_end2.primary_supp_mem.flange_direction
        central_segment = self._grillage.central_segment(segment)

        # End 1
        if flange_dir1 is FlangeDirection.INWARD \
                and any(l_beam_at_end1) and not t_beam_at_end1 \
                and segment_at_end1.primary_supp_mem.rel_dist > 0.5\
                and not central_segment:
            bf_max1 = 0

        if flange_dir1 is FlangeDirection.OUTWARD \
                and any(l_beam_at_end1) and not t_beam_at_end1 \
                and segment_at_end1.primary_supp_mem.rel_dist < 0.5 \
                and not central_segment:
            bf_max1 = 0

        # End 2
        if flange_dir2 is FlangeDirection.INWARD \
                and any(l_beam_at_end2) and not t_beam_at_end2 \
                and segment_at_end2.primary_supp_mem.rel_dist < 0.5 \
                and not central_segment:
            bf_max2 = 0

        if flange_dir2 is FlangeDirection.OUTWARD \
                and any(l_beam_at_end2) and not t_beam_at_end2 \
                and segment_at_end2.primary_supp_mem.rel_dist > 0.5 \
                and not central_segment:
            bf_max2 = 0

        return bf_max1, bf_max2

    def opposite_flange_element_num(self, segment: Segment):
        """
        :param segment: Selected segment of a primary supporting member.
        :return: Number of elements across the flange of a perpendicular
            segment that influences the FE mesh of the selected segment.
            Method returns values for both ends of the selected segment.
            Returns 0 for FB beam type.
        """
        bf_max1, bf_max2 = self.opposite_flange_width(segment)
        if bf_max1 == 0 and bf_max2 == 0:
            return 0, 0
        elif bf_max1 == 0 and bf_max2 != 0:
            return 0, self._num_eaf
        elif bf_max1 != 0 and bf_max2 == 0:
            return self._num_eaf, 0
        else:
            return self._num_eaf, self._num_eaf

    def flange_transition_dim(self, segment: Segment):
        """
        :param segment: Selected segment.
        :return: Transition element dimensions for any segment on both ends.
        """
        direction = segment.primary_supp_mem.direction
        fl_el_width_1, fl_el_width_2 = self.opposite_flange_width(segment)
        plate_list = self._grillage.segment_common_plates(segment)
        plate = plate_list[0]

        if direction == BeamDirection.LONGITUDINAL:
            base_dim = self.get_base_dim_x(plate)
            perp_base_dim = self.get_base_dim_y(plate)
            plate_tr_dim_1, plate_tr_dim_2 = self.get_tr_dim_x(plate)
        else:
            base_dim = self.get_base_dim_y(plate)
            perp_base_dim = self.get_base_dim_x(plate)
            plate_tr_dim_1, plate_tr_dim_2 = self.get_tr_dim_y(plate)

        if plate_tr_dim_1 != 0:
            remaining_dist1 = plate_tr_dim_1 - fl_el_width_1
        else:
            remaining_dist1 = base_dim - fl_el_width_1

        if plate_tr_dim_2 != 0:
            remaining_dist2 = plate_tr_dim_2 - fl_el_width_2
        else:
            remaining_dist2 = base_dim - fl_el_width_2

        n_elem1 = np.floor(remaining_dist1 / base_dim)
        n_elem2 = np.floor(remaining_dist2 / base_dim)

        if n_elem1 == 0:
            if plate_tr_dim_1 == 0:
                fl_tr_1 = base_dim - fl_el_width_1
            else:
                fl_tr_1 = plate_tr_dim_1 - fl_el_width_1
        else:
            if plate_tr_dim_1 == 0:
                fl_tr_1 = n_elem1 * base_dim - fl_el_width_1
            else:
                fl_tr_1 = plate_tr_dim_1 - n_elem1 * base_dim - fl_el_width_1

        if n_elem2 == 0:
            if plate_tr_dim_2 == 0:
                fl_tr_2 = base_dim - fl_el_width_2
            else:
                fl_tr_2 = plate_tr_dim_2 - fl_el_width_2
        else:
            if plate_tr_dim_1 == 0:
                fl_tr_2 = n_elem2 * base_dim - fl_el_width_2
            else:
                fl_tr_2 = plate_tr_dim_2 - n_elem2 * base_dim - fl_el_width_2

        if fl_tr_1 != 0:
            ar = self.element_aspect_ratio(fl_tr_1, perp_base_dim)
            if ar > self._flange_aspect_ratio and remaining_dist1 > base_dim:
                fl_tr_1 += base_dim

        if fl_tr_2 != 0:
            ar = self.element_aspect_ratio(fl_tr_2, perp_base_dim)
            if ar > self._flange_aspect_ratio and remaining_dist2 > base_dim:
                fl_tr_2 += base_dim

        return fl_tr_1, fl_tr_2

    def flange_base_element_num(self, segment: Segment):
        """
        :param segment: Selected segment of a primary supporting member.
        :return: Number of base size elements the axial direction of the
            selected segment.

        Index reference:
            x - number of elements or dimension in the longitudinal direction
            y - number of elements or dimension in the transverse direction
        """
        plate_list = self._grillage.segment_common_plates(segment)
        plate = plate_list[0]
        direction = segment.primary_supp_mem.direction
        segment_len = segment.segment_len() * 1000
        fl_tr_1, fl_tr_2 = self.flange_transition_dim(segment)
        fl_el_width_1, fl_el_width_2 = self.opposite_flange_width(segment)

        if direction == BeamDirection.LONGITUDINAL:
            base_dim = self.get_base_dim_x(plate)
            split_num = self.get_tran_split_element_num(plate)
        else:
            base_dim = self.get_base_dim_y(plate)
            split_num = self.get_long_split_element_num(plate)

        if segment in self.mesh_extent.half_segments.values():
            segment_len = segment_len / 2
            fl_tr_2 = 0
            fl_el_width_2 = 0

        base_range = segment_len - fl_el_width_1 - fl_tr_1 - fl_tr_2
        base_range -= fl_el_width_2 + split_num * base_dim * 0.5
        base_dim_num = np.round(base_range / base_dim)
        return base_dim_num

    def get_flange_transition_num(self, segment: Segment):
        """
        :param segment: Selected segment.
        :return: Number of transition elements on the segment flange.
            Returns value for both ends of the selected segment.
        """
        tr_el_dim_1, tr_el_dim_2 = self.flange_transition_dim(segment)

        if tr_el_dim_1 == 0:
            tr_num_edge1 = 0
        else:
            tr_num_edge1 = 1

        if tr_el_dim_2 == 0:
            tr_num_edge2 = 0
        else:
            tr_num_edge2 = 1

        return tr_num_edge1, tr_num_edge2

    def flange_edge_node_spacing(self, segment: Segment):
        pass

    def calc_plate_start_node_ids(self):
        """
        :return: Calculates starting node ID for all meshed plating zone and
            saves first node ID into start_nodes.
        """
        node_id_counter = 1
        for plate in self._mesh_extent.all_plating_zones.values():
            edge_nodes_x = self.plate_edge_node_spacing_x(plate)
            edge_nodes_y = self.plate_edge_node_spacing_y(plate)
            column_limit = len(edge_nodes_x) + 1
            row_limit = len(edge_nodes_y) + 1
            total = row_limit * column_limit
            self.start_nodes[plate.id] = node_id_counter
            node_id_counter += total

    def calculate_mesh_dimensions(self):
        """
        :return: Calculates all element dimensions and saves them to
            dictionaries mesh_dim_x, mesh_dim_y.
            These values need to be calculated only once and will be used
            for all node and element generation.
        """

        if self._mesh_extent.axis_of_symm_override:
            print("Selected Axis of Symmetry override:",
                  self._mesh_extent.axis_of_symm_override.name)
            print("Automatic symmetry discovery would have selected:",
                  self._mesh_extent.aos_input.name)

        else:
            print("Automatically discovered grillage model symmetry:",
                  self._mesh_extent.aos_input.name)

        self._mesh_extent.grillage_mesh_extent()
        base_dim_x, base_dim_y = self.calc_element_base_size()
        self.mesh_dim_x = self.assign_base_dim_x(base_dim_x)
        self.mesh_dim_y = self.assign_base_dim_y(base_dim_y)
        self.calc_plate_start_node_ids()


class ElementSizeV1(MeshSize):
    def __init__(self, mesh_extent: MeshExtent,
                 min_num_ebs: int,
                 min_num_eweb: int,
                 num_eaf: int,
                 flange_aspect_ratio: float,
                 plate_aspect_ratio: float,
                 des_plate_aspect_ratio: float):
        """
        Class for calculating finite element dimensions
        specific to meshing variant V1.

        Mesh variant V1 reflects primary supporting member flange elements
        onto plating and has the following limitations:

            1.) Flange width of a primary supporting member has to be constant.
            2.) All primary supporting members need to have the same web height.
            3.) Flange element overlap has to be in the same plane.
            4.) Grillage plating can not be defined with any camber.
        """
        super().__init__(mesh_extent, min_num_ebs, min_num_eweb, num_eaf,
                         flange_aspect_ratio, plate_aspect_ratio, des_plate_aspect_ratio)

    def get_reduced_plate_dim(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Reduced plate dimensions based on plate stiffener orientation
            and flange width for mesh variant V1.
        """
        n_eaf = self._num_eaf
        if plate.stiff_dir == BeamDirection.LONGITUDINAL:
            L = plate.plate_longitudinal_dim() * 1000
            dim_xf1 = self.get_flange_el_width(plate.trans_seg1) * n_eaf
            dim_xf2 = self.get_flange_el_width(plate.trans_seg2) * n_eaf
            L_red = L - dim_xf1 - dim_xf2
            return L_red

        elif plate.stiff_dir == BeamDirection.TRANSVERSE:
            B = plate.plate_transverse_dim() * 1000
            dim_yf1 = self.get_flange_el_width(plate.long_seg1) * n_eaf
            dim_yf2 = self.get_flange_el_width(plate.long_seg2) * n_eaf
            B_red = B - dim_yf1 - dim_yf2
            return B_red

    def calc_element_base_size(self):
        """
        Difference in this method specific to mesh variant V1 is in plate
        dimension plate_dim used for calculating base mesh dimensions.
        If flange elements are reflected onto plating mesh,
        a reduced plate dimension should be used.

        :return: Calculates the quad element size based on stiffener spacing
            and maximum allowed aspect ratio for all plating zones. Returns
            dictionaries of x and y base dimensions for all plating zones.
        """
        mesh_dim_x = {}  # Dimension x for all plating zones
        mesh_dim_y = {}  # Dimension y for all plating zones

        for plate in self._mesh_extent.all_plating_zones.values():
            reduced_dim = self.get_reduced_plate_dim(plate)
            plate_dim = self._mesh_extent.get_plate_dim(plate, reduced_dim)
            dim_x, dim_y = self.element_size_plating_zone(plate, plate_dim)
            mesh_dim_x[plate.id] = dim_x
            mesh_dim_y[plate.id] = dim_y

        return mesh_dim_x, mesh_dim_y

    def transition_dim_x(self, plate: Plate):
        """
        Method for local consideration of transition mesh dimensions dim_tr_x
        for each plating zone individually. Specific to meshing variant V1

        If stiffener direction is transverse, there are transition elements
        next to transverse segments.

        :param plate: Selected plating zone.

        n_elem - Number of elements with dimension dim_x that fit inside
            the remaining distance
        """
        dim_tr_x1 = 0
        dim_tr_x2 = 0

        if plate.stiff_dir == BeamDirection.TRANSVERSE:
            dim_x = self.get_base_dim_x(plate)
            dim_y = self.get_base_dim_y(plate)
            stiff_offset = plate.get_equal_stiffener_offset() * 1000
            n_eaf = self._num_eaf

            flange_width1 = self.get_flange_el_width(plate.trans_seg1) * n_eaf
            flange_width2 = self.get_flange_el_width(plate.trans_seg2) * n_eaf
            remaining_dist1 = stiff_offset - flange_width1
            remaining_dist2 = stiff_offset - flange_width2

            if remaining_dist1 < 0 or remaining_dist2 < 0:
                raise NegativeRemainingDistance(plate.id)

            n_elem1 = np.floor(remaining_dist1 / dim_x)
            n_elem2 = np.floor(remaining_dist2 / dim_x)

            if n_elem1 == 0:
                dim_tr_x1 = stiff_offset - flange_width1
            else:
                dim_tr_x1 = stiff_offset - n_elem1 * dim_x - flange_width1

            if n_elem2 == 0:
                dim_tr_x2 = stiff_offset - flange_width2

            else:
                dim_tr_x2 = stiff_offset - n_elem2 * dim_x - flange_width2

            if dim_tr_x1 != 0:
                ar = self.element_aspect_ratio(dim_tr_x1, dim_y)
                if ar > self._plate_aspect_ratio and remaining_dist1 > dim_x:
                    dim_tr_x1 += dim_x

            if dim_tr_x2 != 0:
                ar = self.element_aspect_ratio(dim_tr_x2, dim_y)
                if ar > self._plate_aspect_ratio and remaining_dist2 > dim_x:
                    dim_tr_x2 += dim_x

        return dim_tr_x1, dim_tr_x2

    def transition_dim_y(self, plate: Plate):
        """
        Method for local consideration of transition mesh dimensions dim_tr_y
        for each plating zone individually. Specific to meshing variant V1.

        If stiffener direction is longitudinal, there are transition elements
        next to longitudinal segments.

        :param plate: Selected plating zone.

        n_elem - Number of elements with dimension dim_y that fit inside
            the remaining distance
        """
        dim_tr_y1 = 0
        dim_tr_y2 = 0

        if plate.stiff_dir == BeamDirection.LONGITUDINAL:
            dim_x = self.get_base_dim_x(plate)
            dim_y = self.get_base_dim_y(plate)
            stiff_offset = plate.get_equal_stiffener_offset() * 1000
            n_eaf = self._num_eaf

            flange_width1 = self.get_flange_el_width(plate.long_seg1) * n_eaf
            flange_width2 = self.get_flange_el_width(plate.long_seg2) * n_eaf
            remaining_dist1 = stiff_offset - flange_width1
            remaining_dist2 = stiff_offset - flange_width2

            if remaining_dist1 < 0 or remaining_dist2 < 0:
                raise NegativeRemainingDistance(plate.id)

            n_elem1 = np.floor(remaining_dist1 / dim_y)
            n_elem2 = np.floor(remaining_dist2 / dim_y)

            if n_elem1 == 0:
                dim_tr_y1 = stiff_offset - flange_width1
            else:
                dim_tr_y1 = stiff_offset - n_elem1 * dim_y - flange_width1

            if n_elem2 == 0:
                dim_tr_y2 = stiff_offset - flange_width2

            else:
                dim_tr_y2 = stiff_offset - n_elem2 * dim_y - flange_width2

            if dim_tr_y1 != 0:
                ar = self.element_aspect_ratio(dim_tr_y1, dim_x)
                if ar > self._plate_aspect_ratio and remaining_dist1 > dim_y:
                    dim_tr_y1 += dim_y

            if dim_tr_y2 != 0:
                ar = self.element_aspect_ratio(dim_tr_y2, dim_x)
                if ar > self._plate_aspect_ratio and remaining_dist2 > dim_y:
                    dim_tr_y2 += dim_y

        return dim_tr_y1, dim_tr_y2

    def get_base_element_number(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Number of base size elements along x and y dimension of
            the plating zone.

        Index reference:
            x - number of elements or dimension in the longitudinal direction
            y - number of elements or dimension in the transverse direction
        """
        L = self._mesh_extent.get_long_plate_dim(plate)
        B = self._mesh_extent.get_tran_plate_dim(plate)
        n_eaf = self._num_eaf
        dim_x = self.get_base_dim_x(plate)
        dim_y = self.get_base_dim_y(plate)

        tr_el_dim_x1, tr_el_dim_x2 = self.get_tr_dim_x(plate)
        tr_el_dim_y1, tr_el_dim_y2 = self.get_tr_dim_y(plate)

        fl_dim_x1 = self.get_flange_el_width(plate.trans_seg1) * n_eaf
        fl_dim_x2 = self.get_flange_el_width(plate.trans_seg2) * n_eaf
        fl_dim_y1 = self.get_flange_el_width(plate.long_seg1) * n_eaf
        fl_dim_y2 = self.get_flange_el_width(plate.long_seg2) * n_eaf

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

        dim_x_range = L - tr_el_dim_x1 - tr_el_dim_x2 - fl_dim_x1 - fl_dim_x2
        dim_y_range = B - tr_el_dim_y1 - tr_el_dim_y2 - fl_dim_y1 - fl_dim_y2

        n_dim_x = np.round(dim_x_range / dim_x)  # Number of elements with dim_x
        n_dim_y = np.round(dim_y_range / dim_y)  # Number of elements with dim_y
        return n_dim_x, n_dim_y

    def get_long_split_element_num(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Number of base size elements being split in half by
            longitudinal Axis of Symmetry.
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
        :return: Number of base size elements being split in half by
            transverse Axis of Symmetry.
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

    def plate_edge_node_spacing_x(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Distance between plate nodes along longitudinal edges
         (along x axis), in order, for the selected plating zone.
        """
        fl_dim_x1 = self.get_flange_el_width(plate.trans_seg1)
        fl_dim_x2 = self.get_flange_el_width(plate.trans_seg2)
        base_dim_x = self.get_base_dim_x(plate)
        tr_element_dim_x1, tr_element_dim_x2 = self.get_tr_dim_x(plate)

        tr_el_num_seg_1 = self.get_tran_tr_element_num(plate, plate.trans_seg1)
        tr_el_num_seg_2 = self.get_tran_tr_element_num(plate, plate.trans_seg2)
        flange_el_num_tran_seg1 = self.get_flange_element_num(plate.trans_seg1)
        flange_el_num_tran_seg2 = self.get_flange_element_num(plate.trans_seg2)
        base_mesh_element_num = self.get_base_element_number(plate)[0]
        split_element_number = self.get_tran_split_element_num(plate)

        if plate.id in self._mesh_extent.tran_half_plate_zones or \
                plate.id in self._mesh_extent.quarter_plate_zone:
            flange_el_num_tran_seg2 = 0
            tr_el_num_seg_2 = 0

        x_spacing = {}
        self.save_node_spacing(x_spacing, flange_el_num_tran_seg1, fl_dim_x1)
        self.save_node_spacing(x_spacing, tr_el_num_seg_1, tr_element_dim_x1)
        self.save_node_spacing(x_spacing, base_mesh_element_num, base_dim_x)
        self.save_node_spacing(x_spacing, split_element_number, base_dim_x / 2)
        self.save_node_spacing(x_spacing, tr_el_num_seg_2, tr_element_dim_x2)
        self.save_node_spacing(x_spacing, flange_el_num_tran_seg2, fl_dim_x2)

        return x_spacing

    def plate_edge_node_spacing_y(self, plate: Plate):
        """
        :param plate: Selected plating zone.
        :return: Distance between plate nodes along transverse edges
            (along y axis), in order, for the selected plating zone.
        """
        fl_dim_y1 = self.get_flange_el_width(plate.long_seg1)
        fl_dim_y2 = self.get_flange_el_width(plate.long_seg2)
        base_dim_y = self.get_base_dim_y(plate)
        tr_element_dim_y1, tr_element_dim_y2 = self.get_tr_dim_y(plate)

        tr_el_num_seg_1 = self.get_long_tr_element_num(plate, plate.long_seg1)
        tr_el_num_seg_2 = self.get_long_tr_element_num(plate, plate.long_seg2)
        flange_el_num_long_seg1 = self.get_flange_element_num(plate.long_seg1)
        flange_el_num_long_seg2 = self.get_flange_element_num(plate.long_seg2)
        base_mesh_element_num = self.get_base_element_number(plate)[1]
        split_element_number = self.get_long_split_element_num(plate)

        if plate.id in self._mesh_extent.long_half_plate_zones or \
                plate.id in self._mesh_extent.quarter_plate_zone:
            flange_el_num_long_seg2 = 0
            tr_el_num_seg_2 = 0

        y_spacing = {}

        self.save_node_spacing(y_spacing, flange_el_num_long_seg1, fl_dim_y1)
        self.save_node_spacing(y_spacing, tr_el_num_seg_1, tr_element_dim_y1)
        self.save_node_spacing(y_spacing, base_mesh_element_num, base_dim_y)
        self.save_node_spacing(y_spacing, split_element_number, base_dim_y / 2)
        self.save_node_spacing(y_spacing, tr_el_num_seg_2, tr_element_dim_y2)
        self.save_node_spacing(y_spacing, flange_el_num_long_seg2, fl_dim_y2)

        return y_spacing


class ElementSizeV2(MeshSize):
    def __init__(self, mesh_extent: MeshExtent,
                 min_num_ebs: int,
                 min_num_eweb: int,
                 num_eaf: int,
                 flange_aspect_ratio: float,
                 plate_aspect_ratio: float,
                 des_plate_aspect_ratio: float):
        """
        Class for calculating finite element dimensions
        specific to meshing variant V2.

        Mesh variant V2 has a more uniform plating mesh with transition plate
        elements between edges of plating zones and base mesh elements.
        Allows for different number of edge nodes along the top and bottom
        row of segment web nodes. Creates a transition mesh on the segment web
        using both quad and triangle elements.

            1.) All primary supporting members need to have the same web height.
            2.) Flange element overlap has to be in the same plane.
            3.) Grillage plating can not be defined with any camber.
            4.) Value of num_eaf can not be greater than 1. ???????????????
        """
        super().__init__(mesh_extent, min_num_ebs, min_num_eweb, num_eaf,
                         flange_aspect_ratio, plate_aspect_ratio, des_plate_aspect_ratio)

    def flange_edge_node_spacing(self, segment: Segment):
        """
        :param segment: Selected segmenet.
        :return: Distance between segment nodes in the axial direction of the
            selected segment, in order in the direction of global csy axis.
        """
        bf_max1, bf_max2 = self.opposite_flange_width(segment)
        tr_dim_1, tr_dim_2 = self.flange_transition_dim(segment)
        fl_el_num1, fl_el_num2 = self.opposite_flange_element_num(segment)
        tr_num_1, tr_num_2 = self.get_flange_transition_num(segment)
        direction = segment.primary_supp_mem.direction
        plate_list = self._grillage.segment_common_plates(segment)
        plate = plate_list[0]
        base_el_num = self.flange_base_element_num(segment)
        if direction == BeamDirection.LONGITUDINAL:
            base_dim = self.get_base_dim_x(plate)
        else:
            base_dim = self.get_base_dim_y(plate)

        split_element_number = self.get_tran_split_element_num(plate)

        if segment in self.mesh_extent.half_segments.values():
            fl_el_num2 = 0
            tr_num_2 = 0

        spacing = {}
        self.save_node_spacing(spacing, fl_el_num1, bf_max1)
        self.save_node_spacing(spacing, tr_num_1, tr_dim_1)
        self.save_node_spacing(spacing, base_el_num, base_dim)
        self.save_node_spacing(spacing, split_element_number, base_dim / 2)
        self.save_node_spacing(spacing, tr_num_2, tr_dim_2)
        self.save_node_spacing(spacing, fl_el_num2, bf_max2)
        return spacing


class PlateMesh:
    def __init__(self, mesh_size: MeshSize, plate: Plate):
        """
        Class for generating FE mesh on a selected plating zone.

        :param mesh_size: Selected input mesh size for plate zone mesh generation.
        :param plate: Selected plate for generating Node and Element objects.
            of the selected plating zone based on Axis Of Symmetry.
        :return: Determines node coordinates and generates finite element Node
            and Element objects on the selected plating zone. Returns last node
            and element ID, to continue node and element numeration.
        """
        self._mesh_size = mesh_size
        self._plate = plate

        self._edge_nodes_x = self._mesh_size.plate_edge_node_spacing_x(plate)
        self._edge_nodes_y = self._mesh_size.plate_edge_node_spacing_y(plate)

    def get_mesh_limits(self):
        """
        :return: Row and column limit values for generating plate nodes and
            elements, based on which Axis of Symmetry splits the plating zone.

        row_limit - Number of node rows on the entire plating zone
        column_limit - Number of node columns on the entire plating zone
        """
        row_limit = len(self._edge_nodes_y) + 1
        column_limit = len(self._edge_nodes_x) + 1
        return row_limit, column_limit

    def get_plate_element_property(self, fem: GeoGrillageFEM):
        """
        :return: Quad element GeoFEM Plate property ID used for plating.
        """
        plate_prop_id = self._plate.plate_prop.id
        fem_prop_id = fem.plate_property_IDs[plate_prop_id]
        return fem_prop_id

    def get_stiffener_beam_property(self, fem: GeoGrillageFEM):
        """
        :return: Beam element GeoFEM property ID used for stiffeners.
        """
        beam_prop_id = self._plate.stiff_layout.beam_prop.id
        fem_prop_id = fem.stiff_beam_prop_IDs[beam_prop_id]
        return fem_prop_id

    def get_half_stiffener_beam_property(self, fem: GeoGrillageFEM):
        """
        :return: Beam element GeoFEM property ID used for stiffeners located on
            Axis Of Symmetry.
        """
        beam_prop_id = self._plate.stiff_layout.beam_prop.id
        fem_prop_id = fem.half_stiff_beam_prop_IDs[beam_prop_id]
        return fem_prop_id

    @staticmethod
    def reference_node_ID_array(start_id, row_limit, column_limit):
        """
        :param start_id:
        :param row_limit:
        :param column_limit:
        :return: 2D array of node IDs arranged to represent relative placement
            of nodes on the plating zone. Used as a reference for quad element
            generation.
        """
        total = row_limit * column_limit  # Total number of nodes
        id_list = np.arange(start_id, start_id + total, 1)
        node_id_array = np.reshape(id_list, [row_limit, column_limit])
        return node_id_array

    def generate_plate_nodes(self, fem: GeoGrillageFEM):
        """
        :return: Generates nodes on the entire plating zone and returns last
            node ID to continue node numeration on other plating zones.
            Overlapping nodes along the edges of the plating zone are saved
            into a separate dictionary for node overlap identification.

        row - Row of nodes along x axis.
        column - Column of nodes along y axis.
        spacing_vector - Node coordinates in the local coordinate system.
        """
        ref_node1 = Segment.get_segment_node1(self._plate.long_seg1)
        ref_node2 = Segment.get_segment_node2(self._plate.long_seg1)
        ref_vector = np.subtract(ref_node2, ref_node1)
        unit_ref_vector = ref_vector / np.linalg.norm(ref_vector)
        normal_vector = np.array((0, 0, 1))  # Vector normal to the plating
        perp_vector = np.cross(normal_vector, unit_ref_vector)

        spacing_vector_x = np.zeros(3)
        spacing_vector_y = np.zeros(3)
        dim_y_index = 1
        row_limit, column_limit = self.get_mesh_limits()

        for row in range(0, row_limit):
            dim_x_index = 1
            if row > 0:
                dim_y = self._edge_nodes_y[dim_y_index]
                spacing_vector_y += dim_y * perp_vector
                dim_y_index += 1
            else:
                spacing_vector_y = np.zeros(3)

            for column in range(0, column_limit):
                if column > 0:
                    dim_x = self._edge_nodes_x[dim_x_index]
                    spacing_vector_x += dim_x * unit_ref_vector
                    dim_x_index += 1
                else:
                    spacing_vector_x = np.zeros(3)

                spacing_vector = spacing_vector_x + spacing_vector_y
                node_coords = spacing_vector + ref_node1

                node = fem.add_node(node_coords)
                if row == 0 or row == row_limit - 1:
                    fem.add_node_to_node_overlaps(node)
                if column == 0 or column == column_limit - 1:
                    fem.add_node_to_node_overlaps(node)

    def generate_plate_elements(self, fem: GeoGrillageFEM):
        """
        :return: Generates elements on the entire plating zone and returns last
            element ID to continue element numeration on other plating zones.

        row_limit - Row of elements along x axis.
        column_limit - Column of elements along y axis.
        """
        row_limit, column_limit = self.get_mesh_limits()
        start_id = self._mesh_size.start_nodes[self._plate.id]
        node_id_array = self.reference_node_ID_array(start_id, row_limit, column_limit)
        fem_prop_id = self.get_plate_element_property(fem)

        id_el_nodes = [None] * 4
        for row in range(0, row_limit - 1):
            for column in range(0, column_limit - 1):
                id_el_nodes[0] = node_id_array[row, column]
                id_el_nodes[1] = node_id_array[row, column + 1]
                id_el_nodes[2] = node_id_array[row + 1, column + 1]
                id_el_nodes[3] = node_id_array[row + 1, column]

                elem = fem.add_quad_element(fem_prop_id, id_el_nodes)
                fem.add_to_plate_elements(elem)

    def identify_beam_nodes(self):
        """
        :return: Method identifies rows or columns of nodes, depending on
            stiffener orientation, where ordinary stiffeners are located in
            the reference node ID array. Returns a list of row or column
            indexes in the reference node ID array for beam element generation.
        """
        stiff_num = self._plate.get_stiffener_number()
        stiff_spacing = self._plate.get_stiffener_spacing() * 1000
        stiff_offset = self._plate.get_equal_stiffener_offset() * 1000
        id_list = []
        dist = 0
        stiff_counter = 0
        spacing = stiff_spacing
        if self._plate.stiff_dir == BeamDirection.TRANSVERSE:
            edge_nodes = self._edge_nodes_x.items()
        else:
            edge_nodes = self._edge_nodes_y.items()

        for item in edge_nodes:
            key, val = item
            dist += val
            if np.isclose(dist, stiff_offset):
                id_list.append(key)
                stiff_counter += 1
            if np.isclose(dist, stiff_offset + spacing):
                if stiff_counter == stiff_num:
                    break
                spacing += stiff_spacing
                id_list.append(key)
                stiff_counter += 1

        return id_list

    def generate_beam_elements(self, fem: GeoGrillageFEM):
        """
        :return: Generates beam elements of ordinary stiffeners between nodes
            identified using method identify_beam_nodes. Returns last beam
            element ID to continue beam element numeration on other plating
            zones.
        """
        row_limit, column_limit = self.get_mesh_limits()
        start_id = self._mesh_size.start_nodes[self._plate.id]
        node_id_array = self.reference_node_ID_array(start_id, row_limit, column_limit)
        stiff_dir = self._plate.stiff_dir
        prop_id = self.get_stiffener_beam_property(fem)
        aos_on_stiff = self._mesh_size.mesh_extent.aos_on_stiffener(self._plate)
        node_id_index_list = self.identify_beam_nodes()
        id_el_nodes = [None] * 2
        dir_vector = np.array([0, 0, -1])

        if stiff_dir == BeamDirection.LONGITUDINAL:
            for index in range(0, len(node_id_index_list)):
                if index == len(node_id_index_list) - 1 and aos_on_stiff:
                    prop_id = self.get_half_stiffener_beam_property(fem)

                for i in range(0, column_limit - 1):
                    id_el_nodes[0] = node_id_array[node_id_index_list[index], i]
                    id_el_nodes[1] = node_id_array[node_id_index_list[index], i + 1]
                    fem.add_beam_element(prop_id, id_el_nodes, dir_vector)
        else:
            for index in range(0, len(node_id_index_list)):
                if index == len(node_id_index_list) - 1 and aos_on_stiff:
                    prop_id = self.get_half_stiffener_beam_property(fem)

                for i in range(0, row_limit - 1):
                    id_el_nodes[0] = node_id_array[i, node_id_index_list[index]]
                    id_el_nodes[1] = node_id_array[i + 1, node_id_index_list[index]]
                    fem.add_beam_element(prop_id, id_el_nodes, dir_vector)

    def generate_mesh(self, fem: GeoGrillageFEM):
        """
        :return: Generates all nodes and elements on the selected plating zone.
            Returns last node and element IDs to continue numeration.
        """
        self.generate_plate_nodes(fem)
        self.generate_plate_elements(fem)
        self.generate_beam_elements(fem)


class SegmentMesh:
    def __init__(self, mesh_size: MeshSize, segment: Segment,
                 start_n_id, start_e_id):
        """
        Class for generating FE mesh on a selected segment.

        :param mesh_size: Selected input mesh size for segment mesh generation.
        :param segment: Selected segment for generating Noodes and Elements.
        :param start_n_id: Starting node ID which allows continued numeration.
        :param start_e_id: Starting element ID which allows continued numeration.
        :return: Determines node coordinates and generates finite element Node
            and Element objects on the selected segment. Returns last node and
            element ID, to continue node and element numbering.

        Class contains the following data:

        edge_plate_nodes - Distances between nodes, at the point of connection
            of primary supporting member web with plating, along the length
            of the selected segment.
        edge_flange_nodes - Distances between nodes, at the point of connection
            of primary supporting member web with its flange, along the length
            of the selected segment.
        edge_nodes_z - Distances between nodes, in order along z axis

        web_node_ref_array - Reference array for generating web nodes.
        """
        self._mesh_size = mesh_size
        self._mesh_extent = self._mesh_size.mesh_extent
        self._segment = segment
        self._start_node_id = start_n_id
        self._start_element_id = start_e_id

        self._edge_plate_nodes = {}
        self._edge_flange_nodes = self._mesh_size.flange_edge_node_spacing(segment)
        self._edge_nodes_z = {}

    def get_plate_edge_nodes(self):
        """
        :return: Identifies a plating zone the segment defines and gets
            distances between edge nodes in the appropriate direction, at the
            connection of segment web and plating zone.
            Loads identified dimensions into dictionary edge_plate_nodes,
            which are necessary for web node generation.
            This method is the first step of generating any segment nodes.
        """
        direction = self._segment.primary_supp_mem.direction
        for plate in self._mesh_extent.all_plating_zones.values():
            segment_defines_plate = plate.test_plate_segment(self._segment)
            if segment_defines_plate:
                if direction == BeamDirection.LONGITUDINAL:
                    edge_nodes = self._mesh_size.plate_edge_node_spacing_x(plate)
                else:
                    edge_nodes = self._mesh_size.plate_edge_node_spacing_y(plate)
                self._edge_plate_nodes = edge_nodes
                break

    def get_inward_flange_vector(self):
        """
        :return: Inward flange direction unit vector, based on Primary
            supporting member direction and relative distance.
        """
        direction = self._segment.primary_supp_mem.direction
        rel_dist = self._segment.primary_supp_mem.rel_dist

        if direction == BeamDirection.LONGITUDINAL:
            if rel_dist < 0.5:
                return np.array((0, 1, 0))
            else:
                return np.array((0, -1, 0))

        if direction == BeamDirection.TRANSVERSE:
            if rel_dist < 0.5:
                return np.array((1, 0, 0))
            else:
                return np.array((-1, 0, 0))

    def get_outward_flange_vector(self):
        """
        :return: Outward flange direction unit vector, based on Primary
            supporting member direction and relative distance.
        """
        direction = self._segment.primary_supp_mem.direction
        rel_dist = self._segment.primary_supp_mem.rel_dist

        if direction == BeamDirection.LONGITUDINAL:
            if rel_dist < 0.5:
                return np.array((0, -1, 0))
            else:
                return np.array((0, 1, 0))

        if direction == BeamDirection.TRANSVERSE:
            if rel_dist < 0.5:
                return np.array((-1, 0, 0))
            else:
                return np.array((1, 0, 0))

    def get_web_element_property_id(self, fem: GeoGrillageFEM):
        """
        :return: Quad element plate property ID used for primary supporting
            member segment web elements.
        """
        beam_prop_id = self._segment.beam_prop.id
        if self._mesh_extent.aos_on_segment(self._segment):
            fem_prop_id = fem.half_web_property_IDs[beam_prop_id]
        else:
            fem_prop_id = fem.web_property_IDs[beam_prop_id]
        return fem_prop_id

    def get_flange_element_property_id(self, fem: GeoGrillageFEM):
        """
        :return: Quad element plate property ID used for primary supporting
            member segment flange elements.
        """
        beam_prop_id = self._segment.beam_prop.id
        fem_prop_id = fem.flange_property_IDs[beam_prop_id]
        return fem_prop_id

    def reference_web_node_ID_array(self):
        """
        :return: 2D array of node IDs arranged to represent relative placement
            of nodes on the primary supporting member segment web. Used as a
            reference for quad element generation.
            Version 1 assumes equal number of nodes in each row and quad
            elements with edges parallel to the global coordinate axis.

        row_limit - Number of web nodes along global z axis.
        column_limit - Number of nodes along the local longitudinal segment axis.
        total - Total number of web nodes.
        """
        column_limit = len(self._edge_plate_nodes) + 1
        row_limit = self._mesh_size.min_num_eweb + 1
        total = row_limit * column_limit
        id_list = np.arange(self._start_node_id, self._start_node_id + total, 1)
        web_node_id_array = np.reshape(id_list, [row_limit, column_limit])
        return web_node_id_array

    def plate_node_row_number(self):
        """
        :return: Number of node rows with plate edge spacing along the web.
        """
        eweb = self._mesh_size.min_num_eweb
        if eweb <= 2:
            p_row_limit = 1
        else:
            p_row_limit = eweb - 1
        return p_row_limit

    def flange_node_row_number(self):
        """
        :return: Number of node rows with flange edge spacing along the web.
        """
        eweb = self._mesh_size.min_num_eweb
        if eweb == 1:
            f_row_limit = 1
        else:
            f_row_limit = 2
        return f_row_limit

    def reference_web_node_ID_list(self):
        """
        :return: List of node IDs arranged to represent relative placement
            of nodes on the primary supporting member segment web. Used as a
            reference for quad element generation. For MeshVariantV2.

        row_limit - Number of web nodes along global z axis.
        column_limit - Number of nodes along the local longitudinal segment axis.
        total - Total number of web nodes.
        """
        web_node_id_list = []
        p_column_limit = len(self._edge_plate_nodes) + 1       # Plate col limit
        f_column_limit = len(self._edge_flange_nodes) + 1      # Flange col limit
        p_row_limit = self.plate_node_row_number()             # Plate row limit
        f_row_limit = self.flange_node_row_number()            # Flange row limit

        # print("column_limit_f:", f_column_limit, ", column_limit_p", p_column_limit,
        #       ", plate row_limit:", p_row_limit, ", flange row_limit:", f_row_limit)
        # total = p_column_limit * p_row_limit + f_column_limit * f_row_limit
        # print("Total number of web nodes:", total)

        start_id = self._start_node_id
        for p_row in range(1, p_row_limit + 1):
            id_list = np.arange(start_id, start_id + p_column_limit, 1)
            start_id += p_column_limit
            web_node_id_list.append(id_list)

        for f_row in range(1, f_row_limit + 1):
            id_list = np.arange(start_id, start_id + f_column_limit, 1)
            start_id += f_column_limit
            web_node_id_list.append(id_list)

        return web_node_id_list

    def generate_web_nodes(self, fem: GeoGrillageFEM):
        pass

    def generate_web_elements(self, fem: GeoGrillageFEM):
        pass

    def ref_flange_node_ID_array(self, flange_start_node):
        """
        :return: 2D array of node IDs arranged to represent relative placement
            of nodes on the primary supporting member segment flange of L beam
            type. Used as a reference for quad element generation.
            Method uses common nodes at the connection of web and flange.

        last_web_node_row - List of common node IDs at the connection.
        row_limit - Number of nodes in the direction of flange width.
        column_limit - Number of nodes along the local longitudinal segment axis.
        total - Total number of flange nodes, excluding the middle row.
        """

        web_node_id_list = self.reference_web_node_ID_list()
        last_web_node_row = web_node_id_list[-1]

        column_limit = len(self._edge_flange_nodes) + 1
        row_limit = self._mesh_size.num_eaf + 1
        total = column_limit * (row_limit - 1)
        id_list = np.arange(flange_start_node, flange_start_node + total, 1)
        id_array = np.reshape(id_list, [row_limit - 1, column_limit])
        flange_node_id_array = np.insert(id_array, 0, last_web_node_row, axis=0)
        return flange_node_id_array

    def generate_flange_nodes(self, fem: GeoGrillageFEM,
                              direction: FlangeDirection, start_node_id):
        pass

    def generate_flange_elements(self, fem: GeoGrillageFEM,
                                 start_node_id, start_element_id):
        pass

    def generate_mesh(self, fem: GeoGrillageFEM):
        """
        :return: Generates all nodes and elements on a segment of a primary
            supporting member. Returns last node and element ID to continue
            numeration on other segments.
        """
        nodes = self._start_node_id
        elements = self._start_element_id

        self.get_plate_edge_nodes()
        beam_type = self._segment.beam_prop.beam_type
        flange_dir = self._segment.primary_supp_mem.flange_direction

        if beam_type is BeamType.T:
            web_nodes = self.generate_web_nodes(fem)
            web_elements = self.generate_web_elements(fem)
            flange_nodes = self.generate_flange_nodes(fem, FlangeDirection.INWARD, web_nodes)
            elements = self.generate_flange_elements(fem, web_nodes, web_elements)

            if not self._mesh_extent.aos_on_segment(self._segment):
                nodes = self.generate_flange_nodes(fem, FlangeDirection.OUTWARD, flange_nodes)
                elements = self.generate_flange_elements(fem, flange_nodes, elements)
            else:
                nodes = flange_nodes

        elif beam_type is BeamType.L:
            web_nodes = self.generate_web_nodes(fem)
            elements = self.generate_web_elements(fem)
            nodes = self.generate_flange_nodes(fem, flange_dir, web_nodes)
            elements = self.generate_flange_elements(fem, web_nodes, elements)

        elif beam_type is BeamType.FB:
            nodes = self.generate_web_nodes(fem)
            elements = self.generate_web_elements(fem)

        return nodes, elements


class SegmentMeshV1(SegmentMesh):
    def __init__(self, mesh_size: MeshSize, segment: Segment,
                 start_n_id, start_e_id):
        """
        CLass for segment mesh generation specific to meshing variant V1.

        Distances between flange nodes are equal to plate edge node distances
        on meshing variant V1.
        """
        super().__init__(mesh_size, segment, start_n_id, start_e_id)
        self._mesh_size = mesh_size
        self._segment = segment
        self._start_node_id = start_n_id        # Starting node ID
        self._start_element_id = start_e_id     # Starting element ID

    def generate_web_nodes(self, fem: GeoGrillageFEM):
        """
        :return: Generates nodes on the web of one segment of a primary
            supporting member and returns last node ID to continue node
            numeration on other segments.

        row_limit - Number of web nodes along global z axis.
        column_limit - Number of nodes along the local longitudinal segment axis.
        dim_z - Vertical dimension of every web element.
        ref_node1 - Reference node 1 coordinates in [mm], origin of the local csy.
        ref_vector - Reference vector in the direction of the local csy.
        perpendicular_vector - Vector opposite of global z axis.
        long_spacing_vector - Longitudinal vector in the direction of PSM.
        position_vector - Node position vector in the local coordinate system.
        """
        row_limit = self._mesh_size.min_num_eweb + 1
        column_limit = int(len(self._edge_plate_nodes) + 1)
        dim_z = self._mesh_size.get_web_el_height(self._segment)
        mesh_dim = self._edge_plate_nodes
        eaf = self._mesh_size.num_eaf

        node_id = self._start_node_id
        ref_node1 = Segment.get_segment_node1(self._segment)
        ref_node2 = Segment.get_segment_node2(self._segment)
        ref_vector = np.subtract(ref_node2, ref_node1)
        unit_ref_vector = ref_vector / np.linalg.norm(ref_vector)
        perpendicular_vector = np.array((0, 0, -1))
        long_spacing_vector = np.zeros(3)

        for row in range(0, row_limit):
            vertical_spacing_vector = perpendicular_vector * dim_z * row
            long_dim_index = 1
            for column in range(0, column_limit):
                if column > 0:
                    long_mesh_dim = mesh_dim[long_dim_index]
                    long_spacing_vector += long_mesh_dim * unit_ref_vector
                    long_dim_index += 1
                else:
                    long_spacing_vector = np.zeros(3)

                position_vector = long_spacing_vector + vertical_spacing_vector
                node_coords = position_vector + ref_node1
                node = fem.add_node(node_coords)

                if row == 0:
                    fem.add_node_to_node_overlaps(node)
                if row == row_limit - 1 and column <= eaf:
                    fem.add_node_to_node_overlaps(node)
                if row == row_limit - 1 and column >= column_limit - eaf - 1:
                    fem.add_node_to_node_overlaps(node)
                if column == 0 or column == column_limit - 1:
                    fem.add_node_to_node_overlaps(node)
                node_id += 1

        return node_id

    def ref_flange_node_ID_array(self, flange_start_node):
        """
        :return: 2D array of node IDs arranged to represent relative placement
            of nodes on the primary supporting member segment flange of L beam
            type. Used as a reference for quad element generation.
            Method uses common nodes at the connection of web and flange.

        last_web_node_row - List of common node IDs at the connection.
        row_limit - Number of nodes in the direction of flange width.
        column_limit - Number of nodes along the local longitudinal segment axis.
        total - Total number of flange nodes, excluding the middle row.
        """
        web_node_id_array = self.reference_web_node_ID_array()
        last_web_node_row = web_node_id_array[-1, :]
        column_limit = len(self._edge_plate_nodes) + 1
        n_eaf = self._mesh_size.num_eaf
        row_limit = n_eaf + 1
        total = column_limit * (row_limit - 1)

        id_list = np.arange(flange_start_node, flange_start_node + total, 1)
        id_array = np.reshape(id_list, [row_limit - 1, column_limit])
        flange_node_id_array = np.insert(id_array, n_eaf, last_web_node_row, axis=0)

        return flange_node_id_array

    def generate_flange_nodes(self, fem: GeoGrillageFEM,
                              direction: FlangeDirection, flange_start_node):
        """
        :param fem:
        :param direction: Selected flange direction for node generation.
        :param flange_start_node: Start flange node for continued node
            generation after web nodes.
        :return: Generates nodes on one side of the flange of one segment of a
            primary supporting member. Returns last node ID to continue node
            numeration on other segments.

        row_limit - Number of flange nodes across the width
        column_limit - Number of nodes along the local longitudinal segment axis.
        mesh_dim - Distances between nodes at the connection of web and flange.
        ref_node1 - Reference node 1 coordinates in [mm], origin of the local csy.
        ref_vector - Reference vector in the direction of the local csy.
        perpendicular_vector - Vector opposite of global z axis.
        long_spacing_vector - Longitudinal vector in the direction of PSM.
        position_vector - Node position vector in the local coordinate system.
        """

        if direction is FlangeDirection.INWARD:
            flange_unit_vector = self.get_inward_flange_vector()
        else:
            flange_unit_vector = self.get_outward_flange_vector()
        ref_array = self.ref_flange_node_ID_array(flange_start_node)
        row_limit, column_limit = np.shape(ref_array)
        row_limit -= 1
        eaf = self._mesh_size.num_eaf
        mesh_dim = self._edge_plate_nodes

        ref_node1 = Segment.get_segment_node1(self._segment)
        ref_node2 = Segment.get_segment_node2(self._segment)
        ref_vector = np.subtract(ref_node2, ref_node1)
        unit_ref_vector = ref_vector / np.linalg.norm(ref_vector)

        width_spacing_vector = np.zeros(3)
        long_spacing_vector = np.zeros(3)

        fl_el_width = self._mesh_size.get_flange_el_width(self._segment)
        z_start_offset = ref_node1 * np.array((1, 1, 0))
        width_offset = flange_unit_vector * fl_el_width * eaf
        start_node = z_start_offset + width_offset
        node_id = flange_start_node
        for row in range(0, row_limit):
            if row > 0:
                width_spacing_vector -= flange_unit_vector * fl_el_width
            else:
                width_spacing_vector = np.zeros(3)

            dim_index = 1
            for column in range(0, column_limit):
                if column > 0:
                    long_spacing_vector += mesh_dim[dim_index] * unit_ref_vector
                    dim_index += 1
                else:
                    long_spacing_vector = np.zeros(3)

                position_vector = long_spacing_vector + width_spacing_vector
                node_coords = position_vector + start_node

                node = fem.add_node(node_coords)
                if column >= (column_limit - eaf - 1):
                    fem.add_node_to_node_overlaps(node)
                if column <= eaf:
                    fem.add_node_to_node_overlaps(node)
                node_id += 1

        return node_id

    def generate_web_elements(self, fem: GeoGrillageFEM):
        """
        :return: Generates elements on the entire segment web and returns last
            element ID to continue element numeration on other segments.

        row_limit - Row of elements along x axis.
        column_limit - Column of elements along y axis.
        """
        column_limit = len(self._edge_plate_nodes) + 1
        element_id = fem.id_element_count
        node_id_array = self.reference_web_node_ID_array()
        fem_prop_id = self.get_web_element_property_id(fem)
        id_el_nodes = [None] * 4
        for row in range(0, self._mesh_size.min_num_eweb):
            for column in range(0, column_limit - 1):
                id_el_nodes[0] = node_id_array[row, column]
                id_el_nodes[1] = node_id_array[row, column + 1]
                id_el_nodes[2] = node_id_array[row + 1, column + 1]
                id_el_nodes[3] = node_id_array[row + 1, column]
                fem.add_quad_element(fem_prop_id, id_el_nodes)
                element_id += 1

        return element_id

    def generate_flange_elements(self, fem: GeoGrillageFEM, flange_start_node, start_element_id):
        element_id = start_element_id
        row_limit, column_limit = np.shape(self.ref_flange_node_ID_array(flange_start_node))
        flange_id_array = self.ref_flange_node_ID_array(flange_start_node)
        fem_prop_id = self.get_flange_element_property_id(fem)

        id_el_nodes = [None] * 4
        for row in range(0, row_limit - 1):
            for column in range(0, column_limit - 1):
                id_el_nodes[0] = flange_id_array[row, column]
                id_el_nodes[1] = flange_id_array[row, column + 1]
                id_el_nodes[2] = flange_id_array[row + 1, column + 1]
                id_el_nodes[3] = flange_id_array[row + 1, column]
                elem = fem.add_quad_element(fem_prop_id, id_el_nodes)
                fem.add_element_to_element_overlaps(elem)
                element_id += 1
        return element_id


class SegmentMeshV2(SegmentMesh):
    def __init__(self, mesh_size: MeshSize, segment: Segment,
                 start_n_id, start_e_id):
        """
        CLass for segment mesh generation specific to meshing variant V2.

        Distances between flange and plate edge nodes do not have to be equal.
        Transition mesh on the segment web uses both deformed quad elements and
        triangle elements.
        """
        super().__init__(mesh_size, segment, start_n_id, start_e_id)
        self._mesh_size = mesh_size
        self._segment = segment
        self._start_node_id = start_n_id
        self._start_element_id = start_e_id

    def idenetify_num_of_tris(self, segment: Segment):
        """
        :param segment: Selected segment.
        :return: Number of triangles at the start and end of transition row.

        Number of triangles depends on the beam type of segments in the
        perpendicular direction to the segment being meshed.
        """
        num_end1, num_end2 = self._mesh_size.opposite_flange_element_num(segment)
        n_tri1 = num_end1
        n_tri2 = num_end2

        plate_edge_nodes = self._edge_plate_nodes
        flange_edge_nodes = self._edge_flange_nodes
        model_check = self._mesh_extent.model_check
        model_check.mesh_V2_feasibility(plate_edge_nodes, flange_edge_nodes)

        p_end1 = plate_edge_nodes[1]
        p_end2 = plate_edge_nodes[len(plate_edge_nodes)]
        f_end1 = flange_edge_nodes[1]
        f_end2 = flange_edge_nodes[len(flange_edge_nodes)]

        if np.isclose(p_end1, f_end1) or f_end1 > p_end1:
            n_tri1 = 0

        if np.isclose(p_end2, f_end2) or f_end2 > p_end2:
            n_tri2 = 0

        if segment in self._mesh_size.mesh_extent.half_segments.values():
            n_tri2 = 0

        return n_tri1, n_tri2

    def generate_web_nodes(self, fem: GeoGrillageFEM):
        """
        :return: Generates nodes on the web of one segment of a primary
            supporting member and returns last node ID to continue node
            numeration on other segments.

        row_limit - Number of web nodes along global z axis.
        column_limit - Number of nodes along the local longitudinal segment axis.
        dim_z - Vertical dimension of every web element.
        ref_node1 - Reference node 1 coordinates in [mm], origin of the local csy.
        ref_vector - Reference vector in the direction of the local csy.
        perpendicular_vector - Vector opposite of global z axis.
        long_spacing_vector - Longitudinal vector in the direction of PSM.
        position_vector - Node position vector in the local coordinate system.
        """
        plate_edge_nodes = self._edge_plate_nodes       # Distances between plate nodes
        flange_edge_nodes = self._edge_flange_nodes     # Distances between flange nodes
        p_row_limit = self.plate_node_row_number()      # Plate row limit
        f_row_limit = self.flange_node_row_number()     # Flange row limit
        p_column_limit = int(len(plate_edge_nodes) + 1)   # Plate col limit
        f_column_limit = int(len(flange_edge_nodes) + 1)  # Flange col limit
        dim_z = self._mesh_size.get_web_el_height(self._segment)
        node_id = self._start_node_id
        eaf = self._mesh_size.num_eaf
        ref_node1 = Segment.get_segment_node1(self._segment)
        ref_node2 = Segment.get_segment_node2(self._segment)
        ref_vector = np.subtract(ref_node2, ref_node1)
        unit_ref_vector = ref_vector / np.linalg.norm(ref_vector)
        perpendicular_vector = np.array((0, 0, -1))
        long_spacing_vector = np.zeros(3)

        for row in range(0, p_row_limit):
            vertical_spacing_vector = perpendicular_vector * dim_z * row
            long_dim_index = 1
            for column in range(0, p_column_limit):
                if column > 0:
                    long_mesh_dim = plate_edge_nodes[long_dim_index]
                    long_spacing_vector += long_mesh_dim * unit_ref_vector
                    long_dim_index += 1
                else:
                    long_spacing_vector = np.zeros(3)

                position_vector = long_spacing_vector + vertical_spacing_vector
                node_coords = position_vector + ref_node1

                node = fem.add_node(node_coords)
                if row == 0:
                    fem.add_node_to_node_overlaps(node)
                if column >= (p_column_limit - eaf):
                    fem.add_node_to_node_overlaps(node)
                if column <= eaf - 1:
                    fem.add_node_to_node_overlaps(node)

                node_id += 1

        for row in range(p_row_limit, p_row_limit + f_row_limit):
            vertical_spacing_vector = perpendicular_vector * dim_z * row
            long_dim_index = 1
            for column in range(0, f_column_limit):
                if column > 0:
                    long_mesh_dim = flange_edge_nodes[long_dim_index]
                    long_spacing_vector += long_mesh_dim * unit_ref_vector
                    long_dim_index += 1
                else:
                    long_spacing_vector = np.zeros(3)

                position_vector = long_spacing_vector + vertical_spacing_vector
                node_coords = position_vector + ref_node1

                node = fem.add_node(node_coords)

                if column >= (f_column_limit - eaf - 1):
                    fem.add_node_to_node_overlaps(node)
                if column <= eaf + 1:
                    fem.add_node_to_node_overlaps(node)

                node_id += 1

        return node_id

    def generate_flange_nodes(self, fem: GeoGrillageFEM,
                              direction: FlangeDirection, flange_start_node):
        """
        :param fem:
        :param direction: Selected flange direction for node generation.
        :param flange_start_node: Start flange node for continued node
            generation after web nodes.
        :return: Generates nodes on one side of the flange of one segment of a
            primary supporting member. Returns last node ID to continue node
            numeration on other segments.

        row_limit - Number of flange nodes across the width
        column_limit - Number of nodes along the local longitudinal segment axis.
        mesh_dim - Distances between nodes at the connection of web and flange.
        ref_node1 - Reference node 1 coordinates in [mm], origin of the local csy.
        ref_vector - Reference vector in the direction of the local csy.
        perpendicular_vector - Vector opposite of global z axis.
        long_spacing_vector - Longitudinal vector in the direction of PSM.
        position_vector - Node position vector in the local coordinate system.
        """

        if direction is FlangeDirection.INWARD:
            flange_unit_vector = self.get_inward_flange_vector()
        else:
            flange_unit_vector = self.get_outward_flange_vector()
        ref_array = self.ref_flange_node_ID_array(flange_start_node)

        row_limit, column_limit = np.shape(ref_array)
        row_limit -= 1
        eaf = self._mesh_size.num_eaf
        mesh_dim = self._edge_flange_nodes

        ref_node1 = Segment.get_segment_node1(self._segment)
        ref_node2 = Segment.get_segment_node2(self._segment)
        ref_vector = np.subtract(ref_node2, ref_node1)
        unit_ref_vector = ref_vector / np.linalg.norm(ref_vector)

        long_spacing_vector = np.zeros(3)

        fl_el_dim = self._mesh_size.get_flange_el_width(self._segment)
        fl_el_dim_end1 = self._mesh_size.get_end1_max_flange_width(self._segment)
        fl_el_dim_end2 = self._mesh_size.get_end2_max_flange_width(self._segment)
        z_start_offset = ref_node1 * np.array((1, 1, 0))
        start_node = z_start_offset + flange_unit_vector * fl_el_dim
        half_segments = self._mesh_size.mesh_extent.half_segments.values()
        node_id = flange_start_node

        for row in range(0, row_limit):
            dim_index = 1
            for column in range(0, column_limit):
                if row > 0:
                    width_spacing_vector = flange_unit_vector * fl_el_dim
                else:
                    width_spacing_vector = np.zeros(3)

                if column > 0:
                    long_spacing_vector += mesh_dim[dim_index] * unit_ref_vector
                    dim_index += 1
                else:
                    long_spacing_vector = np.zeros(3)

                if column < 2:
                    width_spacing_vector = flange_unit_vector * \
                                           (fl_el_dim_end1 - fl_el_dim)

                elif column > column_limit - 3 \
                        and self._segment not in half_segments:
                    width_spacing_vector = flange_unit_vector * \
                                           (fl_el_dim_end2 - fl_el_dim)

                position_vector = long_spacing_vector + width_spacing_vector
                node_coords = position_vector + start_node
                node = fem.add_node(node_coords)

                if column >= (column_limit - eaf - 1):
                    fem.add_node_to_node_overlaps(node)
                if column <= eaf + 1:
                    fem.add_node_to_node_overlaps(node)

                node_id += 1

        return node_id

    def top_web_element_row(self, fem: GeoGrillageFEM, start_element_id):
        """
        :param fem:
        :param start_element_id: Starting element ID for the row.
        :return: Generates quad elements on the top row.
        """
        plate_edge_nodes = self._edge_plate_nodes
        ref_node_list = self.reference_web_node_ID_list()
        p_row_limit = self.plate_node_row_number()
        element_id = start_element_id
        fem_prop_id = self.get_web_element_property_id(fem)
        range_end = len(plate_edge_nodes)
        for row in range(0, p_row_limit - 1):
            local_id = 0
            for quad_id in range(local_id, range_end):
                node1 = ref_node_list[row][local_id]
                node2 = ref_node_list[row][local_id + 1]
                node3 = ref_node_list[row + 1][local_id + 1]
                node4 = ref_node_list[row + 1][local_id]
                node_id_list = [node1, node2, node3, node4]
                fem.add_quad_element(fem_prop_id, node_id_list)

                element_id += 1
                local_id += 1
        return element_id

    def add_first_quad_to_row(self, fem: GeoGrillageFEM, row1, row2, l_id, num):
        """
        :param fem: Grillage mesh.
        :param row1: Input node ID list, upper row of nodes.
        :param row2: Input node ID list, lower row of nodes.
        :param l_id: Local ID of the element in the transition row.
        :param num: Number of quad elements to be generated.
        """
        fem_prop_id = self.get_web_element_property_id(fem)
        for quad_id in range(l_id, l_id + num):
            node1 = row1[quad_id]
            node2 = row1[quad_id + 1]
            node3 = row2[quad_id + 1]
            node4 = row2[quad_id]
            node_id_list = [node1, node2, node3, node4]
            fem.add_quad_element(fem_prop_id, node_id_list)

    def add_first_tria_to_row(self, fem: GeoGrillageFEM, row1, row2, local_id):
        fem_prop_id = self.get_web_element_property_id(fem)
        node1 = row1[local_id]
        node2 = row2[local_id + 1]
        node3 = row2[local_id]
        node_id_list = [node1, node2, node3]
        fem.add_tria_element(fem_prop_id, node_id_list)

    def add_second_tria_to_row(self, fem: GeoGrillageFEM, row1, row2, local_id):
        # Triangle at ref node 2 (right) for n_tri1 == 1:
        fem_prop_id = self.get_web_element_property_id(fem)
        node1 = row1[local_id - 1]
        node2 = row2[local_id + 1]
        node3 = row2[local_id]
        node_id_list = [node1, node2, node3]
        fem.add_tria_element(fem_prop_id, node_id_list)

    def add_quad_after_tria(self, fem: GeoGrillageFEM, row1, row2, local_id, el_num):
        # Central (non deformed) quad elements after a first tria element
        fem_prop_id = self.get_web_element_property_id(fem)
        for quad_id in range(local_id, local_id + el_num):
            node1 = row1[quad_id - 1]
            node2 = row1[quad_id]
            node3 = row2[quad_id + 1]
            node4 = row2[quad_id]
            node_id_list = [node1, node2, node3, node4]
            fem.add_quad_element(fem_prop_id, node_id_list)

    def add_last_quad_to_row(self, fem: GeoGrillageFEM, row1, row2, local_id, el_num):
        # Last element: deformed quad for n_tria1 == 1
        fem_prop_id = self.get_web_element_property_id(fem)
        for quad_id in range(local_id, local_id + el_num):
            node1 = row1[local_id - 2]
            node2 = row1[local_id - 1]
            node3 = row2[local_id + 1]
            node4 = row2[local_id]
            node_id_list = [node1, node2, node3, node4]
            fem.add_quad_element(fem_prop_id, node_id_list)

    def add_def_quad_end2(self, fem: GeoGrillageFEM, row1, row2, local_id, el_num):
        # Last element: deformed quad for
        fem_prop_id = self.get_web_element_property_id(fem)
        for quad_id in range(local_id, local_id + el_num):
            node1 = row1[local_id - 1]
            node2 = row1[local_id]
            node3 = row2[local_id + 1]
            node4 = row2[local_id]
            node_id_list = [node1, node2, node3, node4]
            fem.add_quad_element(fem_prop_id, node_id_list)

    def tr_web_element_row(self, fem: GeoGrillageFEM, start_element_id):
        """
        :param fem:
        :param start_element_id: Starting element ID for the row.
        :return: Generates elements on the transition row using quads and tris.

        Method generates elements on the segment transition web row.

        flange_edge_nodes - input distances between flange nodes.
        p_row_limit - number of node rows with plate edge node spacing.
        ref_node_list - list of nodes on the selected segment for each row.
        tr_row1 - input node ID list, upper row of nodes.
        tr_row2 - input node ID list, lower row of nodes.
        n_quad - number of quad elements at the plating.
        n_tri1, n_tri2 - number of triangle elements at edge1 and edge2
        n_elem - total number of elements in transition row.
        n_nd_quad - total number of non deformed quad elements.
        """
        flange_edge_nodes = self._edge_flange_nodes
        p_row_limit = self.plate_node_row_number()
        ref_node_list = self.reference_web_node_ID_list()

        tr_row1 = ref_node_list[p_row_limit - 1]
        tr_row2 = ref_node_list[p_row_limit]

        n_tri1, n_tri2 = self.idenetify_num_of_tris(self._segment)
        n_elem = len(flange_edge_nodes)
        n_nd_quad = n_elem - 2 * n_tri1 - 2 * n_tri2
        end_id = start_element_id + n_elem - 1
        local_id = 0

        if n_tri1 == 1:
            self.add_first_quad_to_row(fem, tr_row1, tr_row2, local_id, 1)
            local_id += 1
            self.add_first_tria_to_row(fem, tr_row1, tr_row2, local_id)
            local_id += 1
            self.add_quad_after_tria(fem, tr_row1, tr_row2, local_id, n_nd_quad)
        else:
            self.add_first_quad_to_row(fem, tr_row1, tr_row2, local_id, n_nd_quad)

        if n_tri2 == 1:
            local_id += n_nd_quad
            if n_tri1 == 1:
                self.add_second_tria_to_row(fem, tr_row1, tr_row2, local_id)
                local_id += 1
                self.add_last_quad_to_row(fem, tr_row1, tr_row2, local_id, 1)
            else:
                self.add_first_tria_to_row(fem, tr_row1, tr_row2, local_id)
                local_id += 1
                self.add_def_quad_end2(fem, tr_row1, tr_row2, local_id, 1)

        return end_id

    def bot_web_element_row(self, fem: GeoGrillageFEM, start_element_id):
        """
        :param fem:
        :param start_element_id: Starting element ID for the row.
        :return: Generates quad elements on the bottom row.
        """
        flange_edge_nodes = self._edge_flange_nodes
        ref_node_list = self.reference_web_node_ID_list()
        p_row_limit = self.plate_node_row_number()          # Plate row limit
        f_row_limit = self.flange_node_row_number()         # Flange row limit
        element_id = start_element_id
        range_end = len(flange_edge_nodes)
        fem_prop_id = self.get_web_element_property_id(fem)

        for row in range(p_row_limit, p_row_limit + f_row_limit - 1):
            local_id = 0
            for quad_id in range(local_id, range_end):
                node1 = ref_node_list[row][local_id]
                node2 = ref_node_list[row][local_id + 1]
                node3 = ref_node_list[row + 1][local_id + 1]
                node4 = ref_node_list[row + 1][local_id]
                node_id_list = [node1, node2, node3, node4]
                fem.add_quad_element(fem_prop_id, node_id_list)
                element_id += 1
                local_id += 1
        return element_id

    def generate_web_elements(self, fem: GeoGrillageFEM):
        """
        Method for generating web elements one row at a time, using a combination
        of quad and triangle elements.

        If there are 3 or more rows of elements, the transition row will always
        be the second row from the flange.
        If there are 2 rows of elements, the transition row is on top.
        If there is 1 row of elements, only the transition row exists.

        :return: Generates elements on the entire segment web and returns last
            element ID to continue element numeration on other segments.
        """
        if self._mesh_size.min_num_eweb == 1:
            end_id = self.tr_web_element_row(fem, self._start_element_id)

        elif self._mesh_size.min_num_eweb == 2:
            end_id = self.tr_web_element_row(fem, self._start_element_id)
            end_id = self.bot_web_element_row(fem, end_id)

        else:
            end_id = self.top_web_element_row(fem, self._start_element_id)
            for row in range(0, self._mesh_size.min_num_eweb - 3):
                end_id = self.top_web_element_row(fem, end_id)
            end_id = self.tr_web_element_row(fem, end_id)
            end_id = self.bot_web_element_row(fem, end_id)

        return end_id

    def generate_flange_elements(self, fem: GeoGrillageFEM, flange_start_node,
                                 start_element_id):
        element_id = start_element_id
        node_id_arr = self.ref_flange_node_ID_array(flange_start_node)
        row_limit, column_limit = np.shape(node_id_arr)
        flange_id_array = self.ref_flange_node_ID_array(flange_start_node)
        fem_prop_id = self.get_flange_element_property_id(fem)

        id_el_nodes = [None] * 4
        for row in range(0, row_limit - 1):
            for column in range(0, column_limit - 1):
                id_el_nodes[0] = flange_id_array[row, column]
                id_el_nodes[1] = flange_id_array[row, column + 1]
                id_el_nodes[2] = flange_id_array[row + 1, column + 1]
                id_el_nodes[3] = flange_id_array[row + 1, column]
                elem = fem.add_quad_element(fem_prop_id, id_el_nodes)
                fem.add_element_to_element_overlaps(elem)
                element_id += 1

        return element_id


class GrillageMesh:
    def __init__(self, grillage: Grillage, axis_of_symm_override: AOS = None):
        """
        Class for generating FE mesh on the entire grillage model.

        :param grillage: Input Grillage model.
        :param axis_of_symm_override: Selected Axis of Symmetry override.
        """
        self._grillage = grillage
        self.mesh_extent = MeshExtent(self._grillage, axis_of_symm_override)

    def generate_FEM_property(self, fem: GeoGrillageFEM):
        self.mesh_extent.generate_FEM_material(fem)
        self.mesh_extent.generate_FEM_plate_property(fem)
        self.mesh_extent.generate_FEM_beam_property(fem)
        self.mesh_extent.generate_half_FEM_beam_property(fem)

    def generate_plate_mesh(self, size: MeshSize, fem: GeoGrillageFEM):
        for plate in self.mesh_extent.full_plate_zones.values():
            pzm = PlateMesh(size, plate)
            pzm.generate_mesh(fem)

        for plate in self.mesh_extent.long_half_plate_zones.values():
            pzm = PlateMesh(size, plate)
            pzm.generate_mesh(fem)

        for plate in self.mesh_extent.tran_half_plate_zones.values():
            pzm = PlateMesh(size, plate)
            pzm.generate_mesh(fem)

        for plate in self.mesh_extent.quarter_plate_zone.values():
            pzm = PlateMesh(size, plate)
            pzm.generate_mesh(fem)

    def generate_psm_mesh(self, size: MeshSize, fem: GeoGrillageFEM):
        pass

    def pinned_bc_node_group(self, fem: GeoGrillageFEM):
        """
        Adds all nodes at the ends of primary supporting members to nodal group
        for pinned boundary conditions. Group ID = 1
        """
        nodes_dict = fem.nodes.values()
        end_nodes = self.mesh_extent.identify_boundary_coords()
        coords = [node.p for node in nodes_dict]
        id_list = [node.id for node in nodes_dict]
        coords_1 = np.expand_dims(coords, 0)
        coords_2 = np.expand_dims(end_nodes, 1)
        boolean_array = np.isclose(coords_1, coords_2, rtol=1e-3).all(-1)

        pinned_boundary_nodes = {}
        index_list = np.where(boolean_array)[1]
        for index in index_list:
            node_id = id_list[index]
            node = fem.nodes[node_id]
            pinned_boundary_nodes[node.id] = node

        fem.add_node_group(1, pinned_boundary_nodes)

    @staticmethod
    def long_symm_bc_node_group(fem: GeoGrillageFEM):
        """
        :return: Dictionary of all nodes on the longitudinal Axis of Symmetry
            for symmetry boundary conditions. Group ID = 2
        """
        nodes_dict = fem.nodes.values()
        y_coords = [node.p[1] for node in nodes_dict]
        id_list = [node.id for node in nodes_dict]

        y_max = np.max(y_coords)
        boolean_array = np.isclose(y_coords, y_max, rtol=1e-2)
        node_index = np.where(boolean_array)
        node_index = np.concatenate(node_index)

        long_symm_nodes = {}
        for index in range(0, len(node_index)):
            node_id = id_list[node_index[index]]
            node = fem.nodes[node_id]
            long_symm_nodes[node.id] = node

        fem.add_node_group(2, long_symm_nodes)

    @staticmethod
    def tran_symm_bc_node_group(fem: GeoGrillageFEM):
        """
        :return: Dictionary of all nodes on the transverse Axis of Symmetry
            for symmetry boundary conditions. Group ID = 3
        """
        nodes_dict = fem.nodes.values()
        x_coords = [node.p[0] for node in nodes_dict]
        id_list = [node.id for node in nodes_dict]

        x_max = np.max(x_coords)
        boolean_array = np.isclose(x_coords, x_max, rtol=1e-2)
        node_index = np.where(boolean_array)
        node_index = np.concatenate(node_index)

        tran_symm_nodes = {}
        for index in range(0, len(node_index)):
            node_id = id_list[node_index[index]]
            node = fem.nodes[node_id]
            tran_symm_nodes[node.id] = node

        fem.add_node_group(3, tran_symm_nodes)

    def generate_boundary_conditions(self, fem: GeoGrillageFEM, symmetry: AOS):
        """
        Method generates all nodal groups for loads and boundary conditions
        based on Axis of Symmetry.
        """
        self.pinned_bc_node_group(fem)
        lc_id = 1
        fem.add_boundary_condition(1, lc_id, [0], [3], 1)

        if symmetry is AOS.LONGITUDINAL:
            self.long_symm_bc_node_group(fem)
            dof = [2, 4, 6]
            dof_val = [0.0, 0.0, 0.0]
            fem.add_boundary_condition(2, lc_id, dof_val, dof, 2)

        elif symmetry is AOS.TRANSVERSE:
            self.tran_symm_bc_node_group(fem)
            dof = [1, 5, 6]
            dof_val = [0, 0, 0]
            fem.add_boundary_condition(3, lc_id, dof_val, dof, 3)

        elif symmetry is AOS.BOTH:
            self.long_symm_bc_node_group(fem)
            self.tran_symm_bc_node_group(fem)
            dof = [2, 4, 6]
            dof_val = [0.0, 0.0, 0.0]
            fem.add_boundary_condition(2, lc_id, dof_val, dof, 2)
            dof = [1, 5, 6]
            dof_val = [0, 0, 0]
            fem.add_boundary_condition(3, lc_id, dof_val, dof, 3)

    @staticmethod
    def generate_pressure_load(fem: GeoGrillageFEM, pressure):
        lc_id = 1
        fem.addLoadCase(lc_id, "Pressure LoadCase")
        fem.add_element_group(4, fem.plate_elements)
        fem.add_pressure_load(1, lc_id, pressure, 4)

    def generate_loadcase(self, fem: GeoGrillageFEM, symmetry, pressure):
        self.generate_boundary_conditions(fem, symmetry)
        self.generate_pressure_load(fem, pressure)

    def generate_grillage_mesh(self, name, ebs, eweb, eaf, far, par, dpar):
        pass


class MeshVariantV1(GrillageMesh):
    def __init__(self, grillage: Grillage, axis_of_symm_override: AOS = None):
        super().__init__(grillage, axis_of_symm_override)

    def generate_psm_mesh(self, size: MeshSize, fem: GeoGrillageFEM):
        n_id = fem.id_node_count
        e_id = fem.id_element_count

        for segment in self.mesh_extent.full_segments.values():
            sm = SegmentMeshV1(size, segment, n_id, e_id)
            n_id, e_id = sm.generate_mesh(fem)

        for segment in self.mesh_extent.half_segments.values():
            sm = SegmentMeshV1(size, segment, n_id, e_id)
            n_id, e_id = sm.generate_mesh(fem)

    def generate_grillage_mesh(self, name, ebs, eweb, eaf, far, par, dpar):
        """
        :param name: Mesh name.
        :param ebs: Number of elements between stiffeners.
        :param eweb: Number of elements representing the web of a primary
            supporting member along its height.
        :param eaf: Number of elements across primary supporting member flange.
        :param far: Maximum PSM flange aspect ratio.
        :param par: Maximum plate and PSM web aspect ratio.
        :param dpar: Desired plating aspect ratio, less than the maximum.
        :return: GeoGrillageFEM object.
        """
        fem = GeoGrillageFEM(name)
        size = ElementSizeV1(self.mesh_extent, ebs, eweb, eaf, far, par, dpar)
        size.calculate_mesh_dimensions()
        self.generate_FEM_property(fem)
        self.generate_plate_mesh(size, fem)
        self.generate_psm_mesh(size, fem)
        fem.merge_coincident_nodes()
        fem.merge_coincident_elements()

        print("Mesh V1 generation complete.")
        return fem


class MeshVariantV2(GrillageMesh):
    def __init__(self, grillage: Grillage, axis_of_symm_override: AOS = None):
        super().__init__(grillage, axis_of_symm_override)

    def generate_psm_mesh(self, size: MeshSize, fem: GeoGrillageFEM):
        n_id = fem.id_node_count
        e_id = fem.id_element_count

        for segment in self.mesh_extent.full_segments.values():
            sm = SegmentMeshV2(size, segment, n_id, e_id)
            n_id, e_id = sm.generate_mesh(fem)

        for segment in self.mesh_extent.half_segments.values():
            sm = SegmentMeshV2(size, segment, n_id, e_id)
            n_id, e_id = sm.generate_mesh(fem)

    def generate_grillage_mesh(self, name, ebs, eweb, eaf, far, par, dpar):
        """
        :param name: Mesh name.
        :param ebs: Number of elements between stiffeners.
        :param eweb: Number of elements representing the web of a primary
            supporting member along its height.
        :param eaf: Number of elements across primary supporting member flange.
        :param far: Maximum PSM flange aspect ratio.
        :param par: Maximum plate and PSM web aspect ratio.
        :param dpar: Desired plating aspect ratio, less than the maximum.
        :return: GeoGrillageFEM object.
        """
        fem = GeoGrillageFEM(name)
        size = ElementSizeV2(self.mesh_extent, ebs, eweb, eaf, far, par, dpar)
        size.calculate_mesh_dimensions()
        self.generate_FEM_property(fem)
        self.generate_plate_mesh(size, fem)
        self.generate_psm_mesh(size, fem)
        fem.merge_coincident_nodes()
        fem.merge_coincident_elements()

        print("Mesh V2 generation complete.")
        return fem
