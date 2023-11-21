from typing import Dict, List
import numpy as np
from collections import Counter


class BeamProperty:
    def __init__(self,id):
        self._id = id

    def getW(self):
        pass

    def getIy(self):
        pass

    def getE(self):
        return 0.0

class IBeamProperty(BeamProperty):
    def __init__(self, id, hw, tw, bf, tf, bp, tp, E):
        super.__init__(id)
        self._hw = hw
        self._tw = tw
        self._bf = bf
        self._tf = tf
        self._bp = bp
        self._tp = tp
        self._E = E

    @property
    def id(self):
        return self._id

    def getW(self):
        hw = self._hw
        tw = self._tw
        bf = self._bf
        tf = self._tf
        bp = self._bp
        tp = self._tp

        A = bf * tf + hw * tw + bp * tp
        z = ((bp * tp) * tp / 2 + hw * tw * (tp + hw / 2) + bf * tf * (tp + hw + tf / 2)) / A
        Iy = (bp * tp ** 3) / 12 + ((z - tp / 2) ** 2) * bp * tp + (tw * hw ** 3) / 12 + (
                ((tp + hw / 2) - z) ** 2) * tw * hw + (bf * tf ** 3) / 12 + (
                     ((tp + hw + tf / 2) - z) ** 2) * bf * tf
        return Iy / ((tp + hw + tf) - z)

    def getIy(self):
        hw = self._hw
        tw = self._tw
        bf = self._bf
        tf = self._tf
        bp = self._bp
        tp = self._tp

        A = bf * tf + hw * tw + bp * tp
        z = ((bp * tp) * tp / 2 + hw * tw * (tp + hw / 2) + bf * tf * (tp + hw + tf / 2)) / A
        Iy = (bp * tp ** 3) / 12 + ((z - tp / 2) ** 2) * bp * tp + (tw * hw ** 3) / 12 + (
                ((tp + hw / 2)-z) ** 2) * tw * hw + (bf * tf ** 3) / 12 + (
                     ((tp + hw + tf / 2) - z) ** 2) * bf * tf
        return Iy

    def getE(self):
        return self._E



class BeamLoad:
    def __init__(self, id: int, q0: float):
        self._id = id
        self._q0 = q0

    def getMx(self, x, L):
        pass

    def get_xcrit(self, L):
        pass

    def get_moment_clamped_beam_end(self, L, is2:bool):
        pass

    def get_beam_deflection(self, x, E, I, L ):
        pass


class BeamLoadNoLoad(BeamLoad):
    def __init__(self, id: int, q0):
        super().__init__(id, q0)
        pass

    def getMx(self, x, L):
        return 0.0

    def get_xcrit(self, L):
        return 0.0

    def get_moment_clamped_beam_end(self, L, is2:bool):
        return 0.0


class BeamLoadConstContLoad(BeamLoad):
    def __init__(self, id: int, q0):
        super().__init__(id, q0)
        pass

    def getMx(self, x, L):
        return (self._q0 * x) / 2. * (L - x)

    def get_xcrit(self, L):
        return L / 2

    def get_moment_clamped_beam_end(self, L, is2:bool):
        return self._q0*L**2.0/12.0

    def get_beam_deflection(self, x, E, I, L ):
        return ((self.q0*x)/(24.0*E*I))*(L**2.0 - 2.0*L*x**2.0+x**3.0)


class BeamLoadTriangleContLoad(BeamLoad):
    def __init__(self, id: int, q0: float, is2: bool = True):
        super().__init__(id, q0)
        self._is2 = is2 # true if q0 is on second node

    def getMx(self, x, L):
        W = self._q0 * L / 2.
        if not self._is2:
            x = L - x
        Mx = W * x * (L ** 2.0 - x ** 2.0) / (3.0 * L ** 2.0)
        return Mx

    def get_xcrit(self, L):
        x = 0.5774 * L
        if not self._is2:
            x = L - x
        return x

    def get_moment_clamped_beam_end(self, L, is2:bool):
        if (is2 and self._is2) or ((not is2) and (not self._is2)):
            return self._q0*L**2.0 / 20.0
        else:
            return self._q0 * L ** 2.0 / 30.0


class Node:
    def __init__(self, id: int, x=0.0, y=0.0, z=0.0):
        self._id = id
        self._cords: np.ndarray = np.array([x, y, z])

    @property
    def id(self):
        return self._id

    @property
    def cords(self):
        return self._cords

    def get_distance(self, node):
        return np.linalg.norm(node.cords - self._cords)

    @staticmethod
    def get_distance(node1, node2):
        return np.linalg.norm(node1.cords - node2.cords)



class AnalyticBeam():
    def __init__(self, id: int, prop=None, node1=None, node2=None):
        self._id = id
        self._nodes: List[Node] = []
        if node1 is not None:
            self._nodes.append(node1)
        if node2 is not None:
            self._nodes.append(node2)
        self._loads:List[BeamLoad] = []
        self._prop: BeamProperty = prop
        self._M1:float = 0.0 # internal moment on the node 1 end calculated by displacement method
        self._M2:float = 0.0 # internal moment on the node 2 end calculated by displacement method
        pass

    @property
    def id(self):
        return self._id

    @property
    def L(self):
        return Node.get_distance(self.node1, self.node2)

    @property
    def node1(self):
        return self._nodes[0]

    @property
    def node2(self):
        return self._nodes[1]

    def have_node(self, node: Node):
        if self.node1 is node:
            return True
        if self.node2 is node:
            return True
        return False

    def add_load(self, load: BeamLoad):
        self._loads.append(load)

    def set_prop(self, prop: BeamProperty):
        self._prop = prop

    def set_load(self, load):
        self._loads = load


    def check_beam(self):
        if len(self._loads) == 0:
            self.add_load(BeamLoadNoLoad())

    def getSigmax_crit(self):
        Mcrit = self.getMx_crit()
        W = self._prop.getW()
        return Mcrit / W

    def getMx_crit(self):
        Mxcrit = 0
        xcand = []
        for load in self._loads:
            xcand.append(load.get_xcrit(self.L))

        for x in xcand:
            Mx = abs(self.get_Mx(x) - self._M1)
            if Mxcrit < Mx:
                Mxcrit = Mx
        if Mxcrit < abs(self._M1):
            Mxcrit = self._M1
        if Mxcrit < abs(self._M2):
            Mxcrit = self._M2
        return Mxcrit

    def get_Mx(self, x):
        Mx = 0
        for load in self._loads:
            Mx = Mx + load.getMx(x, self.L)
        return Mx

    def get_k(self):
        k = 2 * self._prop.getIy() * self._prop.getE() / self.L
        return k

    def get_stifness_for_frame(self, current_node:Node):
        if not self.have_node(current_node):
            return None
        k = self.get_k()
        mb = 0.0
        if self.node1 is current_node:
            # node 1 is current, node 2 is other
            for load in self._loads:
                mb -= load.get_moment_clamped_beam_end(self.L,False)
            return (k, mb, self.node2.id)
        else:
            # node 2 is current, node 1 is other
            for load in self._loads:
                mb += load.get_moment_clamped_beam_end(self.L, True)
            return (k, mb, self.node1.id)

    def calculate_internal_end_moments(self, phi1, phi2):
        k = self.get_k()
        # node 1 is current, node 2 is other
        mb = 0.0
        for load in self._loads:
            mb += load.get_moment_clamped_beam_end(self.L, False)
        self._M1 = k * (2 * phi1 + phi2) + mb
        # node 2 is current, node 1 is other
        mb = 0.0
        for load in self._loads:
            mb += load.get_moment_clamped_beam_end(self.L, True)
        self._M2 = k * (2 * phi2 + phi1) - mb


    # turning a beam into a line for cramer's rule
    def line(self, p1, p2):
        A = (p1[1] - p2[1])
        B = (p2[0] - p1[0])
        C = (p1[0] * p2[1] - p2[0] * p1[1])
        return A, B, -C

    def beam_deflection(self, b, a, x):
        if a < x:
            return ((a * (self.L-x))/(6 * self._prop.getIy() * self._prop.getE())) * (2*self.L*x - x**2 - a**2) #fali l u nazivniku
        else:
            return ((b * x)/(6 * self._prop.getIy() * self._prop.getE())) * (self.L**2 - b**2 - x**2)

    def beam_CL_deflection(self, x):
        for load in self._loads:
            return ((load._q0 * x) / (24.0 * self._prop.getE() * self._prop.getIy())) * ((self.L ** 3.0) - 2.0 * self.L * x ** 2.0 + x ** 3.0)



class WebFrame():
    def __init__(self):
        self._nodes: Dict[int, Node] = {}
        self._beams: Dict[int, AnalyticBeam] = {}
        self._props: Dict[int, BeamProperty] = {}
        self._k = {}
        self._m = {}
        pass

    def assemble_stiff_matrix(self):
        null = None
        self._k.clear()
        self._m.clear()
        for idnode,node in self._nodes.items():
            row_vals = {}
            m = 0.0
            for beam in self._beams.values():
                if not beam.have_node(node):
                    continue
                res = beam.get_stifness_for_frame(node)

                if res is not null:
                    (k, mb, idother) = res
                    row_vals[idnode] = row_vals.get(idnode,0.0)+ 2.0*k
                    row_vals[idother] = row_vals.get(idother, 0.0) + k
                    m+=mb
            if len(row_vals) > 0:
                self._k[idnode] = row_vals
                self._m[idnode] = m

    def calculate_phi(self):
        n = len(self._k)
        k = np.zeros((n,n))
        m = np.zeros(n)
        ids = [-1]*n
        indexes = {}
        index = 0
        for id_row, rowvals in self._k.items():
            indexes[id_row] = index
            ids[index] = id_row
            index += 1

        for id_row,rowvals in self._k.items():
            for id_col,kxx in rowvals.items():
                irow=indexes[id_row]
                icol = indexes[id_col]
                k[irow,icol] = kxx

        for id_row,mrow in self._m.items():
            irow = indexes[id_row]
            m[irow] = mrow


        phi = np.linalg.solve(k,m)
        dphi = {}
        for i in range(n):
            dphi[ids[i]] = phi[i]

        return dphi

    def calculate_M(self):
        phi = self.calculate_phi()
        for beam in self._beams.values():
            beam.calculate_internal_end_moments(phi[beam.node1.id], phi[beam.node2.id])
    def calculate_Sigmacrit(self):
        for beam in self._beams.values():
            sx = beam.getSigmax_crit()
            print('Beam{0}: Sx_crit={1:6.1f}'.format(beam.id, sx))

    def print_Nodal_Moments(self):
        for beam in self._beams.values():
            print('Beam{0}: M1={1:10.3E}, M2={2:10.3E}'.format(beam.id,beam._M1,beam._M2))

    @property
    def props(self):
        return self._props

    @property
    def nodes(self):
        return self._nodes

    @property
    def beams(self):
        return self._beams

    def add_node(self, node: Node):
        self._nodes[node.id] = node

    def add_prop(self, prop: BeamProperty):
        self._props[prop.id] = prop

    def add_beam(self, beam: AnalyticBeam):
        self._beams[beam.id] = beam

    def get_node(self, idnode: int):
        return self._nodes[idnode]

    def get_prop(self, idprop: int):
        return self._props[idprop]

    def get_beam(self, idbeam: int):
        return self._beams[idbeam]

class GrillageAnalysis():
    def __init__(self):
        self._nodes: Dict[int, Node] = {}
        self._beams: Dict[int, AnalyticBeam] = {}
        self._props: Dict[int, BeamProperty] = {}
        self._cross_beams = {}
        self._main_beams = {}
        self._intersections_main = {}
        self._intersections_cross = {}
        self._deflection_in_every_point_main_beams = {}
        self._deflection_in_every_point_cross_beams = {}
        self._K = {}
        self._deflection_CL = {}
        self._reactions = {}
        pass

    def beam_arrangement_cross(self):
        id_beam = 0
        for beam in self._beams.values():
            if type(beam._loads[0]) != BeamLoadNoLoad:
                self._cross_beams[id_beam] = beam
                id_beam += 1

    def beam_arrangement_main(self):
        id_beam = 0
        for beam in self._beams.values():
            if type(beam._loads[0]) == BeamLoadNoLoad:
                self._main_beams[id_beam] = beam
                id_beam += 1

    # cramer's rule,
    def beam_intersection(self, main, cross):
        D = main[0] * cross[1] - main[1] * cross[0]
        Dx = main[2] * cross[1] - main[1] * cross[2]
        Dy = main[0] * cross[2] - main[2] * cross[0]
        if D != 0:
            x = Dx / D
            y = Dy / D
            return (x, y)

    def deflection_of_main_beams(self):
        idx_of_intersection_main = 1
        idx_for_one_x = 1.0
        for keys_main, main in self._main_beams.items():
            if keys_main > 0.0 and keys_main != (len(self._main_beams)-1.0):
                M = main.line(main.node1.cords, main.node2.cords)
                idx = 0.0
                for keys_cross, cross in self._cross_beams.items():
                    if keys_cross > 0.0 and keys_cross != (len(self._cross_beams)-1.0):
                        C = cross.line(cross.node1.cords, cross.node2.cords)
                        x, y = self.beam_intersection(M, C)
                        self._intersections_main[idx] = x, y
                        idx += 1.0
                for x in self._intersections_main.values():
                     idx_for_one_x0 = idx_for_one_x
                     count = 1.0
                     a = main.L / (len(self._cross_beams) - 1.0)
                     deflection_for_one_x_main = {}
                     while count <= len(self._cross_beams) - 2.0:
                          b = main.L - a
                          w_on_cords = main.beam_deflection(b, a, x[0])
                          a = a + main.L/(len(self._cross_beams) - 1.0)
                          count += 1.0
                          deflection_for_one_x_main[idx_for_one_x0] = w_on_cords
                          idx_for_one_x0 += 1.0
                     self._deflection_in_every_point_main_beams[idx_of_intersection_main] = deflection_for_one_x_main
                     idx_of_intersection_main += 1
                idx_for_one_x = idx_for_one_x0
        #print(self._deflection_in_every_point_main_beams)

    def deflection_of_cross_beams(self):
        idx_for_one_x = 1.0
        for keys_cross, cross in self._cross_beams.items():
            if keys_cross > 0 and keys_cross != (len(self._cross_beams) - 1):
                C = cross.line(cross.node1.cords, cross.node2.cords)
                idx = 0.0
                for keys_main, main in self._main_beams.items():
                    if keys_main > 0 and keys_main != (len(self._main_beams) - 1):
                        M = main.line(main.node1.cords, main.node2.cords)
                        x, y = self.beam_intersection(M, C)
                        self._intersections_cross[idx] = x, y
                        idx += 1.0
                for x in self._intersections_cross.values():
                    idx_for_one_x0 = idx_for_one_x
                    count = 1.0
                    a = cross.L / (len(self._main_beams) - 1.0)
                    deflection_for_one_x_cross = {}
                    w_on_CL = cross.beam_CL_deflection(x[1])
                    self._deflection_CL[keys_cross] = w_on_CL
                    while count <= len(self._main_beams) - 2.0:
                        b = (cross.L - a)
                        w_on_R = cross.beam_deflection(b, a, x[1])
                        a = a + cross.L / (len(self._main_beams) - 1.0)
                        count += 1.0
                        deflection_for_one_x_cross[idx_for_one_x0] = w_on_R
                        idx_for_one_x0 += (len(self._cross_beams) - 2.0)
                    self._deflection_in_every_point_cross_beams[keys_cross] = deflection_for_one_x_cross
                    keys_cross = keys_cross + (len(self._cross_beams) - 2)
                idx_for_one_x += 1.0
        print(self._deflection_in_every_point_cross_beams)
        print(self._deflection_CL)

    def assemble_matrix(self):
        idx_row = 1
        for key_main, main_deflection_for_one_x_main in self._deflection_in_every_point_main_beams.items():
             for key_cross, main_deflection_for_one_x_cross in self._deflection_in_every_point_cross_beams.items():
                 if key_main == key_cross:
                     value_main = Counter(main_deflection_for_one_x_main)
                     value_cross = Counter(main_deflection_for_one_x_cross)
                     res = value_main + value_cross
                     self._K[idx_row] = res
             idx_row += 1.0

    def calculate_reaction_on_intersections(self):
        n = len(self._K)
        k = np.zeros((n,n))
        m = np.zeros(n)
        ids = [-1]*n
        indexes = {}
        index = 0
        for id_row, rowvals in self._K.items():
            indexes[id_row] = index
            ids[index] = id_row
            index += 1

        for id_row,rowvals in self._K.items():
            for id_col,kxx in rowvals.items():
                irow=indexes[id_row]
                icol = indexes[id_col]
                k[irow,icol] = kxx

        for id_row,mrow in self._deflection_CL.items():
            irow = indexes[id_row]
            m[irow] = mrow

        R = np.linalg.solve(k,m)
        for i in range(n):
            self._reactions[ids[i]] = R[i]

        return self._reactions

    def deflection_value_main(self, location):
        main_beam_deflection_value = []
        sum_box_main = []
        n = len(self._cross_beams) - 2.0
        bound = 0
        for keys_main, main in self._main_beams.items():
             if keys_main > 0.0 and keys_main != (len(self._main_beams)-1.0):
                 part_Reactions = {k: v for k, v in self._reactions.items() if k <= n and k > bound}
                 a = main.L / (len(self._cross_beams) - 1.0)
                 for _R in part_Reactions.values():
                     b = main.L - a
                     w_value_one = main.beam_deflection(b, a, location) * _R
                     a = a + main.L / (len(self._cross_beams) - 1.0)
                     sum_box_main.append(w_value_one)
                 Sum = sum(sum_box_main)
                 main_beam_deflection_value.append(Sum)
                 n += n
                 bound+= len(self._cross_beams) - 2.0
        return main_beam_deflection_value

    def deflection_value_cross(self,location):
        cross_beam_deflection_value = []
        sum_box_cross =[]
        n = len(self._cross_beams) - 1.0
        bound = 1
        for keys_cross, cross in self._cross_beams.items():
            if keys_cross > 0.0 and keys_cross != (len(self._cross_beams) - 1.0):
                part_Reactions_cross = {key: self._reactions[key] for key in self._reactions.keys() & {bound, n}}
                a = cross.L / (len(self._main_beams) - 1.0)
                for _R_cross in part_Reactions_cross.values():
                    b = (cross.L - a)
                    w_value_one_cross = cross.beam_deflection(b, a, location) * _R_cross
                    a = a + cross.L / (len(self._main_beams) - 1.0)
                    sum_box_cross.append(w_value_one_cross)
                Sum_cross = sum(sum_box_cross)
                sum_Reaction = cross.beam_CL_deflection(location) - Sum_cross
                cross_beam_deflection_value.append(sum_Reaction)
                n += 1
                bound += 1
        return cross_beam_deflection_value

    def calculate_reactions(self):
        self.beam_arrangement_main()
        self.beam_arrangement_cross()
        self.deflection_of_main_beams()
        self.deflection_of_cross_beams()
        self.assemble_matrix()
        R = self.calculate_reaction_on_intersections()
        return R

    def read_file(self, file_path):
        with open(file_path, 'r') as f:
            params = {}
                # params['type'] = 'I'
                # params[b] = 2
            lines = f.readlines()
            iline = 0
            line = lines[iline]
            split = line.split(' ')
            n_nodes = int(split[0])  # 4
            n_props = int(split[1])  # 4
            n_elem = int(split[2])   # 4
            n_loads = int(split[3])  # 4


            for i in range(1, n_nodes + 1):  # 1,2,3,4
                line = lines[i]
                split = line.split(' ')
                node = Node(int(split[0]), float(split[1]), float(split[2]), float(split[3]))
                self.add_node(node)

            for i in range(n_nodes + 1, n_nodes + n_props + 1):  # 5,6,7
                line = lines[i]
                split = line.split(' ')
                fsplit = []
                for idx in range(4, len(split), 2):
                    fsplit.append(float(split[idx]))

                beam_property = IBeamProperty(int(split[0]), fsplit[0], fsplit[1], fsplit[2], fsplit[3], fsplit[4], fsplit[5], fsplit[6])
                self.add_prop(beam_property)

            for i in range(n_nodes + n_props + 1, n_nodes + n_props + n_elem + 1):  # 8,9,10
                line = lines[i]
                split = line.split(' ')
                beam = AnalyticBeam(id=int(split[0]), prop=self.get_prop(int(split[3])), node1=self.get_node(int(split[1])), node2=self.get_node(int(split[2])))
                self.add_beam(beam)

            for i in range(n_nodes + n_props + n_elem + 1, n_nodes + n_props + n_elem + n_loads + 1):  # 11, 12, 13
                line = lines[i]
                split = line.split(' ')
                id_beam = int(split[2])
                if split[4] == 'CL':
                    load = BeamLoadConstContLoad(int(split[0]), float(split[6]))
                elif split[4] == 'TL':
                    is2 = int(split[7])==2
                    load = BeamLoadTriangleContLoad(int(split[0]), float(split[6]),is2)
                else:
                    load = BeamLoadNoLoad(int(split[0]), float(split[6]))
                self.get_beam(id_beam).add_load(load)


    @property
    def props(self):
        return self._props

    @property
    def nodes(self):
        return self._nodes

    @property
    def beams(self):
        return self._beams

    def add_node(self, node: Node):
        self._nodes[node.id] = node

    def add_prop(self, prop: BeamProperty):
        self._props[prop.id] = prop

    def add_beam(self, beam: AnalyticBeam):
        self._beams[beam.id] = beam

    def get_node(self, idnode: int):
        return self._nodes[idnode]

    def get_prop(self, idprop: int):
        return self._props[idprop]

    def get_beam(self, idbeam: int):
        return self._beams[idbeam]


def _test_grillage():
    print('Test Grillage Analysis')
    grillage = GrillageAnalysis()
    grillage.read_file("grillage_analysis_test.txt")
    grillage.calculate_reactions()
    print(grillage.deflection_value_main(2))  # x koordinata lokalna s crteza
    print(grillage.deflection_value_cross(2))  # x koordinata lokalna s crteza
    print(grillage._reactions)
if __name__ == "__main__":
    _test_grillage()