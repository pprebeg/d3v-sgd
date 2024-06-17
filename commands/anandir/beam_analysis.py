from typing import Dict, List,Tuple
import numpy as np
from collections import Counter
from enum import Enum
import matplotlib.pyplot as plt

class BeamBC(Enum):
    SS_SS=1
    C_C=2
    SS_C = 3
    C_SS = 4

class BeamProperty:
    def __init__(self,id):
        self._id = id

    def getW(self):
        pass

    def getIy(self):
        pass

    def getE(self):
        return 0.0
    @property
    def id(self):
        return self._id

class IBeamProperty(BeamProperty):
    def __init__(self, id, hw, tw, bf, tf, bp, tp, E):
        super().__init__(id)
        self._hw = hw
        self._tw = tw
        self._bf = bf
        self._tf = tf
        self._bp = bp
        self._tp = tp
        self._E = E



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

    def getMx(self, x, L, bc=BeamBC.SS_SS):
        pass

    def get_x_Mcrit(self, L, bc=BeamBC.SS_SS):
        pass

    def get_moment_clamped_beam_end(self, L, is2:bool, bc=BeamBC.SS_SS):
        pass

    def get_beam_deflection(self, x, E, I, L , bc=BeamBC.SS_SS):
        pass


class BeamLoadNoLoad(BeamLoad):
    def __init__(self, id: int, q0 = 0.0):
        super().__init__(id, q0)
        pass

    def getMx(self, x, L, bc=BeamBC.SS_SS):
        if isinstance(x, np.ndarray):
            return np.zeros(x.size)
        else:
            return 0.0

    def get_x_Mcrit(self, L, bc=BeamBC.SS_SS):
        return 0.0

    def get_moment_clamped_beam_end(self, L, is2:bool, bc=BeamBC.SS_SS):
        return 0.0
    def get_beam_deflection(self, x, E, I, L , bc=BeamBC.SS_SS):
        return 0.0

class BeamLoadConstContLoad(BeamLoad):
    def __init__(self, id: int, q0):
        super().__init__(id, q0)
        pass
    @property
    def q0(self):
        return self._q0

    def getMx(self, x, L, bc=BeamBC.SS_SS):
        return (self.q0 * x) / 2. * (L - x)

    def get_x_Mcrit(self, L, bc=BeamBC.SS_SS):
        return L / 2

    def get_moment_clamped_beam_end(self, L, is2:bool, bc=BeamBC.SS_SS):
        return self.q0*L**2.0/12.0

    def get_beam_deflection(self, x, E, I, L , bc=BeamBC.SS_SS):
        return ((self.q0*x)/(24.0*E*I))*(L**3.0 - 2.0*L*x**2.0+x**3.0)


class BeamLoadTriangleContLoad(BeamLoad):
    def __init__(self, id: int, q0: float, is2: bool = True):
        super().__init__(id, q0)
        self._is2 = is2 # true if q0 is on second node

    def getMx(self, x, L, bc=BeamBC.SS_SS):
        W = self._q0 * L / 2.
        if not self._is2:
            x = L - x
        Mx = W * x * (L ** 2.0 - x ** 2.0) / (3.0 * L ** 2.0)
        return Mx

    def get_x_Mcrit(self, L, bc=BeamBC.SS_SS):
        x = 0.5774 * L
        if not self._is2:
            x = L - x
        return x

    def get_moment_clamped_beam_end(self, L, is2:bool, bc=BeamBC.SS_SS):
        if (is2 and self._is2) or ((not is2) and (not self._is2)):
            return self._q0*L**2.0 / 20.0
        else:
            return self._q0 * L ** 2.0 / 30.0

class BeamLoadPointLoadAnyPoint(BeamLoad):
    def __init__(self, id: int, P,a):
        super().__init__(id, P)
        self._a = a
        pass

    @property
    def a(self):
        return self._a

    @property
    def P(self):
        return self._q0

    def getMx(self, x, L, bc=BeamBC.SS_SS):
        a = self.a
        b = L - a
        if isinstance(x, (int, float)):
            if x < a:
                Mx = self.P * b * x / L
            else:
                Mx = self.P * a / L * (L - x)
        elif isinstance(x, np.ndarray):
            Mx = np.where(x < a, self.P * b * x / L, self.P * a / L * (L - x))
        else:
            raise ValueError("Unsupported input type for 'x'.")

        return Mx

    def get_x_Mcrit(self, L, bc=BeamBC.SS_SS):
        return self.a

    def get_moment_clamped_beam_end(self, L, is2:bool, bc=BeamBC.SS_SS):
        a = self.a
        b = L - a
        if is2:
            Mx = self.P * b **2.0* a/ L**2.0
        else:
            Mx = self.P * b  * a** 2.0 / L ** 2.0
        return Mx

    @staticmethod
    def get_any_beam_deflection(a,P, x, E, I, L, bc=BeamBC.SS_SS ):
        if bc == BeamBC.SS_SS:
            b = L - a
            if x < a:
                wx = (P * b * x) / (6 * L * E * I) * (L**2 - b**2 - x**2)
            else:
                wx = (P * a * (L - x)) / (6 * L * E * I) * (2*L*x - x**2 - a**2)
            return wx
        else:
            return 0.0
    def get_beam_deflection(self, x, E, I, L, bc=BeamBC.SS_SS ):
        a = self.a
        P = self.P
        return BeamLoadPointLoadAnyPoint.get_any_beam_deflection(a,P, x, E, I, L,bc)

class BeamLoadReaction(BeamLoadPointLoadAnyPoint):
    def __init__(self, id: int, R,a):
        super().__init__(id, R,a)

    @property
    def R(self):
        return self._q0

class Node:
    def __init__(self, id: int, x=0.0, y=0.0, z=0.0,cords:np.ndarray=None):
        self._id = id
        if cords is None:
            self._cords: np.ndarray = np.array([x, y, z])
        else:
            self._cords = cords

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
    npdiag=40
    def __init__(self, id: int, prop=None, node1=None, node2=None):
        self._id = id
        self._nodes: List[Node] = []
        if node1 is not None:
            self._nodes.append(node1)
        if node2 is not None:
            self._nodes.append(node2)
        self._loads:List[BeamLoad] = []
        self._elastic_reaction_loads: List[BeamLoad] = []
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

    def add_elastic_reaction_load(self, load: BeamLoad):
        self._elastic_reaction_loads.append(load)

    def clear_elastic_reaction_loads(self):
        self._elastic_reaction_loads.clear()

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

    def getSigmax_crit_from_diagram(self):
        Mcrit = self.getMx_crit_from_diagram()
        W = self._prop.getW()
        return Mcrit / W

    def getMx_crit(self):
        Mxcrit = 0
        xcand = []
        for load in self._loads:
            xcand.append(load.get_x_Mcrit(self.L))

        for x in xcand:
            Mx = abs(self.get_Mx(x) - self._M1)
            if Mxcrit < Mx:
                Mxcrit = Mx
        if Mxcrit < abs(self._M1):
            Mxcrit = self._M1
        if Mxcrit < abs(self._M2):
            Mxcrit = self._M2
        return Mxcrit

    def getMx_crit_from_diagram(self):
        x = np.linspace(0.0,self.L,AnalyticBeam.npdiag)
        Mx = np.abs(self.get_Mx(x))
        Mxcrit = np.max(Mx)
        return Mxcrit

    def get_Mx(self, x):
        Mx = 0
        for load in self._loads:
            Mx = Mx + load.getMx(x, self.L)
        for load in self._elastic_reaction_loads:
            Mx = Mx + load.getMx(x, self.L)
#        if isinstance(x, np.ndarray):
#            plt.plot(x, Mx)
#            plt.show()
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

    def find_intersection_with_beam(self, beam):
        """
        Returns:
        Intersection point (x, y, z), or None if there is no unique intersection.
        """
        P1, D1 = self.node1.cords, self.node2.cords - self.node1.cords
        P2, D2 = beam.node1.cords, beam.node2.cords - beam.node1.cords

        # Matrix A composed of direction vectors D1 and D2 and the difference between start points
        A = np.column_stack((D1, -D2, P2 - P1))

        # Check if there is a solution
        if np.linalg.matrix_rank(A[:, :2]) < 2 or np.linalg.matrix_rank(A) != np.linalg.matrix_rank(A[:, :2]):
            return None

        # Solve for parameters t and s
        params, _, _, _ = np.linalg.lstsq(A[:, :2], A[:, 2], rcond=None)
        t, s = params[:2]

        # Compute the intersection point
        intersection_point = P1 + t * D1
        return intersection_point

    def beam_deflection_unit_point_load(self, b, a, x):
        unit_defl = BeamLoadPointLoadAnyPoint.get_any_beam_deflection(a,1.0,x,self._prop.getE(),self._prop.getIy(),self.L)
        return unit_defl

    def beam_deflection_original_loads(self, x):
        total=0.0
        for load in self._loads:
            total += load.get_beam_deflection(x, self._prop.getE(), self._prop.getIy(), self.L)
        return total

    def beam_deflection_original_and_elastic(self, x):
        total=0.0
        for load in self._loads:
            total += load.get_beam_deflection(x, self._prop.getE(), self._prop.getIy(), self.L)
        for load in self._elastic_reaction_loads:
            total += load.get_beam_deflection(x, self._prop.getE(), self._prop.getIy(), self.L)
        return total

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
        self._cross_beams:Dict[AnalyticBeam] = {}
        self._main_beams:Dict[AnalyticBeam] = {}
        self._intersections:List[Tuple[int]] = []
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


    def identify_intersections(self):
        self._intersections.clear()
        for keys_main, main in self._main_beams.items():
            for keys_cross, cross in self._cross_beams.items():
                    P = main.find_intersection_with_beam(cross)
                    if P is not None:
                        self._intersections.append((keys_main,keys_cross))
        print(self._intersections)

    def assemble_matrix(self):
        n = len(self._intersections)
        kij = np.zeros((n, n))
        bi = np.zeros(n)
        for i in range(n):
            (idmi,idci) = self._intersections[i]
            mbi:AnalyticBeam = self._main_beams[idmi]
            cbi:AnalyticBeam = self._cross_beams[idci]
            P = mbi.find_intersection_with_beam(cbi)
            xlocm = np.linalg.norm(P - mbi.node1.cords)
            xlocc = np.linalg.norm(P - cbi.node1.cords)
            bi[i] = cbi.beam_deflection_original_loads(xlocc)
            bi[i] +=mbi.beam_deflection_original_loads(xlocm)
            for j in range(n):
                (idmj, idcj) = self._intersections[j]
                if idmi == idmj:
                    cbj:AnalyticBeam = self._cross_beams[idcj]
                    P = mbi.find_intersection_with_beam(cbj)
                    a = np.linalg.norm(P - mbi.node1.cords)
                    b=mbi.L-a
                    kij[i,j] += mbi.beam_deflection_unit_point_load(b,a,xlocm)
                if idci == idcj:
                    mbj:AnalyticBeam = self._main_beams[idmj]
                    P = cbi.find_intersection_with_beam(mbj)
                    a = np.linalg.norm(P - cbi.node1.cords)
                    b = cbi.L - a
                    kij[i,j] += cbi.beam_deflection_unit_point_load(b,a,xlocc)
        return kij,bi

    def check_equivalence_of_displacements_at_intersection(self):
        n = len(self._intersections)
        wm = np.zeros(n)
        wc = np.zeros(n)
        for i in range(n):
            (idmi,idci) = self._intersections[i]
            mbi:AnalyticBeam = self._main_beams[idmi]
            cbi:AnalyticBeam = self._cross_beams[idci]
            P = mbi.find_intersection_with_beam(cbi)
            xlocm = np.linalg.norm(P - mbi.node1.cords)
            xlocc = np.linalg.norm(P - cbi.node1.cords)
            wc[i] = cbi.beam_deflection_original_and_elastic(xlocc)
            wm[i] = mbi.beam_deflection_original_and_elastic(xlocm)
        is_close = np.allclose(wm,wc)
        return is_close,wc,wm,self._intersections

    def calculate_reaction_on_intersections(self,kij,bi):

        R = np.linalg.solve(kij,bi)
        n= len(bi)
        for i in range(n):
            self._reactions[i] = R[i]


    def set_elastic_reaction_loads(self):
        n = len(self._intersections)
        idrl = 0
        for beam in self._beams.values():
            beam.clear_elastic_reaction_loads()
        for i in range(n):
            (idmi, idci) = self._intersections[i]
            mbi: AnalyticBeam = self._main_beams[idmi]
            cbi: AnalyticBeam = self._cross_beams[idci]
            P = mbi.find_intersection_with_beam(cbi)
            am = np.linalg.norm(P - mbi.node1.cords)
            mbi.add_elastic_reaction_load(BeamLoadReaction(idrl, self._reactions[i], am))
            ac = np.linalg.norm(P - cbi.node1.cords)
            cbi.add_elastic_reaction_load(BeamLoadReaction(idrl, -self._reactions[i], ac))


    def calculate_and_apply_elastic_reactions(self):
        self.beam_arrangement_main()
        self.beam_arrangement_cross()
        self.identify_intersections()
        kij,bi = self.assemble_matrix()
        self.calculate_reaction_on_intersections(kij,bi)
        self.set_elastic_reaction_loads()

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
    grillage.calculate_and_apply_elastic_reactions()
    isEq,wc,wm,itsec = grillage.check_equivalence_of_displacements_at_intersection()
    print ('Main beams and cross beams deformation equivalent on interesctions = {:}'.format(isEq))
    print(grillage._reactions)
    for i in range(len(wc)):
        print('Deformation on intersection {0:}, main beam =  {1:.3f}, cross beam = {2:.3f}'.format(itsec[i],wm[i],wc[i]))
    for beam in grillage._main_beams.values():
        Mx = beam.getMx_crit_from_diagram()
        print('Main beam {0:.1f} Mx = {1:.1f}'.format(beam.id, Mx))
    for beam in grillage._cross_beams.values():
        Mx = beam.getMx_crit_from_diagram()
        print('Cross beam {0:.1f} Mx = {1:.1f}'.format(beam.id, Mx))
    for beam in grillage._main_beams.values():
        Sx = beam.getSigmax_crit_from_diagram()
        print('Main beam {0:.1f} Sigmax = {1:.1f}'.format(beam.id, Sx))
    for beam in grillage._cross_beams.values():
        Sx = beam.getSigmax_crit_from_diagram()
        print('Cross beam {0:.1f} Sigmax = {1:.1f}'.format(beam.id, Sx))
if __name__ == "__main__":
    _test_grillage()