from typing import List, Dict
import numpy as np
import math
import  openmesh as om

def rotation_matrix_from_vectors(vec1, vec2):
    """ Find the rotation matrix that aligns vec1 to vec2
    :param vec1: A 3d "source" vector
    :param vec2: A 3d "destination" vector
    :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    """
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    return rotation_matrix

class Cylinder():
    def __init__(self,baseRadius, topRadius, height, sectors,stacks):
        super().__init__()
        self.vertices: List[int] = []
        self.indices: List[int] = []
        self.indices.clear()
        self.baseRadius = baseRadius;
        self.topRadius = topRadius;
        self.height = height;
        self.sectorCount = sectors;
        self.stackCount = stacks;
        self._unitCircleVertices: List[float] =[]
        self._buildUnitCircleVertices()

    def getmesh(self,vbase:np.ndarray, vdir:np.ndarray):
        mesh = om.TriMesh()
        #translate
        tr_mat = np.ones((4,4))
        tr_mat[0:2,3] = vbase
        #rotate
        vdirorigin = np.array([0, 0, 1])
        rot_mat = rotation_matrix_from_vectors(vdirorigin,vdir)
        tr_mat[0:2, 0:2] = rot_mat



    def _buildUnitCircleVertices(self):
        sectorCount = self.sectorCount
        PI = math.pi
        sectorStep = 2 * PI / sectorCount;

        self._unitCircleVertices.clear()
        for i in range (sectorCount+1):
            sectorAngle = i * sectorStep;
            self._unitCircleVertices.append(math.cos(sectorAngle)); # x
            self._unitCircleVertices.append(math.sin(sectorAngle)); # y
            self._unitCircleVertices.append(0); # z

    def _addVertex(self, x, y, z):
        self.vertices.append(np.array([x,y,z]))
        #self.vertices.append(x)
        #self.vertices.append(y)
        #self.vertices.append(z)

    def _clearArrays(self):
        self.vertices.clear()
        self.indices.clear()

    def _buildVerticesSmooth(self):
        baseRadius=self.baseRadius
        topRadius = self.topRadius
        height = self.height
        stackCount = self.stackCount
        self._buildUnitCircleVertices()
        sectorCount = self.sectorCount
        # clear memory of prev arrays
        self._clearArrays()

        x= y = z =0                                 # vertex position
        #float s, t;                                     # texCoord
        radius=0                                   # radius for each stack

        # put vertices of side cylinder to array by scaling unit circle
        for i in range(stackCount):
            z = -(height * 0.5) + float(i) / stackCount * height      # vertex position z
            radius = baseRadius + float(i) / stackCount * (topRadius - baseRadius)     # lerp
            t = 1.0 - float(i) / stackCount;   # top-to-bottom
            k=0
            for j in range(sectorCount+1):
                x = self._unitCircleVertices[k]
                y = self._unitCircleVertices[k + 1]
                self._addVertex(x * radius, y * radius, z)  # position
                k=k+3



        # remember where the base.top vertices start
        baseVertexIndex = len(self.vertices) / 3

        # put vertices of base of cylinder
        z = -height * 0.5
        self._addVertex(0, 0, z)
        j=0
        for i in range(sectorCount):
            x = self._unitCircleVertices[j]
            y = self._unitCircleVertices[j + 1]
            self._addVertex(x * baseRadius, y * baseRadius, z)
            j+=3

        # remember where the base vertices start
        topVertexIndex = len(self.vertices) / 3

        # put vertices of top of cylinder
        z = height * 0.5
        self._addVertex(0, 0, z)
        for i in range(sectorCount):
            x = self._unitCircleVertices[j];
            y = self._unitCircleVertices[j + 1];
            self._addVertex(x * topRadius, y * topRadius, z);
            j += 3

        # put dofs for sides
        for i in range(sectorCount):
            k1 = i * (sectorCount + 1);     # beginning of current stack
            k2 = k1 + sectorCount + 1;      # beginning of next stack

            for j in range(sectorCount):
                # 2 trianles per sector
                self.addIndices(k1, k1 + 1, k2);
                self.addIndices(k2, k1 + 1, k2 + 1);
                k1 += 1
                k2 += 1

        # remember where the base dofs start
        baseIndex = len(self.indices)

        # put dofs for base
        k = baseVertexIndex + 1
        for i in range(sectorCount):
            if(i < (sectorCount - 1)):
                self.addIndices(baseVertexIndex, k + 1, k);
            else:    # last triangle
                self.addIndices(baseVertexIndex, baseVertexIndex + 1, k);
            k+=1


        # remember where the top dofs start
        topIndex = len(self.indices);
        k = topVertexIndex + 1
        for i in range(sectorCount):
            if(i < (sectorCount - 1)):
                self.addIndices(topVertexIndex, k, k + 1);
            else:
                self.addIndices(topVertexIndex, k, topVertexIndex + 1);
            k += 1
    # generate interleaved vertex array as well
class Icosphere:
    def __init__(self, r,num_subdiv):
        self._r = r
        self._alpha = np.arctan(0.5)
        self._indices: List[int] = []
        self._vertices: List[np.ndarray] = []
        self._create()
        max_subdiv = 5
        num_subdiv= min(max_subdiv,num_subdiv)
        for i in range(num_subdiv):
            self._subdivide()

    def _upper_vertex(self, i):
        return self.octahedron_vertex(self._r, self._alpha, i)

    def _lower_vertex(self, i):
        return self.octahedron_vertex(self._r, -self._alpha, i)

    @staticmethod
    def octahedron_vertex(r, alpha, i):
        x = r * np.cos(alpha) * np.cos(np.radians(72.0) * i)
        y = r * np.cos(alpha) * np.sin(np.radians(72.0) * i)
        z = r * np.sin(alpha)
        return np.array([x,y,z])



    def _create(self):
        vertices = self._vertices
        indices = self._indices
        n = np.array([0, 0, self._r])  # north pole
        s = np.array([0, 0, -self._r])  # south pole
        vertices.append(n)
        vertices.append(s)
        # Start from top - 5 triangles
        itop=len(vertices)
        for i in range(6):
            vertices.append(self._upper_vertex(i))
        for i in range(5):
            indices.append([0,itop+i,itop+i+1])
        # Middle part - 10 triangles
        imid = len(vertices)
        for i in range(6):
            vertices.append(self._upper_vertex(i))
            vertices.append(self._lower_vertex(i))
        for i in range(5):
            ii=i*2
            indices.append([imid + ii, imid + ii + 1, imid + ii + 3])
            indices.append([imid + ii, imid + ii + 3, imid + ii + 2])

        # Lower part - 5 triangles
        ilow = len(vertices)
        for i in range(6):
            vertices.append(self._lower_vertex(i))
        for i in range(5):
            indices.append([1,ilow+i+1,ilow+i])


    @staticmethod
    def v3_length(v3):
        return np.sqrt(v3[0] ** 2 + v3[1] ** 2 + v3[2] ** 2)


    def _subdivide(self):

        subdivided_triangles = []
        subdivided_indices = []
        for triangle in self._indices:

            id0=triangle[0]
            id1 = triangle[1]
            id2 = triangle[2]
            p0 = self._vertices[triangle[0]]
            p1 = self._vertices[triangle[1]]
            p2 = self._vertices[triangle[2]]
            v0 = (p0+p1)/2.0 #0
            v0 *= self._r / np.linalg.norm(v0)
            v1 = (p1 + p2) / 2.0  # 0
            v1 *= self._r / np.linalg.norm(v1)
            v2 = (p2 + p0) / 2.0  # 0
            v2 *= self._r / np.linalg.norm(v2)
            idn0 = len(self._vertices)
            self._vertices.append(v0)
            self._vertices.append(v1)
            self._vertices.append(v2)

            subdivided_indices.append([idn0,idn0+1,idn0+2])
            subdivided_indices.append([id0, idn0, idn0 + 2])
            subdivided_indices.append([id1, idn0+1, idn0])
            subdivided_indices.append([id2, idn0 + 2, idn0+1])
        self._indices = subdivided_indices



    def append_mesh_at_center(self,mesh:om.PolyMesh, center:np.ndarray):
        nv = mesh.n_vertices()
        vertices = self._vertices+center
        indices = np.array(self._indices)+nv
        idf_start = mesh.n_faces()
        mesh.add_vertices(vertices)
        mesh.add_faces(indices)
        idf_end = mesh.n_faces()

        return np.arange(idf_start,idf_end)