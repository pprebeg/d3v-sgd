import numpy as np
import openmesh as om
import copy
#import grillage_model as gm
from grillage.grillage_model import BeamDirection as bd
from grillage.grillage_model import Plate,PrimarySuppMem,Grillage,Segment
from grillage.grillage_model import TBeamProperty,BulbBeamProperty,HatBeamProperty,LBeamProperty,FBBeamProperty
from core.geometry import Geometry
from typing import List,Dict


class Resolution_Value_Error(Exception):
    pass
class Invalid_Segment_id_Error(Exception):
    pass
class Invalid_Stiff_dir(Exception):
    pass

def join_tria_meshes(mesh1:om.TriMesh,mesh2:om.TriMesh):
    if mesh1.n_faces() ==0:
        return  mesh2
    elif mesh2.n_faces() ==0:
        return  mesh1
    fvs1 = mesh1.fv_indices()
    points1 = mesh1.points()
    fvs2 = mesh2.fv_indices()
    points2 = mesh2.points()
    npts1 = points1.shape[0]
    fvs2 = fvs2+npts1
    fvs = np.concatenate((fvs1,fvs2), axis=0)
    points = np.concatenate((points1,points2), axis=0)
    return om.TriMesh(points,fvs)

def offset_line(line_points, line_evi, offset=3):
    normal = np.array([0, 1.0, 0])
    epoints = line_points[line_evi]
    ep_dir_vectors = epoints[:, 1] - epoints[:, 0]
    outward_dir = np.cross(ep_dir_vectors, normal)
    outward_dir = outward_dir[1:] + outward_dir[:-1]
    outward_dir = outward_dir / (((outward_dir ** 2).sum(-1) ** 0.5).reshape(-1, 1))  # len = 1
    outward_start = np.array([[-1,0,0]])
    outward_end = np.array([[1,0,0]])
    outward_dir = np.concatenate([outward_start, outward_dir, outward_end], 0)
    inward_dir = -outward_dir
    outward_points = line_points + outward_dir * offset/2
    inward_points = line_points + inward_dir * offset/2
    offset_points = np.concatenate([inward_points, outward_points[::-1]], 0)
    
    return offset_points




def move_mesh(mesh, move_vector):
    return om.TriMesh(mesh.points() + move_vector, mesh.face_vertex_indices())


def write_om_mesh(mesh, fname):
    om.write_mesh(f"{fname}.vtk", mesh)


def soft_merge_meshes(meshes,
                      vh_idx_to_sync_list=None):  # meshes je lista sa meshevima, vh_idx_to_sync_list sa lista isog lena ko meshes, svaka sadržava array sa vh_idx koji zelimo syncat
    points = np.empty((0, 3))
    merged_fvi = np.empty((0, 3))

    if vh_idx_to_sync_list is None:
        for mesh in meshes:

            mesh_fvi = mesh.face_vertex_indices()
            if mesh_fvi.size == 0:
                continue
            merged_fvi = np.append(merged_fvi, mesh_fvi + points.shape[0],
                                   axis=0)  # +points.shape[0] je tu da poreda face_vertex_indices sa njihovim indexom u novom arrayu
            points = np.append(points, mesh.points(), axis=0)

        return om.TriMesh(points, merged_fvi)

    else:
        synced_vh_idx = []
        for i in range(len(meshes)):
            mesh = meshes[i]
            mesh_fvi = mesh.face_vertex_indices()
            merged_fvi = np.append(merged_fvi, mesh_fvi + points.shape[0],
                                   axis=0)  # +points.shape[0] je tu da poreda face_vertex_indices sa njihovim indexom u novom arrayu
            synced_vh_idx.append(vh_idx_to_sync_list[i] + points.shape[0])
            points = np.append(points, mesh.points(), axis=0)

        return (om.TriMesh(points, merged_fvi), synced_vh_idx)

#resulting evi will have outward normals
def extrude_outline(outline_points,spacing, direction):
    n_outline_points = outline_points.shape[0]
    n_copies = spacing.shape[0]
    
    a1 = np.expand_dims(outline_points,0)
    extruded_points = np.concatenate([a1,]*n_copies, axis = 0)
    #replace x or y values of points with those of spacing:
    extruded_points[:,:,direction] = np.broadcast_to(spacing.reshape(-1,1), (n_copies,n_outline_points))   
    extruded_points = extruded_points.reshape(-1,3)
    
    #generate extruded evi:
    #face_set1:
    b1 = np.arange(n_outline_points).reshape(-1,1)
    a1 = np.roll(b1, -1, 0)
    c1 = a1 + n_outline_points 
    extruded_fvi1 = np.concatenate([a1,b1,c1],axis = -1)
    
    #face_set2:
    c2 = np.arange(n_outline_points).reshape(-1,1)
    a2 = c2+ n_outline_points 
    b2 = np.roll(a2, -1, 0)
    extruded_fvi2 = np.concatenate([a2,b2,c2],axis = -1)
    
    extruded_fvi = np.expand_dims(np.append(extruded_fvi1, extruded_fvi2, 0), 0)
    extruded_fvi = np.concatenate([extruded_fvi,]*n_copies, 0)
    extruded_fvi = extruded_fvi + (np.arange(n_copies)*n_outline_points).reshape(n_copies,1,1)
    extruded_fvi = extruded_fvi.reshape(-1,3)
    
    return (extruded_points,extruded_fvi)

#returns points in xz plane, and fvi
def get_Tshape_data(shape_prop, t, corr_add = None):
    # /1000 because original is in [mm]      
    if corr_add is None:
        hw = shape_prop.hw   # Web height
        tw = shape_prop.tw     # Web thickness
        bf = shape_prop.bf  # Flange width
        tf = shape_prop.tf   # Flange thickness
    else:
        hw = shape_prop.hw_net(corr_add, t)   # Web height
        tw = shape_prop.tw_net(corr_add)     # Web thickness
        bf = shape_prop.bf_net(corr_add)  # Flange width
        tf = shape_prop.tf_net(corr_add)   # Flange thickness
    
    
    
    face_points = np.array([
        [tw/2,0,0],
        [tw/2,0,t/2+hw],
        [bf/2,0,t/2+hw],
        [bf/2,0,t/2+hw+tf],
        [-bf/2,0,t/2+hw+tf],
        [-bf/2,0,t/2+hw],
        [-tw/2,0,t/2+hw],
        [-tw/2,0,0],
        ])
    
    face_fvi = np.array([
        [0,1,6],
        [6,7,0],
        [4,5,6],
        [4,6,1],
        [4,1,3],
        [3,1,2],
        ])
    return (face_points, face_fvi)
    

def get_HPshape_data(shape_prop, t, corr_add = None):
    # /1000 because original is in [mm]      
    if corr_add is None:
        hw = shape_prop.hw_HP   # Web height
        tw = shape_prop.tw_HP     # Web thickness
    else:
        hw = shape_prop.hw_ekv_net(corr_add)   # Web height
        tw = shape_prop.tw_ekv_net(corr_add)     # Web thickness
        
    
    face_points = np.array([
        [0.5,0,0.7665],
        [2.7309,0,0.8833],
        [3.0732,0,0.9067],
        [3.3231,0,0.9317],
        [3.398,0,0.9566],
        [3.2731,0,0.9766],
        [3.0984,0,0.9933],
        [2.87,0,1],
        [-0.5,0,1],
        [-0.5,0,0],
        [0.5,0,0],
        ])
    
    #scale points 
    face_points = face_points * np.array([tw, 1 ,hw + t/2])

    face_fvi = np.array([
        [0,9,10],
        [8,9,0],
        [7,8,0],
        [7,0,1],
        [6,7,1],
        [1,2,6],
        [5,6,2],
        [2,3,4],
        [5,2,4],
        ])
    
    return (face_points, face_fvi) 

def get_HATshape_data(shape_prop, t, corr_add = None):
    if corr_add is None:
        h = shape_prop.h  
        tw = shape_prop.t  
        bf = shape_prop.bf  
    else:
        h = shape_prop.h_net(corr_add)  
        tw = shape_prop.t_net(corr_add)
        bf = shape_prop.bf_net(corr_add)  
    
    fi = shape_prop.fi * 2*np.pi/360        #turn into radians
    
    
    h_mid = h+t/2+tw/2
    web_direction_vector = np.array([np.cos(fi), 0, np.sin(fi)])
    len_to_0 = h_mid / web_direction_vector[2]
    
    
    point1 = np.array([[-bf/2,0,h_mid]])
    point2 = np.array([[bf/2,0,h_mid]])
    point0 = point1 + -web_direction_vector*len_to_0
    point3 = point2 + np.array([1,1,-1])*web_direction_vector*len_to_0
    
    web_line_points = np.concatenate([point0,point1,point2,point3], axis = 0)
    web_line_evi = np.array([
        [0,1],
        [1,2],
        [2,3],
        ])
    
    
    face_points = offset_line(web_line_points, web_line_evi, offset = tw)
    face_fvi = np.array([
        [0,1,6],
        [6,7,0],
        [1,2,6],
        [2,5,6],
        [2,3,5],
        [3,4,5],
        ])
    

    return (face_points, face_fvi) 


def get_Lshape_data(shape_prop, t, corr_add = None):
    # /1000 because original is in [mm]      
    if corr_add is None:
        hw = shape_prop.hw   # Web height
        tw = shape_prop.tw     # Web thickness
        bf = shape_prop.bf  # Flange width
        tf = shape_prop.tf   # Flange thickness
    else:
        hw = shape_prop.hw_net(corr_add, t)   # Web height
        tw = shape_prop.tw_net(corr_add)     # Web thickness
        bf = shape_prop.bf_net(corr_add)  # Flange width
        tf = shape_prop.tf_net(corr_add)   # Flange thickness
    
    face_points = np.array([
        [tw/2,0,0],
        [tw/2,0,t/2],
        [tw/2,0,t/2+hw],
        [bf-tw/2,0,t/2+hw],
        [bf-tw/2,0,t/2+hw+tf],
        [-tw/2,0,t/2+hw+tf],
        [-tw/2,0,t/2],
        [-tw/2,0,0],
        ])
    
    face_fvi = np.array([
        [7,0,1],
        [1,6,7],
        [6,1,2],
        [2,5,6],
        [2,4,5],
        [2,3,4],
        ])
      

    return (face_points, face_fvi)

def get_FBshape_data(shape_prop, t, corr_add = None):
    # /1000 because original is in [mm]      
    if corr_add is None:
        hw = shape_prop.hw   # Web height
        tw = shape_prop.tw     # Web thickness
    else:
        hw = shape_prop.hw_net(corr_add, t)   # Web height
        tw = shape_prop.tw_net(corr_add)     # Web thickness
        
    
    
    face_points = np.array([
        [-tw/2,0,0],
        [tw/2,0,0],
        [tw/2,0,t/2],
        [tw/2,0,t/2+hw],
        [-tw/2,0,t/2+hw],
        [-tw/2,0,t/2]
        ])
    
    face_fvi = np.array([
        [0,1,2],
        [2,5,0],
        [5,2,4],
        [2,3,4],
        ])
    return (face_points, face_fvi)

def extrude_shape(shape_prop, t, orientation, spacing, position, corr_add = None, L_flange_direction = None):
    if isinstance(shape_prop, TBeamProperty):
        data = get_Tshape_data(shape_prop, t, corr_add) 

    if isinstance(shape_prop, LBeamProperty):
        data = get_Lshape_data(shape_prop, t, corr_add) 

    if isinstance(shape_prop, FBBeamProperty):
        data = get_FBshape_data(shape_prop, t, corr_add) 

    if isinstance(shape_prop, BulbBeamProperty):
        data = get_HPshape_data(shape_prop, t, corr_add) 

    if isinstance(shape_prop, HatBeamProperty):
        data = get_HATshape_data(shape_prop, t, corr_add)
    #more shapes here:
    
    
    face_points = data[0]
    
    start_face_fvi = data[1]
    

    end_face_fvi = start_face_fvi + (spacing.shape[0]-1)*face_points.shape[0]
    end_face_fvi = end_face_fvi[:,::-1]   #reverse direction
    
    if orientation is bd.LONGITUDINAL:
        # rotate points:
        face_points[:,1] = -face_points[:,0]
        face_points[:,0] = 0
        direction_ind = 0
        
        #special rotation for L beam orientation:        
        if np.equal(L_flange_direction, np.array([0,1,0])).all():   #if flange is in y+ direction, y- direction is default for longitudinal beams! 
            face_points[:,1] = -face_points[:,1] 
    
        
        
    elif orientation is bd.TRANSVERSE:
        direction_ind = 1
        
        #special rotation for L beam orientation: 
        if np.equal(L_flange_direction, np.array([-1,0,0])).all():   #if flange is in x- direction, x+ direction is default for transverse beams! 
            face_points[:,0] = -face_points[:,0] 

           
            

    #translate points:    
    face_points = face_points + position    
    data = extrude_outline(face_points, spacing, direction_ind)
    
    extruded_fvi = data[1]
    extruded_points = data[0]
    
    extruded_fvi = np.concatenate([start_face_fvi, extruded_fvi, end_face_fvi])
    mesh = om.TriMesh(extruded_points, extruded_fvi)
    return mesh



class BeamSegmentMesher():
    def __init__(self, beam, element_size):
      
        # self.element_size = 0.6 # [m]    dodaj uvijet ako je element size veci od beam lena
        self.element_size = element_size # [m]    dodaj uvijet ako je element size veci od beam lena
        self.beam = beam
        self.primary_supp_mem = self.beam._primary_supp_mem
        
        try:
            self.corr_add = self.primary_supp_mem.grillage.corrosion_addition()[1]
        except:
            self.corr_add = None
            
        self.beam_prop = self.beam._beam_prop
        if  isinstance(self.beam_prop, LBeamProperty):
            self.L_flange_direction = self.primary_supp_mem.flange_direction_vector 
        else:
            self.L_flange_direction = None
            
        self.beam_dir = self.primary_supp_mem._direction
        self.beam_len = self.beam.segment_len() * 1000
        self.beam_position = self.beam.get_segment_node1()
        self.t = self.beam.get_attplate()[1] 
        self.attatched_plating = self.primary_supp_mem.grillage.segment_common_plates(self.beam)    #list [p1, p2]
        self.make_spacing()      # points at which beam mesh will be 
        
    def make_spacing(self):    
        if self.beam_dir is bd.LONGITUDINAL:    #if beam is in x direction:
            direction_id = 0
        elif self.beam_dir is bd.TRANSVERSE:
            direction_id = 1
        
        n_of_spaces = self.beam_len/self.element_size
        if n_of_spaces.is_integer() == False:   #if beam cannot be divided into segments of desiered size, determine leftover len and put half of it ant start and end of beam 
            extra_space = np.array([(self.beam_len - np.floor(n_of_spaces)*self.element_size)/2])
            spacing = np.array([self.element_size]*int(np.floor(n_of_spaces)))
            spacing = np.concatenate([np.array([0]), extra_space, spacing, extra_space])
        else:
            spacing = np.array([self.element_size]*int(np.floor(n_of_spaces)))
            spacing = np.concatenate([np.array([0]), spacing])
        
        spacing = np.broadcast_to(spacing, (spacing.shape[0], spacing.shape[0]))
        
        spacing = np.tril(spacing).sum(-1)
        spacing_points = np.expand_dims(self.beam_position,0)
        spacing_points = np.concatenate([spacing_points]*spacing.shape[0], 0)
        spacing = spacing_points[:, direction_id] + spacing
        
        
        #add spacing points for attatched plate stiffeners:
        additional_stiff_spacing = np.empty((0))
        for plate in self.attatched_plating:
            try:    #if attatched plate has no stiff direction, skip plate
                plate_stiff_dir = plate.stiff_dir
            except:
                continue
            
            if plate_stiff_dir is not self.beam_dir:    #if beam segment and attatched stiff orientation is different, beam segment must be divided with stiffener start points
                stiff_N = int(plate.get_stiffener_number())
                stiff_start_points = []
                for i in range(stiff_N):
                    stiff_start_point = plate.get_stiff_coords(i+1)[0]
                    stiff_start_points.append(stiff_start_point)
                stiff_start_points = np.asarray(stiff_start_points)
                
                if plate_stiff_dir is bd.LONGITUDINAL:    
                    stiff_spacing = stiff_start_points[:,1]
                
                elif plate_stiff_dir is bd.TRANSVERSE:
                    stiff_spacing = stiff_start_points[:,0]
                
                additional_stiff_spacing = np.append(additional_stiff_spacing, stiff_spacing)
        
        if additional_stiff_spacing.shape[0] > 0:   # if spacing was added, append new spacing to beam spacing and sort them
            spacing = np.append(spacing, additional_stiff_spacing)
            spacing = np.sort(spacing)
            spacing = np.unique(spacing)        #remove duplicates incase of coincident nodes to prevent errors
            
        self.spacing = spacing    
        

        
    def make_beam_segment_mesh(self):
        
        beam_mesh = extrude_shape(self.beam_prop, self.t, self.beam_dir, self.spacing, self.beam_position, self.corr_add, self.L_flange_direction)
        
        return beam_mesh
      
class PlateMesher():
    def __init__(self, plate, element_size):
        # self.element_size = 0.6
        self.element_size = element_size
        self.plate = plate

        try:
            self.corr_add = self.plate.long_seg1.primary_supp_mem.grillage.corrosion_addition()[1]
            self.t = self.plate.plate_prop.tp_net(self.corr_add, self.plate.plate_prop.tp)
        except:
            self.corr_add = None
            self.t = self.plate.plate_prop.tp 
            
        #plate and stiff properties:
        self.stiff_N = int(self.plate.get_stiffener_number())
        self.t = self.plate.plate_prop.tp                                            #plate thickness [mm] so /1000
        self.xlen = self.plate.plate_longitudinal_dim()*1000                                #x plate lenght    [m]
        self.ylen = self.plate.plate_transverse_dim()*1000                                 #y plate lenght    [m]        
        self.plate_pos = self.plate._segments[0].get_segment_node1()
       
        self.stiff_orientation = self.plate.stiff_dir 
        self.stiff_property = self.plate.stiff_layout._beam_prop
        self.stiff_type = type(self.stiff_property)                    #T, HP, HAT
        self.make_spacing()
        
        
    def make_spacing(self):    
        #x spacing :
        n_of_x_spaces = self.xlen/self.element_size
        if n_of_x_spaces.is_integer() == False:   #if beam cannot be divided into segments of desiered size, determine leftover len and put half of it ant start and end of beam 
            x_extra_space = np.array([(self.xlen - np.floor(n_of_x_spaces)*self.element_size)/2])
            x_spacing = np.array([self.element_size]*int(np.floor(n_of_x_spaces)))
            x_spacing = np.concatenate([np.array([0]), x_extra_space, x_spacing, x_extra_space])
        else:
            x_spacing = np.array([self.element_size]*int(np.floor(n_of_x_spaces)))
            x_spacing = np.concatenate([np.array([0]), x_spacing])
        
        x_spacing = np.broadcast_to(x_spacing, (x_spacing.shape[0], x_spacing.shape[0]))
        
        x_spacing = np.tril(x_spacing).sum(-1)
        x_spacing_points = np.expand_dims(self.plate_pos, 0)
        x_spacing_points = np.concatenate([x_spacing_points]*x_spacing.shape[0], 0)
        x_spacing = x_spacing_points[:, 0] + x_spacing
        
        
        #y spacing
        n_of_y_spaces = self.ylen/self.element_size
        if n_of_y_spaces.is_integer() == False:   #if beam cannot be divided into segments of desiered size, determine leftover len and put half of it ant start and end of beam 
            y_extra_space = np.array([(self.ylen - np.floor(n_of_y_spaces)*self.element_size)/2])
            y_spacing = np.array([self.element_size]*int(np.floor(n_of_y_spaces)))
            y_spacing = np.concatenate([np.array([0]), y_extra_space, y_spacing, y_extra_space])
        else:
            y_spacing = np.array([self.element_size]*int(np.floor(n_of_y_spaces)))
            y_spacing = np.concatenate([np.array([0]), y_spacing])
        
        y_spacing = np.broadcast_to(y_spacing, (y_spacing.shape[0], y_spacing.shape[0]))
        
        y_spacing = np.tril(y_spacing).sum(-1)
        y_spacing_points = np.expand_dims(self.plate_pos, 0)
        y_spacing_points = np.concatenate([y_spacing_points]*y_spacing.shape[0], 0)
        y_spacing = y_spacing_points[:, 1] + y_spacing
        
        #get stiffner points:
        stiff_start_points = []
        for i in range(self.stiff_N):
            stiff_start_point = self.plate.get_stiff_coords(i+1)[0]
            stiff_start_points.append(stiff_start_point)
        stiff_start_points = np.asarray(stiff_start_points)
        
        if self.stiff_orientation is bd.LONGITUDINAL:    #if beam is in x direction:
            stiff_spacing = stiff_start_points[:,1]
            y_spacing = np.append(y_spacing, stiff_spacing)	
            y_spacing = np.sort(y_spacing)
            y_spacing = np.unique(y_spacing)        #remove duplicates incase of coincident nodes to prevent errors
            
        elif self.stiff_orientation is bd.TRANSVERSE:
            stiff_spacing = stiff_start_points[:,0]
            x_spacing = np.append(x_spacing, stiff_spacing)
            x_spacing = np.sort(x_spacing)
            x_spacing = np.unique(x_spacing)
        
        self.x_spacing = x_spacing
        self.y_spacing = y_spacing
        self.stiff_start_points = stiff_start_points
    
        
    def make_plate_mesh(self):
        xn = self.x_spacing.shape[0] - 1 #not really -1 but everywhere it is used in code it needs to be -1
        yn  = self.y_spacing.shape[0]
        plate_points = np.asarray(np.meshgrid(self.x_spacing, self.y_spacing, np.array([-self.t/2,self.t/2]))).T.reshape(-1,3) 
        #translate points up: 
        plate_points[:,2] = plate_points[:,2] + self.plate_pos[2]
        #make top fvi:
        b1 = np.arange(yn-1).reshape(-1,1)
        a1 = b1 + 1
        c1 = a1 + yn
        plate_fvi1 = np.concatenate([a1,b1,c1],axis = -1)
        
        c2 = np.arange(yn-1).reshape(-1,1)
        a2 = c2 + yn
        b2 = a2 + 1
        plate_fvi2 = np.concatenate([a2,b2,c2],axis = -1)
        
        plate_top_fvi = np.expand_dims(np.append(plate_fvi1, plate_fvi2, 0),0) 
        plate_top_fvi = np.concatenate([plate_top_fvi,]*xn, 0)
        plate_top_fvi = plate_top_fvi + (np.arange(xn)*yn).reshape(xn,1,1)
        plate_top_fvi = plate_top_fvi.reshape(-1,3)
        
        # make bot fvi
        plate_bot_fvi = plate_top_fvi + ((xn+1)*yn) 
        plate_bot_fvi = plate_bot_fvi[:,::-1]   #reverse bot plate orientation
        plate_fvi = np.append(plate_top_fvi, plate_bot_fvi,0)
        
        
        #make edge fvi
        #top edge indices:
        e1t = np.arange(yn)
        e2t = e1t + xn*yn
        e3t = np.arange(xn+1)*yn
        e4t = e3t + (yn - 1)
        
        #bot edge indices:
        e1b = e1t + ((xn+1)*yn)
        e2b = e2t + ((xn+1)*yn)
        e3b = e3t + ((xn+1)*yn)
        e4b = e4t + ((xn+1)*yn)
        
        e1 = self.stitch_plate_edge(e1t, e1b)
        e2 = self.stitch_plate_edge(e2t, e2b)[:,::-1]    #reverse fvi orientation
        e3 = self.stitch_plate_edge(e3t, e3b)[:,::-1]
        e4 = self.stitch_plate_edge(e4t, e4b)
        
        edge_fvi = np.concatenate([e1, e2, e3, e4], 0)
        plate_fvi = np.append(plate_fvi, edge_fvi, 0)
        
        mesh = om.TriMesh(plate_points, plate_fvi)

        
        return mesh
        
    def stitch_plate_edge(self, et,eb):   #input: edge indices top, edge indices bot
        etn = et.shape[0]
        ei = np.arange(etn-1)        #edge indices
        
        #stitch1
        a1 = et[ei].reshape(-1,1)
        b1 = et[ei + 1].reshape(-1,1)
        c1 = eb[ei].reshape(-1,1)
        fvi1 = np.concatenate([a1,b1,c1], -1)
        
        #stitch2
        a2 = eb[ei+1].reshape(-1,1)
        b2 = eb[ei].reshape(-1,1)
        c2 = et[ei + 1].reshape(-1,1)
        fvi2 = np.concatenate([a2,b2,c2], -1)
        
        fvi = np.append(fvi1,fvi2,0)
        return fvi   



    def make_stiffener_mesh(self):
        stiff_meshes = []
        start_pos = self.stiff_start_points[0]
        other_pos = list(self.stiff_start_points[1:])
        if self.stiff_orientation is bd.LONGITUDINAL:
            spacing = self.x_spacing
        elif self.stiff_orientation is bd.TRANSVERSE:
            spacing = self.y_spacing
            
            
        stiff_start_mesh = extrude_shape(self.stiff_property, self.t, self.stiff_orientation, spacing, start_pos)
        stiff_meshes.append(stiff_start_mesh)
        
        for pos in other_pos:
            move_delta = pos - start_pos
            mesh = move_mesh(copy.copy(stiff_start_mesh), move_delta)
            stiff_meshes.append(mesh)
        
        return stiff_meshes
    
      
    
    
    

class GrillageBaseGeometry(Geometry):
    def __init__(self, name = '', element_size = 500):
        self.element_size = element_size        
        super().__init__(name)
        self._gen_mesh_or_subgeometry()
        pass

    def regenerateMesh(self):
        self._gen_mesh_or_subgeometry()

    def _gen_mesh_or_subgeometry(self):
        pass
    def _appendmesh(self,mesh:om.TriMesh):
        new_mesh = join_tria_meshes(mesh,self.mesh)
        for sg in self.sub_geometry:
            new_mesh = sg._appendmesh(new_mesh)
        return new_mesh

class PlateGeometry(GrillageBaseGeometry):
    def __init__(self, plate: Plate, name = '', element_size = 500):
        self._plate = plate
        super().__init__(name, element_size)

    def _gen_mesh_or_subgeometry(self):
        self.sub_geometry.clear()
        mesher = PlateMesher(self._plate, self.element_size)
        meshs = mesher.make_plate_mesh()
        self.mesh = meshs

class StiffenersGeometry(GrillageBaseGeometry):
    def __init__(self, plate: Plate, name = '', element_size = 500):
        self._plate = plate
        super().__init__(name,element_size)

    def _gen_mesh_or_subgeometry(self):
        self.sub_geometry.clear()
        mesher = PlateMesher(self._plate, self.element_size)
        meshs = mesher.make_stiffener_mesh()
        new_mesh = om.TriMesh()
        for mesh in meshs:
            new_mesh = join_tria_meshes(new_mesh,mesh)
        self.mesh = new_mesh


class StiffenedPlateGeometry(GrillageBaseGeometry):
    def __init__(self, plate: Plate, name = '', element_size = 500):
        self._plate = plate
        super().__init__(name,element_size)

    def _gen_mesh_or_subgeometry(self):
        self.mesh = om.TriMesh()
        self.sub_geometry.clear()
        self.sub_geometry.append(PlateGeometry(self._plate, 'Plate',
                                                           self.element_size))
        self.sub_geometry.append(StiffenersGeometry(self._plate, 'Stiffeners',
                                                           self.element_size))


class PlatingGeometry(GrillageBaseGeometry):
    def __init__(self, plates: Plate, name = '', element_size = 500):
        self._plates = plates
        super().__init__(name,element_size)

    def _gen_mesh_or_subgeometry(self):
        self.mesh = om.TriMesh()
        self.sub_geometry.clear()
        for id,plate in self._plates.items():
            self.sub_geometry.append(StiffenedPlateGeometry(plate, 'Stiff_plate_'+str(id),
                                                           self.element_size))

class BeamSegmentGeometry(GrillageBaseGeometry):
    def __init__(self, segment: Segment, name = '', element_size = 500):
        self._segment = segment
        super().__init__(name,element_size)

    def _gen_mesh_or_subgeometry(self):
        self.sub_geometry.clear()
        mesher = BeamSegmentMesher(self._segment, self.element_size)
        self.mesh = mesher.make_beam_segment_mesh()

class PrimarySuppMemGeometry(GrillageBaseGeometry):
    def __init__(self, psm: PrimarySuppMem, name = '', element_size = 500):
        self._psm:PrimarySuppMem = psm
        super().__init__(name, element_size)

    def _gen_mesh_or_subgeometry(self):
        self.mesh = om.TriMesh()
        self.sub_geometry.clear()
        for beam_segment in self._psm._segments:
            self.sub_geometry.append(BeamSegmentGeometry(beam_segment, 'Segment_'+str(beam_segment.id),
                                                           self.element_size))


class Longs_or_Trans_Geometry(GrillageBaseGeometry):
    def __init__(self, items: Dict[int,PrimarySuppMem], name = '', element_size = 500):
        self._items: Dict[int, PrimarySuppMem] = items
        super().__init__(name,element_size)

    def _gen_mesh_or_subgeometry(self):
        self.mesh = om.TriMesh()
        self.sub_geometry.clear()
        for id, psm in self._items.items():
            self.sub_geometry.append(PrimarySuppMemGeometry(psm,'Beam_'+str(id),
                                                           self.element_size))

class LongitudinalsGeometry(Longs_or_Trans_Geometry):

    def __init__(self, items: Dict[int,PrimarySuppMem], name = '', element_size = 500):
        super().__init__(items,name,element_size)

class TransversalsGeometry(Longs_or_Trans_Geometry):

    def __init__(self, items: Dict[int,PrimarySuppMem], name = '', element_size = 500):
        super().__init__(items,name,element_size)


class GrillageGeometry(GrillageBaseGeometry):

    def __init__(self, grillage: Grillage, name = '',element_size = 600):
        self._grill = grillage
        super().__init__(name,element_size)

    @property
    def grillage(self):
        return self._grill

    def _gen_mesh_or_subgeometry(self):
        try:
            self.mesh=om.TriMesh()
            self.sub_geometry.clear()
            self.sub_geometry.append(LongitudinalsGeometry(self._grill.longitudinal_members(),'Longitudinal Beams',
                                                            self.element_size))
            self.sub_geometry.append(TransversalsGeometry(self._grill.transverse_members(),'Transverse Beams',
                                                            self.element_size))
            self.sub_geometry.append(PlatingGeometry(self._grill.plating(), 'Plating',
                                                            self.element_size))
            self.mesh = self._appendmesh(self.mesh)
            
        except BaseException as error:
            print('An exception occurred during visualization mesh generation: {}'.format(error))
        except:
            print('Unknown exception occurred during visualization mesh generation')




















































