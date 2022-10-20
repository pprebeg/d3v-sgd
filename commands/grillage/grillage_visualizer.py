import numpy as np
import openmesh as om
import pygmsh as pg 
import copy
#import grillage_model as gm
from grillage.grillage_model import BeamDirection as bd
from grillage.grillage_model import Plate,PrimarySuppMem,Grillage,Segment
from grillage.grillage_model import TBeamProperty,BulbBeamProperty,HatBeamProperty
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
	# boundary_vi = [0,5]
	normal = np.array([0, 1.0, 0])
	epoints = line_points[line_evi]
	ep_dir_vectors = epoints[:, 1] - epoints[:, 0]
	outward_dir = np.cross(ep_dir_vectors, normal)
	outward_dir = outward_dir / (((outward_dir ** 2).sum(-1) ** 0.5).reshape(-1, 1))  # len = 1
	offset_vectors_middle = outward_dir[1:] + outward_dir[:-1]
	# offset_vectors_middle = offset_vectors_middle/(((offset_vectors_middle**2).sum(-1)**0.5).reshape(-1,1))
	offset_points_middle_top = line_points[1:-1] + offset_vectors_middle * offset
	offset_points_middle_bot = line_points[1:-1] - offset_vectors_middle * offset
	# boundary_offset_point11 = line_points[0:1] + np.array([0,0,1.0]) * offset
	# boundary_offset_point12 = line_points[0:1] - np.array([0,0,1.0]) * offset
	# boundary_offset_point21 = line_points[5:6] + np.array([0,0,1.0]) * offset
	# boundary_offset_point22 = line_points[5:6] - np.array([0,0,1.0]) * offset
	boundary_offset_point11 = copy.copy(line_points[0:1])
	boundary_offset_point12 = copy.copy(line_points[0:1])
	boundary_offset_point21 = copy.copy(line_points[5:6])
	boundary_offset_point22 = copy.copy(line_points[5:6])
	boundary_offset_point11[:, 2] = offset_points_middle_top[0, 2]
	boundary_offset_point12[:, 2] = offset_points_middle_bot[0, 2]
	boundary_offset_point21[:, 2] = offset_points_middle_top[-1, 2]
	boundary_offset_point22[:, 2] = offset_points_middle_bot[-1, 2]

	offset_points_bot = np.concatenate((boundary_offset_point12, offset_points_middle_bot, boundary_offset_point22),
									   axis=0)
	offset_points_top = np.concatenate(
		(boundary_offset_point21, offset_points_middle_top[::-1], boundary_offset_point11), axis=0)

	offset_points = np.concatenate((offset_points_bot, offset_points_top), axis=0)
	offset_points = offset_points + np.array([0, 0, np.abs(offset_points_middle_bot[-1, 2])])  # profile centering

	return offset_points

def pg_to_om_mesh(mesh):
	pg_mesh_points = np.array(mesh.points)
	pg_mesh_fvi = np.array(mesh.cells_dict["triangle"])
	return om.TriMesh(pg_mesh_points, pg_mesh_fvi)


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
			# print(mesh)
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

class BeamSegmentMesher():
	def __init__(self, beam, mesh_resolution = 0.1, plate_n_of_subdivision = 0):
		self.mesh_resolution = mesh_resolution
		self.plate_n_of_subdivision = plate_n_of_subdivision
		self.beam_prop = beam._beam_prop
		self.primary_supp_mem = beam._primary_supp_mem
		self.beam_dir = self.primary_supp_mem._direction
		self.beam_len = beam.segment_len()
		self.beam_position = beam.get_segment_node1()
		# print(f"beam pos: {self.beam_position}")
		# self.beam_node2 = beam.get_segment_node2()
		# self.beam_vector = self.beam_node2 - self.beam_node1
		self.t = beam.get_attplate()[1] / 1000
		
	

	def make_beam_segment_mesh(self):
		# /1000 because original is in [mm]	  
		hw = self.beam_prop.hw / 1000  # Web height
		tw = self.beam_prop.tw / 1000	# Web thickness
		bf = self.beam_prop.bf / 1000 # Flange width
		tf = self.beam_prop.tf / 1000  # Flange thickness
		
		beam_points_transverse = np.array([
			[tw/2,0,0],
			[tw/2,0,self.t/2+hw],
			[bf/2,0,self.t/2+hw],
			[bf/2,0,self.t/2+hw+tf],
			[-bf/2,0,self.t/2+hw+tf],
			[-bf/2,0,self.t/2+hw],
			[-tw/2,0,self.t/2+hw],
			[-tw/2,0,0],
			])
		
		
		if self.beam_dir == bd.LONGITUDINAL:
			extrude_vector = (self.beam_len,0.0,0.0)
			# rotate points:
			beam_points_longitudinal = copy.copy(beam_points_transverse)
			beam_points_longitudinal[:,1] = -beam_points_longitudinal[:,0]
			beam_points_longitudinal[:,0] = 0
			beam_points = beam_points_longitudinal

		elif self.beam_dir == bd.TRANSVERSE:
			extrude_vector = (0.0,self.beam_len,0.0)
			beam_points = beam_points_transverse
			
			
		beam_lines = [
			[0,1],
			[1,2],
			[2,3],
			[3,4],
			[4,5],
			[5,6],
			[6,7],
			[7,0],
			]
		
		with pg.occ.geometry.Geometry() as geom:
			pg_points = []
			for point in list(beam_points):
				pg_points.append(geom.add_point(point, self.mesh_resolution))			#point size is 5000 to ensure that only triangle points are used for face generation
			
			pg_lines = []
			for line in beam_lines: 	
				pg_lines.append(geom.add_line(pg_points[line[0]], pg_points[line[1]])),
			
			pg_loop = geom.add_curve_loop([*pg_lines])
			
			pg_surface = geom.add_plane_surface(pg_loop)
			

			volume = geom.extrude(pg_surface, extrude_vector)

			mesh = geom.generate_mesh(algorithm = 9)
			mesh = pg_to_om_mesh(mesh)
			
			mesh = move_mesh(mesh, self.beam_position)
			
			return mesh

	
class PlateMesher():
	def __init__(self, plate, mesh_resolution = 0.1, plate_n_of_subdivision = 0):
		self.mesh_resolution =  mesh_resolution		#DO NOT GO BELOW 0.05 UNLESS YOU HAVE STRONG CPU
		self.plate_n_of_subdivision = plate_n_of_subdivision		# 0 for no subdivisions
		self.plate = plate
		
		# empty lists containing meshes 
		self.plate_segments_meshes = []		#starts from 0,0 fills x first
		self.stiff_segments_meshes = []
		
		
		#plate and stiff properties:
		self.stiff_N = int(self.plate.get_stiffener_number())
		self.t = self.plate.plate_prop.tp/1000											#plate thickness [mm] so /1000
		self.xlen = self.plate.plate_longitudinal_dim()								#x plate lenght	[m]
		self.ylen = self.plate.plate_transverse_dim()								#y plate lenght	[m]		
		# self.plate_pos = self.plate.get_reference_segment().get_segment_node1()
		
		# if self.plate._ref_edge == gm.Ref.EDGE1:
		self.plate_pos = self.plate._segments[0].get_segment_node1()
		# elif self.plate._ref_edge == gm.Ref.EDGE2:
		# self.plate_pos = self.plate._segments[2]
		
		
		# self.plate_pos = self.plate.get_reference_segment().get_segment_node2()
		# print(f"plate post: {self.plate_pos}")
		# print(f"stiff cords: {self.plate.get_stiff_coords(self.stiff_N)}")
		
		self.stiff_orientation = self.plate.stiff_dir 
		self.stiff_property = self.plate.stiff_layout._beam_prop
		self.stiff_type = type(self.stiff_property)					#T, HP, HAT
		self.stiff_spacing = self.plate.get_stiffener_spacing()
		
		# min n of plate segments is equal to n_stiff
		self.n_of_plate_segments_in_row = int(((self.stiff_N + 1) * (self.plate_n_of_subdivision + 1)))
		self.x_plate_segment_len = self.xlen / self.n_of_plate_segments_in_row
		self.y_plate_segment_len = self.ylen / self.n_of_plate_segments_in_row
		# print(self.ylen, self.y_plate_segment_len, self.stiff_N, self.stiff_spacing)
		# print(self.xlen, self.plate_pos)
		
		
		# return soft_merge_meshes(self.stiff_segments_meshes+self.plate_segments_meshes )
		# plate_mesh = soft_merge_meshes(self.stiff_segments_meshes+self.plate_segments_meshes )
		# self.write_om_mesh(plate_mesh, "segments_meshes!!!!!!!!!!!!!")

	

	def make_plate_mesh(self):		
		
		plate_segment_mesh_template = self.make_plate_segment_mesh()
		
		# takes segment template, copies it and moves it to correct position to make plate
		for i in range(self.n_of_plate_segments_in_row):
			y_step = i*self.y_plate_segment_len
			for j in range(self.n_of_plate_segments_in_row):
				x_step = j*self.x_plate_segment_len
				move_vector = np.array([x_step, y_step, -self.t/2]) 
				plate_segment_mesh = copy.copy(plate_segment_mesh_template)
				plate_segment_mesh = move_mesh(plate_segment_mesh, move_vector)
				self.plate_segments_meshes.append(plate_segment_mesh)
		
		#return self.plate_segments_meshes
		plate_mesh = soft_merge_meshes(self.plate_segments_meshes)
		# self.write_om_mesh(plate_mesh, "plate")
		return plate_mesh

	def make_plate_segment_mesh(self):
		plate_dims = np.array([self.x_plate_segment_len, self.y_plate_segment_len, self.t])
		with pg.occ.geometry.Geometry() as geom:
			geom.add_box(self.plate_pos, plate_dims, mesh_size = self.mesh_resolution)
			mesh = geom.generate_mesh(algorithm = 9)		#9    ,6,8
			mesh = pg_to_om_mesh(mesh)
			return mesh

	def make_stiffener_mesh(self):
			
		if self.stiff_type is TBeamProperty:
			stiff_segment_mesh_template = self.make_Tstiff_segment_mesh()
		if self.stiff_type is BulbBeamProperty:
			stiff_segment_mesh_template = self.make_HPstiff_mesh()
		if self.stiff_type is HatBeamProperty:
			stiff_segment_mesh_template = self.make_HATstiff_mesh()
			
		if self.stiff_orientation == bd.LONGITUDINAL:
			for i in range(self.stiff_N):
				y_step = (i+1)*self.stiff_spacing
				for j in range(self.n_of_plate_segments_in_row):
					x_step = j*self.x_plate_segment_len
					move_vector = np.array([x_step, y_step, 0]) + self.plate_pos
					stiff_segment_mesh = copy.copy(stiff_segment_mesh_template)
					stiff_segment_mesh = move_mesh(stiff_segment_mesh, move_vector)
					self.stiff_segments_meshes.append(stiff_segment_mesh)
		
		if self.stiff_orientation == bd.TRANSVERSE:
			for i in range(self.stiff_N):
				x_step = (i+1)*self.stiff_spacing
				for j in range(self.n_of_plate_segments_in_row):
					y_step = j*self.y_plate_segment_len
					move_vector = np.array([x_step, y_step, 0]) + self.plate_pos
					stiff_segment_mesh = copy.copy(stiff_segment_mesh_template)
					stiff_segment_mesh = move_mesh(stiff_segment_mesh, move_vector)
					self.stiff_segments_meshes.append(stiff_segment_mesh)
		
		return self.stiff_segments_meshes
	
	
	def make_Tstiff_segment_mesh(self):
		# /1000 because original is in [mm]	  
		hw = self.stiff_property.hw / 1000  # Web height
		tw = self.stiff_property.tw / 1000	# Web thickness
		bf = self.stiff_property.bf / 1000 # Flange width
		tf = self.stiff_property.tf / 1000  # Flange thickness
		
		Tstiff_points_transverse = np.array([
			[tw/2,0,0],
			[tw/2,0,self.t/2+hw],
			[bf/2,0,self.t/2+hw],
			[bf/2,0,self.t/2+hw+tf],
			[-bf/2,0,self.t/2+hw+tf],
			[-bf/2,0,self.t/2+hw],
			[-tw/2,0,self.t/2+hw],
			[-tw/2,0,0],
			])
		
		
		if self.stiff_orientation == bd.LONGITUDINAL:
			extrude_vector = (self.x_plate_segment_len,0.0,0.0)
			# rotate points:
			Tstiff_points_longitudinal = copy.copy(Tstiff_points_transverse)
			Tstiff_points_longitudinal[:,1] = -Tstiff_points_longitudinal[:,0]
			Tstiff_points_longitudinal[:,0] = 0
			Tstiff_points = Tstiff_points_longitudinal

		elif self.stiff_orientation == bd.TRANSVERSE:
			extrude_vector = (0.0,self.y_plate_segment_len,0.0)
			Tstiff_points = Tstiff_points_transverse
			
			
		Tstiff_lines = [
			[0,1],
			[1,2],
			[2,3],
			[3,4],
			[4,5],
			[5,6],
			[6,7],
			[7,0],
			]
		
		with pg.occ.geometry.Geometry() as geom:
			pg_points = []
			for point in list(Tstiff_points):
				pg_points.append(geom.add_point(point, self.mesh_resolution))			
			
			pg_lines = []
			for line in Tstiff_lines: 	
				pg_lines.append(geom.add_line(pg_points[line[0]], pg_points[line[1]])),
			
			pg_loop = geom.add_curve_loop([*pg_lines])
			
			pg_surface = geom.add_plane_surface(pg_loop)
			

			volume = geom.extrude(pg_surface, extrude_vector)

			mesh = geom.generate_mesh(algorithm = 9)
			mesh = pg_to_om_mesh(mesh)
			return mesh
		
		
	def make_HPstiff_mesh(self):
		# /1000 because original is in [mm]	  
		hw = self.stiff_property.hw_HP / 1000  # Web height
		tw = self.stiff_property.tw_HP / 1000	# Web thickness
			
		HPstiff_points_transverse = np.array([
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
		
		#scale and move points 
		HPstiff_points_transverse = HPstiff_points_transverse * np.array([tw, 1 ,hw + self.t/2])
		
		
		if self.stiff_orientation == bd.LONGITUDINAL:
			extrude_vector = (self.x_plate_segment_len,0.0,0.0)
			#rotate points:
			HPstiff_points_longitudinal = copy.copy(HPstiff_points_transverse)
			HPstiff_points_longitudinal[:,1] = -HPstiff_points_longitudinal[:,0]
			HPstiff_points_longitudinal[:,0] = 0
			HPstiff_points = HPstiff_points_longitudinal
		
		
		
		elif self.stiff_orientation == bd.TRANSVERSE:
			extrude_vector = (0.0,self.y_plate_segment_len,0.0)
			HPstiff_points = HPstiff_points_transverse
			
			
		HPstiff_lines = [
			[0,1],
			[1,2],
			[2,3],
			[3,4],
			[4,5],
			[5,6],
			[6,7],
			[7,8],
			[8,9],
			[9,10],
			[10,0],
			]
		
		with pg.occ.geometry.Geometry() as geom:
			pg_points = []
			for point in list(HPstiff_points):
				pg_points.append(geom.add_point(point, self.mesh_resolution))			
			
			pg_lines = []
			for line in HPstiff_lines: 	
				pg_lines.append(geom.add_line(pg_points[line[0]], pg_points[line[1]])),
			
			pg_loop = geom.add_curve_loop([*pg_lines])
			
			pg_surface = geom.add_plane_surface(pg_loop)
			

			volume = geom.extrude(pg_surface, extrude_vector)

			mesh = geom.generate_mesh(algorithm = 9)
			mesh = pg_to_om_mesh(mesh)
			return mesh
		
	
	def make_HATstiff_mesh(self):
		h = self.stiff_property.h / 1000 
		tw = self.stiff_property.t / 1000 
		t = self.t
		bf = self.stiff_property.bf / 1000 
		fi = self.stiff_property.fi * 2*np.pi/360		#turn into radians
		
		
		
		web_line_points = np.array([
			[-bf-(h+tw/2+t/4)*np.sin(fi),0,tw/2+t/4],
			[-bf/2-(h+tw/2+t/4)*np.sin(fi),0,tw/2+t/4],
			[-bf/2,0,1.5*tw+t/2+h],
			[bf/2,0,1.5*tw+t/2+h],
			[bf/2+(h+tw/2+t/4)*np.sin(fi),0,tw/2+t/4],
			[bf+(h+tw/2+t/4)*np.sin(fi),0,tw/2+t/4],
			]) 
		
		web_line_evi = np.array([
			[0,1],
			[1,2],
			[2,3],
			[3,4],
			[4,5],
			])
		
		HATstiff_points_transverse = offset_line(web_line_points, web_line_evi, offset = tw)
		# HATstiff_points_transverse = HATstiff_points_transverse[::-1]
		HATstiff_lines = np.array([
			[0,1],
			[1,2],
			[2,3],
			[3,4],
			[4,5],
			[5,6],
			[6,7],
			[7,8],
			[8,9],
			[9,10],
			[10,11],
			[11,0],
			])
		
		if self.stiff_orientation == bd.LONGITUDINAL:
			extrude_vector = (self.x_plate_segment_len,0.0,0.0)
			#rotate points:
			HATstiff_points_longitudinal = copy.copy(HATstiff_points_transverse)
			HATstiff_points_longitudinal[:,1] = -HATstiff_points_longitudinal[:,0]
			HATstiff_points_longitudinal[:,0] = 0
			HATstiff_points = HATstiff_points_longitudinal
		
		
		
		elif self.stiff_orientation == bd.TRANSVERSE:
			extrude_vector = (0.0,self.y_plate_segment_len,0.0)
			HATstiff_points = HATstiff_points_transverse
		
		
		with pg.occ.geometry.Geometry() as geom:
			pg_points = []
			for point in list(HATstiff_points):
				pg_points.append(geom.add_point(point, self.mesh_resolution))			
				
			pg_lines = []
			for line in HATstiff_lines: 	
				pg_lines.append(geom.add_line(pg_points[line[0]], pg_points[line[1]])),
			
			pg_loop = geom.add_curve_loop([*pg_lines])
			
			pg_surface = geom.add_plane_surface(pg_loop)
		

			volume = geom.extrude(pg_surface, extrude_vector)

			mesh = geom.generate_mesh(algorithm = 9)

			mesh = pg_to_om_mesh(mesh)
			return mesh
	










# print(hc_var_1.plating[1].stiff_dir)
# print(hc_var_1.plating[1].stiff_layout)
# print(hc_var_1.plating[1].get_stiffener_spacing())

		

			
			# if stiff_type is gm.StiffPropertyHP:
				# stiff_wt = stiff_prop.tw_HP / 1000				#/1000 because original is in [mm]
			# if stiff_type is gm.TBeamProperty:
				# stiff_wt = stiff_prop.tw / 1000					#/1000 because original is in [mm]
				
		
			# if stiff_dir == gm.BeamOrientation.LONGITUDINAL:
				# longitudinal_spacing = stiff_spacing
				# longitudinal_orientation_trans_plate_dims = plate.plate_transverse_dim()
				# plate_dims = np.array([])		# [x_dim, y_dim, thickness]
						
			# elif stiff_dir == gm.BeamOrientation.TRANSVERSE:
				# transverse_spacing = stiff_spacing
				# transverse_orientation_long_plate_dims = plate.plate_longitudinal_dim()
				
			# else:
				# raise Invalid_Stiff_dir("Invalid stiff dir")
		
		
			
		
		# self.tp = self.plating[plate_id].plate_prop.tp/1000											#plate thickness [mm] so /1000
		# self.xlen = self.plating[plate_id].plate_longitudinal_dim()								#x plate lenght	[m]
		# self.ylen = self.plating[plate_id].plate_transverse_dim()								#y plate lenght	[m]		
		# self._B_overall                # Overall width, m
		# self._N_longitudinal = N_longitudinal       # Number of longitudinal primary supporting members
		# print(self.plating[plate_id].plate_longitudinal_dim(),self._B_overall,self._L_overall)
	
		# mesh1 = self.make_plate_mesh()


class GrillageBaseGeometry(Geometry):
	def __init__(self, name = '', mesh_resolution = 0.3, plate_n_of_subdivision = 0):
		self.mesh_resolution =  mesh_resolution		#DO NOT GO BELOW 0.05 UNLESS YOU HAVE STRONG CPU
		self.plate_n_of_subdivision = plate_n_of_subdivision		# 0 for no subdivision
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
	def __init__(self, plate: Plate, name = '', mesh_resolution = 0.1, plate_n_of_subdivision = 0):
		self._plate = plate
		super().__init__(name,mesh_resolution,plate_n_of_subdivision)

	def _gen_mesh_or_subgeometry(self):
		self.sub_geometry.clear()
		mesher = PlateMesher(self._plate, self.mesh_resolution, self.plate_n_of_subdivision)
		meshs = mesher.make_plate_mesh()
		self.mesh = meshs

class StiffenersGeometry(GrillageBaseGeometry):
	def __init__(self, plate: Plate, name = '', mesh_resolution = 0.1, plate_n_of_subdivision = 0):
		self._plate = plate
		super().__init__(name,mesh_resolution,plate_n_of_subdivision)

	def _gen_mesh_or_subgeometry(self):
		self.sub_geometry.clear()
		mesher = PlateMesher(self._plate, self.mesh_resolution, self.plate_n_of_subdivision)
		meshs = mesher.make_stiffener_mesh()
		new_mesh = om.TriMesh()
		for mesh in meshs:
			new_mesh = join_tria_meshes(new_mesh,mesh)
		self.mesh = new_mesh


class StiffenedPlateGeometry(GrillageBaseGeometry):
	def __init__(self, plate: Plate, name = '', mesh_resolution = 0.1, plate_n_of_subdivision = 0):
		self._plate = plate
		super().__init__(name,mesh_resolution,plate_n_of_subdivision)

	def _gen_mesh_or_subgeometry(self):
		self.mesh = om.TriMesh()
		self.sub_geometry.clear()
		self.sub_geometry.append(PlateGeometry(self._plate, 'Plate',
														   self.mesh_resolution,self.plate_n_of_subdivision))
		self.sub_geometry.append(StiffenersGeometry(self._plate, 'Stiffeners',
														   self.mesh_resolution,self.plate_n_of_subdivision))


class PlatingGeometry(GrillageBaseGeometry):
	def __init__(self, plates: Plate, name = '', mesh_resolution = 0.1, plate_n_of_subdivision = 0):
		self._plates = plates
		super().__init__(name,mesh_resolution,plate_n_of_subdivision)

	def _gen_mesh_or_subgeometry(self):
		self.mesh = om.TriMesh()
		self.sub_geometry.clear()
		for id,plate in self._plates.items():
			self.sub_geometry.append(StiffenedPlateGeometry(plate, 'Stiff_plate_'+str(id),
														   self.mesh_resolution,self.plate_n_of_subdivision))

class BeamSegmentGeometry(GrillageBaseGeometry):
	def __init__(self, segment: Segment, name = '', mesh_resolution = 0.1, plate_n_of_subdivision = 0):
		self._segment = segment
		super().__init__(name,mesh_resolution,plate_n_of_subdivision)

	def _gen_mesh_or_subgeometry(self):
		self.sub_geometry.clear()
		mesher = BeamSegmentMesher(self._segment, self.mesh_resolution, self.plate_n_of_subdivision)
		self.mesh = mesher.make_beam_segment_mesh()

class PrimarySuppMemGeometry(GrillageBaseGeometry):
	def __init__(self, psm: PrimarySuppMem, name = '', mesh_resolution = 0.1, plate_n_of_subdivision = 0):
		self._psm:PrimarySuppMem = psm
		super().__init__(name,mesh_resolution,plate_n_of_subdivision)

	def _gen_mesh_or_subgeometry(self):
		self.mesh = om.TriMesh()
		self.sub_geometry.clear()
		for beam_segment in self._psm._segments:
			self.sub_geometry.append(BeamSegmentGeometry(beam_segment, 'Segment_'+str(beam_segment.id),
														   self.mesh_resolution,self.plate_n_of_subdivision))


class Longs_or_Trans_Geometry(GrillageBaseGeometry):
	def __init__(self, items: Dict[int,PrimarySuppMem], name = '', mesh_resolution = 0.1, plate_n_of_subdivision = 0):
		self._items: Dict[int, PrimarySuppMem] = items
		super().__init__(name,mesh_resolution,plate_n_of_subdivision)

	def _gen_mesh_or_subgeometry(self):
		self.mesh = om.TriMesh()
		self.sub_geometry.clear()
		for id, psm in self._items.items():
			self.sub_geometry.append(PrimarySuppMemGeometry(psm,'Beam_'+str(id),
														   self.mesh_resolution,self.plate_n_of_subdivision))

class LongitudinalsGeometry(Longs_or_Trans_Geometry):

	def __init__(self, items: Dict[int,PrimarySuppMem], name = '', mesh_resolution = 0.1, plate_n_of_subdivision = 0):
		super().__init__(items,name,mesh_resolution,plate_n_of_subdivision)

class TransversalsGeometry(Longs_or_Trans_Geometry):

	def __init__(self, items: Dict[int,PrimarySuppMem], name = '', mesh_resolution = 0.1, plate_n_of_subdivision = 0):
		super().__init__(items,name,mesh_resolution,plate_n_of_subdivision)


class GrillageGeometry(GrillageBaseGeometry):

	def __init__(self, grillage: Grillage, name = '', mesh_resolution = 20, plate_n_of_subdivision = 0):
		self._grill = grillage
		super().__init__(name,mesh_resolution,plate_n_of_subdivision)

	def _gen_mesh_or_subgeometry(self):
		try:
			self.mesh=om.TriMesh()
			self.sub_geometry.clear()
			self.sub_geometry.append(LongitudinalsGeometry(self._grill.longitudinal_members(),'Longitudinal Beams',
														   self.mesh_resolution,self.plate_n_of_subdivision))
			self.sub_geometry.append(TransversalsGeometry(self._grill.transverse_members(),'Transverse Beams',
														   self.mesh_resolution,self.plate_n_of_subdivision))
			self.sub_geometry.append(PlatingGeometry(self._grill.plating(), 'Plating',
														   self.mesh_resolution,self.plate_n_of_subdivision))
			self.mesh = self._appendmesh(self.mesh)
		except BaseException as error:
			print('An exception occurred during visualization mesh generation: {}'.format(error))
		except:
			print('Unknown exception occurred during visualization mesh generation')




















































