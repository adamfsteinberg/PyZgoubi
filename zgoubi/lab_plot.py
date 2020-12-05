"""Module for drawing plots of a lattice in lab coordinates

"""

from __future__ import division, print_function
from math import *
import numpy as np
from zgoubi.core import zlog
import scipy.interpolate
import scipy.spatial
import os


null_elements = "END FAISCEAU FAISCNL FAISTORE MCOBJET OBJET PARTICUL REBELOTE MARKER".split()
rect_elements = "DRIFT MULTIPOL QUADRUPO BEND TOSCA CAVITE ELMULT".split()

def get_param(element, name, fallback=None):
	if hasattr(element, "plot_hints") and name in element.plot_hints:
		return element.plot_hints[name]
	elif hasattr(element, name):
		return element.get(name)
	elif fallback is not None:
		return fallback
	raise ValueError("Element: %s Type: %s has no atrribute or plot hint for %s"%(element.label1, element._zgoubi_name,name))



class LabPlotElement(object):
	def __init__(self, z_element, prev_coord, prev_angle, boro, sector_width, line):
		self.z_element = z_element
		self.prev_coord = list(prev_coord)
		self.prev_angle = prev_angle

		self.entry_coord = list(prev_coord) # coord of the entry of the ref line (exit of previous element)
		self.entry_angle = prev_angle # angle of the entry of the ref line
		self.element_type =  z_element._zgoubi_name

		zlog.debug("LabPlotElement(%s %s %s)"%(z_element._zgoubi_name, prev_coord, prev_angle))

		self.exit_coord = list(self.entry_coord) # coord of the exit of the ref line
		self.exit_angle = self.entry_angle #  angle of the exit of the ref line

		if self.element_type == "TOSCA":
			if hasattr(self.z_element, "plot_hints") and "AT" in self.z_element.plot_hints:
				print("polar")
				# Then is is a polar TOSCA, so change type
				self.element_type = "TOSCAp"
			elif hasattr(self.z_element, "plot_hints") and "XL" in self.z_element.plot_hints:
				# Then is rectagular TOSCA
				pass
			else:
				raise ValueError("TOSCA does not supply enough information for plotting. Use set_plot_hint() to set AT or XL")

		if self.element_type in null_elements:
			# does not effect locations
			pass

		elif self.element_type in rect_elements:
			# step by length, no angle change
			if self.element_type == "CAVITE":
				self.length = 0
				if self.z_element.IOPT == 10: # IOPT=10 has length in m
					self.length = get_param(self.z_element,"L") *100
			else:
				self.length = get_param(self.z_element,"XL")

			self.exit_coord = self.transform(self.length,0)
			self.entrance_wedge_angle = 0
			self.exit_wedge_angle = 0

			if self.element_type in "MULTIPOL QUADRUPO".split():
				self.width_p =  self.z_element.R_0
				self.width_m =  self.z_element.R_0
			else:
				width = get_param(self.z_element,"width", 20)
				try:
					self.width_p = float(width)/2
					self.width_m = float(width)/2
				except TypeError:
					self.width_p = float(width[0])
					self.width_m = float(width[1])

			if self.element_type in "MULTIPOL QUADRUPO".split():
				if self.z_element.KPOS == 1:
					pass
				elif self.z_element.KPOS == 2:
					if self.z_element.XCE != 0 or self.z_element.ALE != 0:
						raise ValueError("Only %s with KPOS=2 with YCE is currently implemented"%self.element_type)
					self.entry_coord[0] += self.z_element.YCE * -sin(self.entry_angle)
					self.entry_coord[1] += self.z_element.YCE * cos(self.entry_angle)
				elif self.z_element.KPOS == 3:
					if self.z_element.ALE == 0 or self.z_element.XCE != 0 or self.z_element.YCE != 0:
						raise ValueError("Only %s with KPOS=3 with ALE is currently implemented"%self.element_type)
					self.entry_angle += self.z_element.ALE
					self.exit_angle += self.z_element.ALE
					self.exit_coord = self.transform(self.length,0)
				else:
					raise ValueError("Only %s with KPOS=1 or 2 is currently implemented"%self.element_type)
			elif self.element_type in "BEND":
				if boro is None : raise ValueError("Must set boro for %s"%self.element_type)
				angle = asin(self.z_element.B1 * self.z_element.XL / 2 / boro)
				self.entrance_wedge_angle = self.z_element.W_E - angle
				self.exit_wedge_angle = self.z_element.W_S - angle 
				if self.z_element.KPOS  == 1:
					pass
				elif self.z_element.KPOS  == 3:
					self.entry_angle -= angle
					self.exit_angle -= 2*angle
					self.exit_coord = self.transform(self.length,0)
				else:
					raise ValueError("Only %s with KPOS=1 or 3 is currently implemented"%self.element_type)

		elif self.element_type == "CHANGREF":
			# change reference element
			# angle
			self.exit_angle += radians(self.z_element.ALE)
			# X
			self.exit_coord[0] += self.z_element.XCE * cos(self.exit_angle)
			self.exit_coord[1] += self.z_element.XCE * sin(self.exit_angle)
			# Y
			self.exit_coord[0] += self.z_element.YCE * -sin(self.exit_angle)
			self.exit_coord[1] += self.z_element.YCE * cos(self.exit_angle)


		elif self.element_type in ["DIPOLE", "DIPOLES", "FFAG", "POLARMES", "TOSCAp"]:
			if self.element_type in ["POLARMES"]:
				# to get opening angle read the list of angles in the field map
				self.fmap_file_path = self.z_element.FNAME
				try:
					fmap_file = open(self.fmap_file_path)
				except IOError:
					for fname in line.input_files:
						if os.path.basename(fname) == self.z_element.FNAME:
							self.fmap_file_path = fname
							fmap_file = open(self.fmap_file_path)
				for n in range(4): dummy = fmap_file.readline()
				phi_line = fmap_file.readline()
				phi_line = phi_line.split()
				self.dip_at = float(phi_line[-1])
			else:
				self.dip_at = radians(get_param(self.z_element,"AT"))

			self.dip_re = self.z_element.RE # radius at entry
			self.dip_rs = self.z_element.RS # radius at exit
			self.dip_te = self.z_element.TE # radius at entry
			self.dip_ts = self.z_element.TS # radius at exit

			width = min(self.dip_re, 250)*2 # default
			if sector_width:
				width = sector_width # override with function argument
			width = get_param(self.z_element,"width",width) # override with plot_hint
			try:
				self.width_p = float(width)/2
				self.width_m = float(width)/2
			except TypeError:
				self.width_p = float(width[0])
				self.width_m = float(width[1])

			if self.dip_te != 0 or self.dip_ts != 0:
				raise ValueError("Non zero TE or TS not implemented for %s"%self.element_type)

			# assume local coords have origin at entry, with Y in middle of magnet pointing and machine center
			self.entry_angle += self.dip_at / 2
			self.exit_angle -= self.dip_at
			self.sector_center_coord = list(self.entry_coord)
			self.sector_center_coord[0] += self.dip_re * sin(self.entry_angle - self.dip_at/2)
			self.sector_center_coord[1] -= self.dip_re * cos(self.entry_angle - self.dip_at/2)

			self.exit_coord[0] = self.sector_center_coord[0] + self.dip_rs * sin(-self.exit_angle)
			self.exit_coord[1] = self.sector_center_coord[1] + self.dip_rs * cos(-self.exit_angle)

			self.sub_boundaries = []
			if self.element_type in ["DIPOLES", "FFAG"]:
				for sub_element in self.z_element.subelements:
					e = radians(sub_element["ACN"] - sub_element["OMEGA_E"])
					s = radians(sub_element["ACN"] - sub_element["OMEGA_S"])
					ea = sub_element["THETA_E"]
					sa = sub_element["THETA_S"]

					self.sub_boundaries.append(dict(e=e, s=s, ea=ea, sa=sa))

		else:
			raise ValueError("Can't handle element "+ self.element_type)
	
	def transform(self, x, y):
		if self.element_type not in ["DIPOLE", "DIPOLES", "POLARMES", "FFAG","TOSCAp"]:
			x0, y0 = self.entry_coord # FIXME how to handle transform in changref
			a0 = self.entry_angle

			x1 = x0 + x * cos(a0) - y * sin(a0)
			y1 = y0 + y * cos(a0) + x * sin(a0)
		else:
			# coords are in polar
			x0, y0 = self.sector_center_coord
			a0 = self.prev_angle
			if self.element_type == "TOSCAp":
				a0-=self.dip_at/2

			x1 = x0 + y * sin(-a0 + x)
			y1 = y0 + y * cos(-a0 + x)
		return [x1, y1]

	def draw_ref_line(self, lpd, style):
		t = self.transform
		if self.element_type in rect_elements:
			points = [t(0,0), t(self.length,0)]
			xs, ys = zip(*points)
			lpd.draw_line(xs, ys, **style["reference"])
		if self.element_type in ["DIPOLE", "DIPOLES", "FFAG", "POLARMES","TOSCAp"]:
			re = self.dip_re
			a = self.dip_at
			ang_shift = 0
			if self.element_type == "TOSCAp":
				ang_shift = -self.dip_at/2
			arcsteps = 20
			points = []
			for a1 in np.linspace(0,a,arcsteps):
				points.append(t(a1+ang_shift, re))
			xs, ys = zip(*points)
			lpd.draw_line(xs, ys, **style["reference"])

	
	def draw_outline(self, lpd, style):
		t = self.transform
		if ((self.element_type in rect_elements) and (self.width_p != 0 and self.width_m != 0 )) :
			if self.entrance_wedge_angle == 0 and self.exit_wedge_angle == 0:
				points = [t(0,self.width_p),
						  t(0,-self.width_m),
						  t(self.length,-self.width_m),
						  t(self.length,self.width_p),
						  t(0,self.width_p)]
			elif abs(self.entrance_wedge_angle) > pi/2:
				raise ValueError("entrance_wedge_angle of %s greater than pi/2"%self.element_type)
			elif abs(self.exit_wedge_angle) > pi/2:
				raise ValueError("exit_wedge_angle of %s greater than pi/2"%self.element_type)
			else:
				entry_offset_p = self.width_p * sin(self.entrance_wedge_angle)
				exit_offset_p = self.width_p * sin(self.exit_wedge_angle)
				entry_offset_m = self.width_m * sin(self.entrance_wedge_angle)
				exit_offset_m = self.width_m * sin(self.exit_wedge_angle)
				points = [t(0+entry_offset_p,self.width_p),
						  t(0-entry_offset_m,-self.width_m),
						  t(self.length+exit_offset_m,-self.width_m),
						  t(self.length-exit_offset_p,self.width_p),
						  t(0+entry_offset_p,self.width_p)]

			xs, ys = zip(*points)
			lpd.draw_line(xs, ys, **style["magnet_outline"])

		if self.element_type in ["DIPOLE", "DIPOLES","POLARMES", "FFAG","TOSCAp"]:
			# in polar
			re = self.dip_re
			wp = self.width_p
			wm = -self.width_m
			a = self.dip_at
			ang_shift = 0
			if self.element_type == "TOSCAp":
				ang_shift = -self.dip_at/2
			arcsteps = 20
			points = [ t(0+ang_shift, re + wm),
			           t(0+ang_shift, re + wp)]
			for a1 in np.linspace(0,a,arcsteps):
				points.append(t(a1+ang_shift, re + wp))
			points.append(t(a+ang_shift, re + wp))
			points.append(t(a+ang_shift, re + wm))
			for a1 in np.linspace(a,0,arcsteps):
				points.append(t(a1+ang_shift, re + wm))

			xs, ys = zip(*points)

			if self.element_type in ["DIPOLES", "FFAG"]:
				lpd.draw_line(xs, ys, **style["element_outline"])
			else:
				lpd.draw_line(xs, ys, **style["magnet_outline"])

			if self.element_type in ["DIPOLES", "FFAG"]:
				#sub element boundaries 
				for seb in self.sub_boundaries:
					if seb['ea'] or seb['sa']:
						zlog.warn("Drawing sub magnet boundaries with THETA != 0 not implemented")
						continue
					points = [ t(seb['e'], re + wm),
							   t(seb['e'], re + wp)]
					for a1 in np.linspace(seb['e'],seb['s'],arcsteps):
						points.append(t(a1, re + wp))
					points.append(t(seb['s'], re + wp))
					points.append(t(seb['s'], re + wm))
					for a1 in np.linspace(seb['s'],seb['e'],arcsteps):
						points.append(t(a1, re + wm))

					xs, ys = zip(*points)
					lpd.draw_line(xs, ys, **style["magnet_outline"])



		
class LabPlotDrawer(object):
	def __init__(self, mode="matplotlib", aspect="equal", plot_extents=None):
		self.mode = mode
		self.aspect = aspect
		self.plot_extents = plot_extents
		
		if self.mode == "matplotlib":
			global matplotlib,plt, Line2D
			import matplotlib
			import matplotlib.pyplot as plt
			from matplotlib.lines import Line2D

			self.fig = plt.figure()
			self.fig.clf()
			self.ax = self.fig.add_subplot(111)
			self.ax.set_aspect(self.aspect, adjustable='datalim')
			if plot_extents is not None:
				self.ax.set_aspect(self.aspect, adjustable='box')
		else:
			raise ValueError("Can't handle mode "+ self.mode)

	
	def draw_line(self, xs, ys, color="k", linestyle="-", linewidth=1):
		xs = np.array(xs)
		ys = np.array(ys)
		#if np.any(xs < -10): raise ValueError
		if self.mode == "matplotlib":
			self.ax.plot(xs, ys, color=color, linestyle=linestyle, linewidth=linewidth)
			
		else: ValueError("Can't handle mode "+ self.mode)
	
	def draw_label(self, x, y, l, marker=""):
		if marker != "":
			self.ax.plot([x],[y], marker)
		self.ax.annotate(str(l), (x,y))
	
	def draw_im(self, im, extent, cm, vmin, vmax, colorbar=True, colorbar_label=""):
		i = self.ax.imshow(im, extent=extent, origin='lower', cmap=plt.get_cmap(cm),vmin=vmin, vmax=vmax, aspect=self.aspect)
		if colorbar:
			cbar = self.fig.colorbar(i)
			if colorbar_label:
				cbar.ax.set_ylabel(colorbar_label)

	def finish(self):
		if self.mode == "matplotlib":
			plt.xlabel("Lab x (cm)")
			plt.ylabel("Lab y (cm)")
			version = matplotlib.__version__.split(".")
			if int(version[0]) >= 1 and int(version[1]) >= 1:
				plt.tight_layout()
			if self.plot_extents:
				self.ax.set_xlim(self.plot_extents[0:2])
				self.ax.set_ylim(self.plot_extents[2:4])

	def show(self):
		if self.mode == "matplotlib":
			plt.show()
		else: ValueError("Can't handle mode "+ self.mode)

	def save(self,fname):
		if self.mode == "matplotlib":
			if self.plot_extents is not None:
				plt.savefig(fname, bbox_inches="tight")
			else:
				plt.savefig(fname)
		else: ValueError("Can't handle mode "+ self.mode)


class LabPlot(object):
	"""A plotter for beam lines and tracks.
	
	"""
	def __init__(self, line, boro=None, sector_width=None, aspect="equal", style=None):
		"""Creates a new plot from the line.
		If using an element that adjusts it shape based on BORO, then it must be passed in.
		If sector_width is a number it used for the width of sector elements, alternatively if sector_width is a list of 2 values they are used as the distance to the boundary in the positive and negative directions

		"""
		
		self.line = line
		self.boro = boro
		self.sector_width = sector_width
		self.aspect = aspect
		self.style = {}
		self.noel_offset = 0
		self.style["track"] = dict(color="r", linestyle="-", linewidth=0.1)
		self.style["magnet_outline"] = dict(color="b", linestyle="-", linewidth=1)
		self.style["element_outline"] = dict(color="b", linestyle=":", linewidth=1)
		self.style["reference"] = dict(color="k", linestyle="-", linewidth=1)

		if style:
			self.set_style(style)
		
		self._scan_line()
		self.lpd = None

	def set_noel_offset(self, offset):
		"If the line passed to LabPlot is missing some initial elements line used for tracking, then set the offset to the number of missing elements so that the element numbers can be synced"
		self.noel_offset = offset

	def _scan_line(self):
		"""Scan through the line, and make a note of where all the elements are.
		
		"""
		self.elements = []
		self.element_label1 = []
		self.tracks = []
		self.mag_tracks = []
		self.field_map_data = []
		self.duped_labels = []
		angle = 0
		position = [0, 0]

		for elem in self.line.elements():
			if hasattr(elem, 'label1'):
				label = elem.label1
				if label in self.element_label1:
					self.duped_labels.append(label)
					#zlog.warn("Repeated label '%s'"%label)
			else:
				label = ""
			lpelem = LabPlotElement(elem, position, angle, boro=self.boro, sector_width=self.sector_width, line=self.line)
			self.elements.append(lpelem)
			self.element_label1.append(label)
			angle = lpelem.exit_angle
			position = lpelem.exit_coord

	def set_style(self, style):
		"""Style is a nested dictionary, that updates the defaults, eg::

		   style = {"track":{"color":"g"}}

		Elements to be styled include "track", "magnet_outline", "element_outline", "reference". Each can take a "color", "linestyle" and "linewidth", using matplotlib notations.

		"""
		for k,v in style.items():
			self.style[k].update(v)

	def draw(self, draw_tracks=True, draw_field_points=False, draw_field_midplane=False, field_component='z', field_steps=100, field_int_mode="kd", plot_extents=None):
		"""Draw the plot in memory, then use :py:meth:`show()` to display to screen or :py:meth:`save()` to save to file.

		plot_extents: list of extents [left, right, bottom, top]
		
		"""
		if self.lpd is None:
			self.lpd = LabPlotDrawer(aspect=self.aspect, plot_extents=plot_extents)

		if field_component not in ['x','y','z']:
			raise ValueError("field_component should be 'y', 'z' or 'x'")

		if draw_field_midplane:
			if field_int_mode=="griddata":
				if not hasattr(scipy.interpolate, "griddata"):
					print("LabPlot.draw() with field_int_mode='griddata' requires scipy > 0.9, try field_int_mode='kd'")
					raise
				label = r"$B_%s$ (kG)"%field_component
				field_map_data = np.array(self.field_map_data)
				points = field_map_data[:,0:3:2]
				if field_component == 'y':
					values = field_map_data[:,3].reshape([-1])
				elif field_component == 'z':
					values = field_map_data[:,4].reshape([-1])
				elif field_component == 'x':
					values = field_map_data[:,5].reshape([-1])

				# http://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.griddata.html
				# interested in y-x plane
				ymin,ymax,xmin,xmax = points[:,0].min(), points[:,0].max(), points[:,1].min(), points[:,1].max()
				nsteps = field_steps*1j # j because of how mgrid works
				grid_y, grid_x = np.mgrid[ymin:ymax:nsteps, xmin:xmax:nsteps]
				
				# check that paths don't deviate from z=0 plane
				if np.abs(field_map_data[:,1]).max() > 1e-6:
					zlog.warn("Some field points are not at z=0. Plot will be of projection onto z=0")

				int_field = scipy.interpolate.griddata(points, values, (grid_y, grid_x))
				#print int_field.nanmax(), int_field.nanmin(),  np.abs(int_field).nanmax()
				vmax = np.nanmax(np.abs(int_field))
				self.lpd.draw_im(int_field.T, (ymin,ymax,xmin,xmax), cm='bwr', vmin=-vmax, vmax=vmax, colorbar_label=label)
			elif field_int_mode=="kd":			
				label = r"$B_%s$ (kG)"%field_component
				field_map_data = np.array(self.field_map_data)
				points = field_map_data[:,0:3:2]
				if field_component == 'y':
					values = field_map_data[:,3].reshape([-1])
				elif field_component == 'z':
					values = field_map_data[:,4].reshape([-1])
				elif field_component == 'x':
					values = field_map_data[:,5].reshape([-1])

				xmin,xmax,ymin,ymax = points[:,0].min(), points[:,0].max(), points[:,1].min(), points[:,1].max()

				#print points.shape
				#print points
				#print values.shape
				#print xmin,xmax,ymin,ymax

				nxsteps = int(field_steps * 2)
				xstep_size = (xmax-xmin) / nxsteps
				nysteps = int(max(1, (ymax-ymin) / xstep_size))

				kd = scipy.spatial.cKDTree(points)

				field_map = np.zeros([nxsteps, nysteps])
				for nx,x in enumerate(np.linspace(xmin,xmax,nxsteps)):
					for ny,y in enumerate(np.linspace(ymin,ymax,nysteps)):
						# get nearby points
						kd_f, kd_i = kd.query([x,y],k=5, distance_upper_bound=xstep_size*1.2)
						# remove index outside range, then mean not points found in distance_upper_bound
						kd_i = kd_i[kd_i < values.size]
						# remove zero fields
						near_values = values[kd_i][ values[kd_i] != 0 ]
						if near_values.size:
							field_map[nx,ny] = near_values.mean()
				vmax = np.nanmax(np.abs(field_map))
				self.lpd.draw_im(field_map.T, (xmin,xmax,ymin,ymax), cm='bwr', vmin=-vmax, vmax=vmax, colorbar_label=label)
			else:
				raise ValueError("field_int_mode must be griddata or kd")



		for elem in self.elements:
			elem.draw_ref_line(self.lpd, self.style)
			elem.draw_outline(self.lpd, self.style)

		if draw_tracks:
			for track in self.tracks:
				xs, ys, dummy, dummy, dummy = zip(*track)
				self.lpd.draw_line(xs, ys, **self.style["track"])

		if draw_field_points:
			for track in self.mag_tracks:
				for p in track:
					xs, ys, by, bz, bx = p
					if by is None: continue
					if field_component == 'y':
						self.lpd.draw_label(xs, ys, "%.3g"%by, 'rx')
					elif field_component == 'z':
						self.lpd.draw_label(xs, ys, "%.3g"%bz, 'rx')
					elif field_component == 'x':
						self.lpd.draw_label(xs, ys, "%.3g"%bx, 'rx')

		self.lpd.finish()
	
	def show(self):
		"Display plot to screen, call :py:meth:`draw()` first."
		self.lpd.show()

	def save(self,fname):
		"Save plot to file, call :py:meth:`draw()` first."
		self.lpd.save(fname)


	def add_tracks(self, ftrack=None, ptrack=None, draw=1, field=1):
		"Add tracks from plt or fai files."
		#tracks = []
 		# find the list of particles and laps/passes
		pids = set()
		passes = set()
		noels = set()
		if ftrack is not None:
			pids  |= set(np.unique(ftrack['ID']))
			passes  |= set(np.unique(ftrack['PASS']))
			noels  |= set(np.unique(ftrack['NOEL']))
		if ptrack is not None:
			pids |= set(np.unique(ptrack['ID']))
			passes  |= set(np.unique(ptrack['PASS']))
			noels  |= set(np.unique(ptrack['NOEL']))
		if ftrack is None and ptrack is None:
			raise ValueError("Must pass a fai track or plt track (or both)")

		pids = sorted(list(pids))
		passes = sorted(list(passes))
		noels = sorted(list(noels))

		# dummy tracks if not passed
		dummy = np.zeros([0], dtype=[('ID',int),('PASS',int),('NOEL',int),('element_label1','S10')])
		if ftrack is None: ftrack = dummy
		if ptrack is None: ptrack = dummy


		# per particle, per lap, per element
		for pid in pids:
			ftrack_p = ftrack[ftrack['ID'] == pid]
			ptrack_p = ptrack[ptrack['ID'] == pid]
			for ppass in passes:
				ftrack_pp = ftrack_p[ftrack_p['PASS'] == ppass]
				ptrack_pp = ptrack_p[ptrack_p['PASS'] == ppass]
				this_track = []
				for noel in noels:
					ftrack_ppn = ftrack_pp[ftrack_pp['NOEL'] == noel]
					ptrack_ppn = ptrack_pp[ptrack_pp['NOEL'] == noel]
					if len(ftrack_ppn) == 0 and len(ptrack_ppn) == 0:
						continue

					el_ind = noel - self.noel_offset - 1

					for t in ptrack_ppn:
						if t['IEX'] != 1: break
						y = t['Y']
						x = t['X']
						z = t['Z']
						if isnan(x) or isnan(y) or isnan(z):
							continue
						by, bz, bx = t['BY'], t['BZ'], t['BX']
						if abs(by) < 1e-50 and abs(bz) < 1e-50 and abs(bx) < 1e-50:
							by, bz, bx = 0, 0, 0 # ignore tiny fields
						xt, yt = self.elements[el_ind].transform(x,y)
						this_track.append([xt,yt, by, bz, bx])
						if field:
							self.field_map_data.append([xt,z,yt,by, bz, bx])
					for t in ftrack_ppn:
						if t['IEX'] != 1: break
						# fai has no x coord, and takes label from element before it
						y = t['Y']
						x = 0

						xt, yt = self.elements[el_ind+1].transform(x,y)
						this_track.append([xt,yt,None,None,None])

				#print this_track
				if len(this_track) > 0 :
					if draw:
						self.tracks.append(this_track)
					if field:
						self.mag_tracks.append(this_track)


					



