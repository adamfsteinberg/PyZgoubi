#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""A bunch object to hold the coordinates for many particles

Note that all values are in SI units, m, rad, eV, s
"""

from __future__ import division, print_function
from math import *
import numpy
from zgoubi import rel_conv
from zgoubi import io
from zgoubi.core import zlog, dep_warn
import struct
import inspect
import itertools
import warnings

#from zgoubi.utils import *


class Bunch(object):
	"""Object to store a bunch of particles efficiently using numpy.
	All values are in SI units, m, rad, eV, s
	Can be use to create a bunch of particles with all there coordinates set to zero (appart from the D coordinate which is set to 1)::
	
		my_bunch = Bunch(nparticles=100, ke=1e6, mass=PROTON_MASS, charge=1)
	
	It can also be used to create a bunch from some existing coordinates stored in a numpy array::
	
		 my_bunch = Bunch(ke=1e6, mass=PROTON_MASS, charge=1, particles=existing_coords)
	
	"""
	min_data_def = [
	('D', numpy.float64), # these coorspond to zgoubi D,Y,T,Z,P,S, but in SI units
	('Y', numpy.float64),
	('T', numpy.float64),
	('Z', numpy.float64),
	('P', numpy.float64),
	('S', numpy.float64),
	('tof', numpy.float64), # these are for accumulating
	('X', numpy.float64),
	]

	def __init__(self, nparticles=0, ke=None, rigidity=0, mass=0, charge=1, particles=None):
		"""Bunch constructor.

		"""
		if particles is not None:
			self.coords = particles
		else:
			nparticles = int(nparticles)
			self.coords = numpy.zeros(nparticles, self.min_data_def)
			self.coords['D'] = 1
		self.mass = mass
		self.charge = charge
		self.rigidity = rigidity
		if ke is not None:
			self.set_bunch_ke(ke)

	def split_bunch(self, max_particles, n_slices):
		"Split a bunch into n_slices smaller bunches, or more if they would have too many particles in."
		if ceil(len(self.coords) / n_slices) > max_particles:
			n_slices = ceil(len(self.coords)/max_particles)

		rigidity = self.get_bunch_rigidity()
		for pslice in numpy.array_split(self.coords, n_slices):
			if pslice.size != 0:
				yield Bunch(rigidity=rigidity, mass=self.mass, charge=self.charge,
			            particles=pslice)

	def __str__(self):
		out = "Bunch:\n"
		out += "\t"+str(len(self)) + " paricles\n"
		out += "\t"+str(self.get_bunch_ke()) + " eV\n"

		return out

	def set_bunch_ke(self, ke):
		"Set bunch kinetic energy"
		if self.mass == 0:
			raise ValueError("Particle mass can't be Zero")
		if self.charge == 0:
			raise ValueError("Particle charge can't be Zero")
		self.rigidity = rel_conv.ke_to_rigidity(mass=self.mass, ke=ke, charge=self.charge)

	def get_bunch_ke(self):
		"Get bunch kinetic energy"
		if self.mass == 0:
			zlog.warn("Particle mass is Zero. Set the mass before getting the KE")
			return 0
		if self.charge == 0:
			zlog.warn("Particle charge is Zero. Set the charge before getting the KE")
			return 0
		return rel_conv.rigidity_to_ke(mass=self.mass, rigidity=self.rigidity, charge=self.charge)

	def set_bunch_rigidity(self, rigidity):
		"Set bunch rigidity"
		self.rigidity = rigidity

	def get_bunch_rigidity(self):
		"Get bunch rigidity"
		return self.rigidity

	def particles(self):
		"Returns the numpy array that holds the coordinates, as a structured numpy array"
		return self.coords

	def raw_particles(self):
		"""Returns the numpy array that holds the coordinates, as a 2d numpy array and a list of comumn names. It is possible that the column names may change or reorder, so if you use this you may want to protect against it with something like
		::
		
			assert(pbunch.raw_particles()[1] == ["D","Y","T","Z","P","S","TOF","X"])

		"""
		return self.coords.view((numpy.float64, len(self.coords.dtype.names))), self.coords.dtype.names

	def get_min_BORO(self):
		"Returns the minimum rigidity of the bunch"
		#min_BORO = ke_to_rigidity(min(self.coords['KE']), self.mass)
		min_BORO = self.rigidity * self.coords['D'].min()
		return min_BORO

	@staticmethod
	def gen_halo_x_xp_y_yp(npart, emit_y, emit_z, beta_y, beta_z, alpha_y, alpha_z, seed=None,
			               ke=None, rigidity=0, mass=0, charge=1):
		"""Generate a halo bunch, i.e. an elipse outline in x-xp (Y-T) and in y-yp (Z-P). S and D are set to 0 and 1 respectively.
		example::
		
			my_bunch = Bunch.gen_halo_x_xp_y_yp(1000, 1e-3, 1e-3, 4, 5, 1e-3, 2e-2, ke=50e6, mass=PROTON_MASS, charge=1)
		
		creates a halo bunch called my_bunch with 1000 particles of the given parameters.

		"""

		if emit_y < 0 or emit_z < 0 or beta_y < 0 or beta_z < 0:
			print("Emittance or beta can't be negative")
			print("emit_y, emit_z, beta_y, beta_z, alpha_y, alpha_z")
			print(emit_y, emit_z, beta_y, beta_z, alpha_y, alpha_z)
			raise ValueError

		npart = int(npart)

		ry = sqrt(emit_y) 
		rz = sqrt(emit_z) 

		if seed is not None:
			numpy.random.seed(seed)

		u1 = numpy.random.random_sample([npart]) * pi * 2
		u2 = numpy.random.random_sample([npart]) * pi * 2

		coords = numpy.zeros([npart, 6], numpy.float64)

		coords[:, 0] = ry * numpy.cos(u1)
		coords[:, 1] = ry * numpy.sin(u1)
		coords[:, 2] = rz * numpy.cos(u2)
		coords[:, 3] = rz * numpy.sin(u2)

		matrix = Bunch._twiss_matrix(beta_y, beta_z, alpha_y, alpha_z)
		
		for n, coord in enumerate(coords):
			#	new_coord = numpy.dot(coord, matrix)
			coords[n] = numpy.dot(matrix, coord)
		
		bunch = numpy.zeros([npart], Bunch.min_data_def)
		
		bunch['Y'] = coords[:, 0]
		bunch['T'] = coords[:, 1]
		bunch['Z'] = coords[:, 2]
		bunch['P'] = coords[:, 3]
		bunch['D'] = 1
		
		return Bunch(ke=ke, rigidity=rigidity, mass=mass, charge=charge, particles=bunch)

	@staticmethod
	def gen_kv_x_xp_y_yp(npart, emit_y, emit_z, beta_y, beta_z, alpha_y, alpha_z, seed=None,
			               ke=None, rigidity=0, mass=0, charge=1):
		"""Generate a uniform (KV) bunch, i.e. a surface of a 4D hypershere in x-xp-y-yp (Y-T-Z-P). S and D are set to 0 and 1 respectively.
		example::
		
			my_bunch = Bunch.gen_kv_x_xp_y_yp(1000, 1e-3, 1e-3, 4, 5, 1e-3, 2e-2, ke=50e6, mass=PROTON_MASS, charge=1)
		
		creates a KV bunch called my_bunch with 1000 particles of the given parameters.

		"""

		if emit_y < 0 or emit_z < 0 or beta_y < 0 or beta_z < 0:
			print("Emittance or beta can't be negative")
			print("emit_y, emit_z, beta_y, beta_z, alpha_y, alpha_z")
			print(emit_y, emit_z, beta_y, beta_z, alpha_y, alpha_z)
			raise ValueError

		if seed is not None:
			numpy.random.seed(seed)

		npart = int(npart)

		ry = sqrt(emit_y)
		rz = sqrt(emit_z)

		#u1 = numpy.random.random_sample([npart]) * pi * 2
		#u2 = numpy.random.random_sample([npart]) * pi * 2

		coords = numpy.zeros([npart, 6], numpy.float64)
	
		# From Chris Prior.
		y = numpy.random.uniform(-1, 1, [npart])
		y2 = numpy.arccos(y) / 2
		phi = numpy.random.uniform(-pi, pi, [npart])
		the = numpy.random.uniform(-pi, pi, [npart])

		coords[:, 0] = ry * numpy.cos(y2) * numpy.cos(phi)
		coords[:, 1] = ry * numpy.cos(y2) * numpy.sin(phi)
		coords[:, 2] = rz * numpy.sin(y2) * numpy.cos(the)
		coords[:, 3] = rz * numpy.sin(y2) * numpy.sin(the)

		matrix = Bunch._twiss_matrix(beta_y, beta_z, alpha_y, alpha_z)
		
		for n, coord in enumerate(coords):
			#	new_coord = numpy.dot(coord, matrix)
			coords[n] = numpy.dot(matrix, coord)
		
		bunch = numpy.zeros([npart], Bunch.min_data_def)
		
		bunch['Y'] = coords[:, 0]
		bunch['T'] = coords[:, 1]
		bunch['Z'] = coords[:, 2]
		bunch['P'] = coords[:, 3]
		bunch['D'] = 1
		
		return Bunch(ke=ke, rigidity=rigidity, mass=mass, charge=charge, particles=bunch)
	@staticmethod
	def gen_waterbag_x_xp_y_yp(npart, emit_y, emit_z, beta_y, beta_z, alpha_y, alpha_z, seed=None,
			               ke=None, rigidity=0, mass=0, charge=1):
		"""Generate a waterbag bunch, i.e. a filled hypersphere in x-xp-y-yp (Y-T-Z-P). S and D are set to 0 and 1 respectively.
		example::
		
			my_bunch = Bunch.gen_waterbag_x_xp_y_yp(1000, 1e-3, 1e-3, 4, 5, 1e-3, 2e-2, ke=50e6, mass=PROTON_MASS, charge=1)
		
		creates a waterbag bunch called my_bunch with 1000 particles of the given parameters.

		"""

		if emit_y < 0 or emit_z < 0 or beta_y < 0 or beta_z < 0:
			print("Emittance or beta can't be negative")
			print("emit_y, emit_z, beta_y, beta_z, alpha_y, alpha_z")
			print(emit_y, emit_z, beta_y, beta_z, alpha_y, alpha_z)
			raise ValueError

		if seed is not None:
			numpy.random.seed(seed)

		npart = int(npart)

		ry = sqrt(emit_y)
		rz = sqrt(emit_z)

		#u1 = numpy.random.random_sample([npart]) * pi * 2
		#u2 = numpy.random.random_sample([npart]) * pi * 2

		coords = numpy.zeros([npart, 6], numpy.float64)
		
		while True:
			# fill a hypercube with 3.5 times as many particles as needed
			a1 = numpy.random.uniform(-1, 1, [int(npart*3.5)])
			a2 = numpy.random.uniform(-1, 1, [int(npart*3.5)])
			a3 = numpy.random.uniform(-1, 1, [int(npart*3.5)])
			a4 = numpy.random.uniform(-1, 1, [int(npart*3.5)])
			
			# discard ones that fall out of hypesphere
			r = a1**2+a2**2+a3**2+a4**2
			keep = r <= 1
			
			a1 = a1[keep]
			a2 = a2[keep]
			a3 = a3[keep]
			a4 = a4[keep]

			#print len(a1)
			
			# it is possible that to many particles will be rejected
			if len(a1) >= npart:
				break

		

		coords[:, 0] = ry * a1[0:npart]
		coords[:, 1] = ry * a2[0:npart]
		coords[:, 2] = rz * a3[0:npart]
		coords[:, 3] = rz * a4[0:npart]

		matrix = Bunch._twiss_matrix(beta_y, beta_z, alpha_y, alpha_z)
		
		for n, coord in enumerate(coords):
			#	new_coord = numpy.dot(coord, matrix)
			coords[n] = numpy.dot(matrix, coord)
		
		bunch = numpy.zeros([npart], Bunch.min_data_def)
		
		bunch['Y'] = coords[:, 0]
		bunch['T'] = coords[:, 1]
		bunch['Z'] = coords[:, 2]
		bunch['P'] = coords[:, 3]
		bunch['D'] = 1
		
		return Bunch(ke=ke, rigidity=rigidity, mass=mass, charge=charge, particles=bunch)

	@staticmethod
	def gen_gauss_x_xp_y_yp(npart, emit_y, emit_z, beta_y, beta_z, alpha_y, alpha_z, seed=None,
			               ke=None, rigidity=0, mass=0, charge=1):
		"""Generate a Gaussian bunch in x-xp (Y-T) and in y-yp (Z-P). S and D are set to 0 and 1 respectively.
		example::
		
			my_bunch = Bunch.gen_kv_x_xp_y_yp(1000, 1e-3, 1e-3, 4, 5, 1e-3, 2e-2, ke=50e6, mass=PROTON_MASS, charge=1)
		
		creates a Gaussian bunch called my_bunch with 1000 particles of the given parameters.

		"""

		if emit_y < 0 or emit_z < 0 or beta_y < 0 or beta_z < 0:
			print("Emittance or beta can't be negative")
			print("emit_y, emit_z, beta_y, beta_z, alpha_y, alpha_z")
			print(emit_y, emit_z, beta_y, beta_z, alpha_y, alpha_z)
			raise ValueError

		if seed is not None:
			numpy.random.seed(seed)

		npart = int(npart)

		coords = numpy.zeros([npart, 6], numpy.float64)
		coords[:, 0:4] = numpy.random.normal(0, 0.5, [npart, 4])
		coords[:, 0] *= sqrt(emit_y)
		coords[:, 1] *= sqrt(emit_y)
		coords[:, 2] *= sqrt(emit_z)
		coords[:, 3] *= sqrt(emit_z)


		matrix = Bunch._twiss_matrix(beta_y, beta_z, alpha_y, alpha_z)
		
		for n, coord in enumerate(coords):
			coords[n] = numpy.dot(matrix, coord)
		
		bunch = numpy.zeros([npart], Bunch.min_data_def)

		bunch['Y'] = coords[:, 0]
		bunch['T'] = coords[:, 1]
		bunch['Z'] = coords[:, 2]
		bunch['P'] = coords[:, 3]
		bunch['D'] = 1
	
		return Bunch(ke=ke, rigidity=rigidity, mass=mass, charge=charge, particles=bunch)

	@staticmethod
	def gen_gauss_x_xp_y_yp_s_dp(npart, emit_y, emit_z, beta_y, beta_z, alpha_y, alpha_z, mom_spread=0, bunch_length=0, disp=0, disp_prime=0, seed=None, ke=None, rigidity=0, mass=0, charge=1):
		"""Generate a Gaussian bunch in transverse and longitudinal phase space
		emit_y, emit_z : horizontal and vertical plane geometric emittance (1 sigma)
		beta_y, beta_z : horizontal and vertical betatron function
		alpha_y, alpha_z : horizontal and vertical alpha
		mom_spead : sigma of momentum spread (percentage) 
		bunch_length : sigma of bunch length (metres)
		disp, disp_prime : dispersion and dispersion prime (D') in horizontal plane """


		if emit_y < 0 or emit_z < 0 or beta_y < 0 or beta_z < 0:
			print("Emittance or beta can't be negative")
			print("emit_y, emit_z, beta_y, beta_z, alpha_y, alpha_z")
			print(emit_y, emit_z, beta_y, beta_z, alpha_y, alpha_z)
			raise ValueError


		if seed is not None:
			numpy.random.seed(seed)

		npart=int(npart)

		#generate momentum distribution
		if mom_spread > 0.0:
			mom_dist = numpy.random.normal(1.0, mom_spread/100, npart)
		#generate longitudinal coordinate distribution (bunch length)
		if bunch_length > 0.0:
			s_dist = numpy.random.normal(0.0, bunch_length, npart)


		#use numpy.random.multivariate_normal. 
		#set off-diagonal terms in covariance matrix to zero so that Y is uncorrelated with T (and Z with P)
		cov_yt = [[emit_y, 0], [0, emit_y]]
		ryrt = numpy.random.multivariate_normal((0, 0), cov_yt, (1, npart))[0]
		cov_zp = [[emit_z, 0], [0, emit_z]]
		rzrp = numpy.random.multivariate_normal((0, 0), cov_zp, (1, npart))[0]

		coords = numpy.zeros([npart, 6], numpy.float64)
		#coords2 = numpy.zeros([npart, 6], numpy.float64)

		coords[:, 0] = ryrt[:, 0]
		coords[:, 1] = ryrt[:, 1]
		coords[:, 2] = rzrp[:, 0]
		coords[:, 3] = rzrp[:, 1]
		if mom_spread > 0.0:
			#Last column of coords is delta_p.
			coords[:, 5] = mom_dist-1

		matrix = Bunch._twiss_matrix(beta_y, beta_z, alpha_y, alpha_z, disp, disp_prime)

		for n, coord in enumerate(coords):
			#	new_coord = numpy.dot(coord, matrix)
			coords[n] = numpy.dot(matrix, coord)
		
		bunch = numpy.zeros([npart], Bunch.min_data_def)
		
		bunch['Y'] = coords[:, 0]
		bunch['T'] = coords[:, 1]
		bunch['Z'] = coords[:, 2]
		bunch['P'] = coords[:, 3]
		if mom_spread > 0.0:
			bunch['D'] = mom_dist
		else:
			bunch['D'] = 1

		if bunch_length > 0.0:
			bunch['S'] = s_dist
		else:
			bunch['S'] = 0.0
		
		return Bunch(ke=ke, rigidity=rigidity, mass=mass, charge=charge, particles=bunch)

	@staticmethod
	def _twiss_matrix(beta_y, beta_z, alpha_y, alpha_z, disp=0, disp_prime=0):
		"Create a matrix that will convert a spherical distribution in to one with the correct twiss values"
		B = numpy.eye(6)
		B[0, 0] = sqrt(beta_y)
		B[1, 1] = 1 / sqrt(beta_y)
		B[2, 2] = sqrt(beta_z)
		B[3, 3] = 1 / sqrt(beta_z)

		A = numpy.eye(6)
		A[1, 0] = -alpha_y
		A[3, 2] = -alpha_z

		#M = numpy.dot(A, B)
		M = numpy.dot(B, A)

		#add dispersion and dispersion_prime terms
		M[0, 5] = disp
		M[1, 5] = disp_prime

		return M


	@staticmethod
	def read_YTZPSD(fname, ke=None, rigidity=0, mass=0, charge=1):
		"""Read in a bunch from a file. assumes columns are Y, T, Z, P, X, D, separated by white space::
			my_bunch = Bunch.read_YTZPSD("mybunch.dat", ke=1e9, mass=ELECTRON_MASS, charge=-1)
		
		"""
		dist = numpy.loadtxt(fname)
		if ((dist.size % 6) != 0):
			raise ValueError("Number of values in %s not a multiple of 6"%fname)
		nparts = dist.size // 6
		dist = dist.reshape(nparts, 6)
		coords = numpy.zeros(nparts, Bunch.min_data_def)
		#self.coords['KE'] = ke
		coords['Y'] = dist[:, 0]
		coords['T'] = dist[:, 1]
		coords['Z'] = dist[:, 2]
		coords['P'] = dist[:, 3]
		coords['X'] = dist[:, 4]
		coords['D'] = dist[:, 5]
		return Bunch(ke=ke, rigidity=rigidity, mass=mass, charge=charge, particles=coords)

	def check_bunch(self):
		"Check that the bunch is not empty, and contains finite values"
		if self.coords.size == 0:
			zlog.error("Empty Bunch. Called by %s()", inspect.stack()[1][3])
			return False
		if not numpy.all(numpy.isfinite(self.raw_particles()[0])):
			zlog.error("Non finite coordinates in bunch. Called by %s()", inspect.stack()[1][3])
			return False
	
		return True


	def write_YTZPSD(self, fname, binary=False):
		"Output a bunch, compatible with read_YTZPSD"
		self.check_bunch()
		
		# it ought to be possible to do this:
		#numpy.savetxt(fh, self.coords[['Y','T','Z','P','S','D']])
		# but the field end up in the wrong order, see http://thread.gmane.org/gmane.comp.python.numeric.general/36933

		if binary:
			fh = open(fname, "wb")
			# this is quite optimised
			# rather than call the general io.write_fortran_record(), use a fast special case
			#header
			header_len = 270 # header length changed from 80  to 270 in Zgoubi svn483
			                 # but old and new version will tolerate the longer header
			io.write_fortran_record(fh, "Binary bunch coordinates from pyzgoubi".ljust(header_len))
			io.write_fortran_record(fh, "Y,T,Z,P,S,D".ljust(header_len))
			io.write_fortran_record(fh, " "*header_len)
			io.write_fortran_record(fh, " "*header_len)
			# record length is always the same
			rec_len_r = struct.pack("i", 6*8)
			rec_len_r2 = rec_len_r + rec_len_r

			fh.write(rec_len_r) # write length once before, twice after each record, and then truncate one at the end
			for p in self.coords.view((numpy.float64, 8)): #  tostring is slightly faster if we ignore the dtypew
			#for p in self.coords:
				ps = p.tostring() # tostring is slow, but quicker than manipulating p and using tofile
				fh.write(ps[8:48] + ps[:8] + rec_len_r2) # writing YTZPSD even though array is DYTZPStofX
			fh.seek(-len(rec_len_r), 1)
			fh.truncate()

		else:
			fh = open(fname, "w")
			nparts = len(self.coords)
			dist = numpy.zeros([nparts, 6])
			dist[:, 0] = self.coords['Y']
			dist[:, 1] = self.coords['T']
			dist[:, 2] = self.coords['Z']
			dist[:, 3] = self.coords['P']
			dist[:, 4] = self.coords['S']
			dist[:, 5] = self.coords['D']
			#dist = dist.reshape(nparts * 2, 3)
			fh.write("# ASCII bunch coordinates from pyzgoubi\n#Y,T,Z,P,S,D\n\n\n")
			numpy.savetxt(fh, dist)
		fh.close()

	
	def get_widths(self):
		"Returns the width of the bunch in each dimension Y,T,Z,P,S,D"
		self.check_bunch()
		y_width = numpy.max(self.coords['Y']) - numpy.min(self.coords['Y'])
		t_width = numpy.max(self.coords['T']) - numpy.min(self.coords['T'])
		z_width = numpy.max(self.coords['Z']) - numpy.min(self.coords['Z'])
		p_width = numpy.max(self.coords['P']) - numpy.min(self.coords['P'])
		s_width = numpy.max(self.coords['S']) - numpy.min(self.coords['S'])
		d_width = numpy.max(self.coords['D']) - numpy.min(self.coords['D'])
		return (y_width, t_width, z_width, p_width, s_width, d_width)

	def get_widths_rms(self):
		"Returns the rms width of the bunch in each dimension Y,T,Z,P,S,D"
		self.check_bunch()
		y_width = self.coords['Y'].std()
		t_width = self.coords['T'].std()
		z_width = self.coords['Z'].std()
		p_width = self.coords['P'].std()
		s_width = self.coords['S'].std()
		d_width = self.coords['D'].std()
		return (y_width, t_width, z_width, p_width, s_width, d_width)

	def get_centers(self):
		"Returns the center of the bunch in each dimension Y,T,Z,P,S,D"
		self.check_bunch()
		y_mean = numpy.mean(self.coords['Y'])
		t_mean = numpy.mean(self.coords['T'])
		z_mean = numpy.mean(self.coords['Z'])
		p_mean = numpy.mean(self.coords['P'])
		s_mean = numpy.mean(self.coords['S'])
		d_mean = numpy.mean(self.coords['D'])
		return (y_mean, t_mean, z_mean, p_mean, s_mean, d_mean)

	def __len__(self):
		"Returns length of bunch. Use len(my_bunch)"
		return len(self.coords)

	def plot(self, fname=None, lims=None, add_bunch=None, fmt=None, longitudinal=True):
		"""Plot a bunch, if no file name give plot is displayed on screen. lims can be used to force axis limits eg [lY,lT,lZ,lP,lX,lD] would plot limit plot from -lY to +lY in Y, etc. Additional bunches can be passed, as add_bunch, to overlay onto the same plot.
		fmt can be a list of formats in matplotlib style, eg ['rx', 'bo']
		"""
		import pylab
		if fmt is None:
			fmt = [',']

		self.check_bunch()
		bunches = [self]
		if add_bunch is not None:
			try:
				# if add_bunch is iterable, append all is memebers
				for a_bunch in add_bunch:
					bunches.append(a_bunch)
			except TypeError:
				# otherwise just append it
				bunches.append(add_bunch)

		coordsz = ['Y', 'T', 'Z', 'P', 'X', 'D']
		coords = ["x", "x'", "y", "y'", "s", "p"]
		plot_specs = [None,
				(0, 2, "x-y (Y-Z)"),
				(2, 3, "y-y' (Z-P)"),
				(0, 1, "x-x' (Y-T)"),
				(4, 5, "s-p (X-D)"),
				]

		for n in range(1, 5):
			if n == 4 and longitudinal == False: continue
			x, y, title = plot_specs[n]
			
			pylab.subplot(2, 2, n)
			pylab.grid()
			for abunch, f in zip(bunches, itertools.cycle(fmt)):
				pylab.plot(abunch.coords[coordsz[x]], abunch.coords[coordsz[y]], f)
				pylab.xlabel(coords[x])
				pylab.ylabel(coords[y])
			if lims is not None and n != 4:
				pylab.xlim(-lims[x], lims[x])
				pylab.ylim(-lims[y], lims[y])
			#pylab.title(title)

		if fname is None:
			pylab.show()
		else:
			pylab.savefig(fname, dpi=300)
			pylab.clf()

	def get_emittance(self):
		"return emittance h and v in m rad. Uses the bunch full width, so should only be used for a hard edge distribution"
		self.check_bunch()
		centers = self.get_centers()
		Ys = self.coords['Y'] - centers[0] # work relative to center
		Ts = self.coords['T'] - centers[1]
		Zs = self.coords['Z'] - centers[2]
		Ps = self.coords['P'] - centers[3]

		r_yt = numpy.sqrt(Ys**2 + Ts**2)
		theta_yt = numpy.arctan2(Ys, Ts)
		major_angle = theta_yt[r_yt.argmax()]
		theta_yt -= major_angle
		rot_y = r_yt * numpy.sin(theta_yt)
		rot_t = r_yt * numpy.cos(theta_yt)
		emittance_h = rot_y.max()*rot_t.max()

		
		r_zp = numpy.sqrt(Zs**2 + Ps**2)
		theta_zp = numpy.arctan2(Zs, Ps)
		major_angle = theta_zp[r_zp.argmax()]
		theta_zp -= major_angle
		rot_z = r_zp * numpy.sin(theta_zp)
		rot_p = r_zp * numpy.cos(theta_zp)
		emittance_v = rot_z.max() * rot_p.max()
		#print "Emittance (h, v):", emittance_h, emittance_v
		return (emittance_h, emittance_v)

	def get_emittance_rms(self):
		"return emittance h and v in m rad. Uses the RMS quantities"
		self.check_bunch()
		centers = self.get_centers()
		Ys = self.coords['Y'] - centers[0] # work relative to center
		Ts = self.coords['T'] - centers[1]
		Zs = self.coords['Z'] - centers[2]
		Ps = self.coords['P'] - centers[3]

		mxs = (Ys**2).mean()
		mxps = (Ts**2).mean()
		mxxp = (Ys*Ts).mean()
		e_h_rms = sqrt(mxs * mxps - mxxp**2)

		mys = (Zs**2).mean()
		myps = (Ps**2).mean()
		myyp = (Zs*Ps).mean()
		e_v_rms = sqrt(mys * myps - myyp**2)
		return (e_h_rms, e_v_rms)

	def get_twiss(self, emittance):
		"Returns the twiss valuse Beta_h, Alpha_h, Beta_v, Alpha_v, calculated from bunch width extent"
		self.check_bunch()
		# emittance may be a single number, or tuple (emittance_h, emittance_v)
		try:
			emittance_h, emittance_v = emittance
		except TypeError:
			emittance_h = emittance_v = emittance
		widths = self.get_widths()
		centers = self.get_centers()
		Ys = self.coords['Y'] - centers[0] # work relative to center
		Ts = self.coords['T'] - centers[1]
		Zs = self.coords['Z'] - centers[2]
		Ps = self.coords['P'] - centers[3]
		beta_h = (widths[0] / 2)**2 / emittance_h
		beta_v = (widths[2] / 2)**2 / emittance_v
		#print "beta", beta_h, beta_v
		#gamma_h = (widths[1] / 2)**2 / emittance_h
		#gamma_v = (widths[3] / 2)**2 / emittance_v
		# get T of particle with bigest Y
		y_p_max = Ts[Ys.argmax()]
		y_p_min = Ts[Ys.argmin()]
		alpha_h = (y_p_min - y_p_max) / 2 / sqrt(emittance_h / beta_h)

		z_p_max = Ps[Zs.argmax()]
		z_p_min = Ps[Zs.argmin()]
		alpha_v = (z_p_min - z_p_max) / 2 / sqrt(emittance_v / beta_v)
		#print "gamma", gamma_h, gamma_v	

		# following is dangerous, due do stat fluctuations causing roots of negative numbers
		#	alpha_h = sqrt((gamma_h*beta_h)-1 )
		#	alpha_v = sqrt((gamma_v*beta_v)-1 )

		#print "twiss (bh,ah,bv,av):", beta_h, alpha_h, beta_v, alpha_v
		return beta_h, alpha_h, beta_v, alpha_v


	def get_twiss_rms(self, emittance):
		"Returns the rms twiss valuse Beta_h, Alpha_h, Beta_v, Alpha_v, calculated from bunch rms width"
		self.check_bunch()
		# emittance may be a single number, or tuple (emittance_h, emittance_v)
		try:
			emittance_h, emittance_v = emittance
		except TypeError:
			emittance_h = emittance_v = emittance

		centers = self.get_centers()
		Ys = self.coords['Y'] - centers[0] # work relative to center
		Ts = self.coords['T'] - centers[1]
		Zs = self.coords['Z'] - centers[2]
		Ps = self.coords['P'] - centers[3]

		beta_h = (Ys**2).mean()/emittance_h
		alpha_h = -(Ys*Ts).mean()/emittance_h
		beta_v = (Zs**2).mean()/emittance_v
		alpha_v = -(Zs*Ps).mean()/emittance_v

		return beta_h, alpha_h, beta_v, alpha_v


