#!/usr/bin/env python

"A bunch object to hold the coordinates for many particles"

from math import *
import numpy
import pylab
from zgoubi import rel_conv

#from zgoubi.utils import *


class Bunch(object):
	"""Object to store a bunch of particles efficiently using numpy.
	All values are in SI units, m, rad, eV, s

	"""
	data_def = [
	('D', numpy.float64), # these coorspond to zgoubi D,Y,T,Z,P,S, but in SI units
	('Y', numpy.float64),
	('T', numpy.float64),
	('Z', numpy.float64),
	('P', numpy.float64),
	('S', numpy.float64),
	('tof', numpy.float64), # these are for accumulating 
	('X', numpy.float64),
	]
	def __init__(self, nparticles=0, ke=0, rigidity=0, mass=0, charge=1):
		self.coords = numpy.zeros(nparticles, self.data_def)
		self.mass = mass
		self.charge = charge
		self.rigidity = rigidity 
		if ke != 0:
			self.set_bunch_ke(ke)
		self.coords['D'] = 1

	
	def set_bunch_ke(self, ke):
		"Set bunch kinetic energy"
		if self.mass == 0:
			raise ValueError, "Particle mass can't be Zero"
		if self.charge == 0:
			raise ValueError, "Particle charge can't be Zero"
		self.rigidity = rel_conv.ke_to_rigidity(mass=self.mass, ke=ke, charge=self.charge)

	def get_bunch_ke(self):
		"Get bunch kinetic energy"
		if self.mass == 0:
			raise ValueError, "Particle mass can't be Zero"
		if self.charge == 0:
			raise ValueError, "Particle charge can't be Zero"
		return rel_conv.rigidity_to_ke(mass=self.mass, rigidity=self.rigidity, charge=self.charge)

	def set_bunch_rigidity(self, rigidity):
		"Set bunch rigidity"
		self.rigidity = rigidity

	def get_bunch_rigidity(self):
		"Get bunch rigidity"
		return self.rigidity

	def particles(self):
		"Returns the numpy array that holds the coordinates"
		return self.coords

	def get_min_BORO(self):
		"Returns the minimum rigidity of the bunch"
		#min_BORO = ke_to_rigidity(min(self.coords['KE']), self.mass)
		min_BORO = self.rigidity * self.coords['D'].min()
		return min_BORO

	def gen_halo_x_xp_y_yp(self, npart, emit_y, emit_z, beta_y, beta_z, alpha_y, alpha_z, seed=None):
		"Generate a halo bunch, i.e. an elipse in x-xp (Y-T) and in y-yp (Z-P)"

		#r = numpy.random.random_sample([npart])
		ry = sqrt(emit_y) 
		rz = sqrt(emit_z) 

		if seed != None:
			numpy.random.seed(seed)

		u1 = numpy.random.random_sample([npart]) * pi * 2
		u2 = numpy.random.random_sample([npart]) * pi * 2
		#u3 = numpy.random.random_sample([npart]) * pi * 2
		#u4 = numpy.random.random_sample([npart]) * pi * 2

		coords = numpy.zeros([npart, 6], numpy.float64)
		coords2 = numpy.zeros([npart, 6], numpy.float64)

		coords[:, 0] = ry * numpy.cos(u1)
		coords[:, 1] = ry * numpy.sin(u1)
		coords[:, 2] = rz * numpy.cos(u2)
		coords[:, 3] = rz * numpy.sin(u2)

		matrix = self._twiss_matrix(beta_y, beta_z, alpha_y, alpha_z)
		
		for n, coord in enumerate(coords):
			#	new_coord = numpy.dot(coord, matrix)
			new_coord = numpy.dot(matrix, coord)
			coords2[n] = new_coord 
		
		bunch =  numpy.zeros([npart], self.data_def)
		
		bunch['Y'] = coords2[:, 0]
		bunch['T'] = coords2[:, 1]
		bunch['Z'] = coords2[:, 2]
		bunch['P'] = coords2[:, 3]
		
	#	bunch['Y'] = coords[:,0]
	#	bunch['T'] = coords[:,1]
	#	bunch['Z'] = coords[:,2]
	#	bunch['P'] = coords[:,3]
		self.coords = bunch

	def _twiss_matrix(self, beta_y, beta_z, alpha_y, alpha_z):
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
		return M


	def read_YTZPSD(self, fname, ke, lim=None):
		"Read in a bunch from a file. assumes columns are Y, T, Z, P, X, D (though D is ignored), separated by white space"
		dist = numpy.loadtxt(fname)
		nparts = dist.size / 6
		dist = dist.reshape(nparts, 6)
		if lim:
			dist = dist[:lim]
		self.coords = numpy.zeros(nparts, self.data_def)
		#self.coords['KE'] = ke
		self.coords['Y'] = dist[:, 0]
		self.coords['T'] = dist[:, 1]
		self.coords['Z'] = dist[:, 2]
		self.coords['P'] = dist[:, 3]
		self.coords['X'] = dist[:, 4]
		self.coords['D'] = dist[:, 5]

	def write_YTZPSD(self, fname):
		"Output a bunch, compatible with read_YTZPSD"
		nparts = len(self.coords)
		dist = numpy.zeros([nparts, 6])
		#FIXME should set D from KE
		dist[:, 0] = self.coords['Y']
		dist[:, 1] = self.coords['T']
		dist[:, 2] = self.coords['Z']
		dist[:, 3] = self.coords['P']
		dist[:, 4] = self.coords['X']
		dist[:,5] = self.coords['D']
		#dist[:, 5] = 1
		#dist = dist.reshape(nparts * 2, 3)
		fh = open(fname, "w")
		fh.write("# bunch\n\n\n\n")
		numpy.savetxt(fh, dist)

	
	def get_widths(self):
		"Returns the width of the bunch in each dimension Y,T,Z,P,X,D (D not calculated yet)"
		y_width = max(self.coords['Y']) - min(self.coords['Y'])
		t_width = max(self.coords['T']) - min(self.coords['T'])
		z_width = max(self.coords['Z']) - min(self.coords['Z'])
		p_width = max(self.coords['P']) - min(self.coords['P'])
		x_width = max(self.coords['X']) - min(self.coords['X'])
		#d_width = max(self.coords['D']) - min(self.coords['D'])
		# FIXME should set from KE
		d_width = 0
		return (y_width, t_width, z_width, p_width, x_width, d_width)

	def __len__(self):
		"Returns length of bunch. Use len(my_bunch)"
		return len(self.coords)

	def plot(self, fname=None, lims=None):
		"Plot a bunch, if no file name give plot is displayed on screen."
		pylab.subplot(2, 2, 1)
		pylab.grid()
		pylab.plot(self.coords['Y'], self.coords['Z'], ',')
		if lims != None:
			pylab.xlim(-lims[0], lims[0])
			pylab.ylim(-lims[2], lims[2])
		plotname = "X-Y (Y-T)"
		pylab.title(plotname)

		
		pylab.subplot(2, 2, 2)
		pylab.grid()
		pylab.plot(self.coords['Y'], self.coords['T'], ',')
		if lims != None:
			pylab.xlim(-lims[0], lims[0])
			pylab.ylim(-lims[1], lims[1])
		plotname = "X-XP (Y-T"
		pylab.title(plotname)

		pylab.subplot(2, 2, 3)
		pylab.grid()
		pylab.plot(self.coords['Z'], self.coords['P'], ',')
		if lims != None:
			pylab.xlim(-lims[2], lims[2])
			pylab.ylim(-lims[3], lims[3])
		plotname = "Y-YP (z-P)"
		pylab.title(plotname)

		pylab.subplot(2, 2, 4)
		pylab.grid()
		pylab.plot(self.coords['X'], self.coords['D'], ',')
		if False:# lims != None:
			pylab.xlim(-lims[4], lims[4])
			pylab.ylim(-lims[5], lims[5])
		plotname = "X-D"
		pylab.title(plotname)
		if fname == None:
			pylab.show()
		else:
			pylab.savefig(fname, dpi=300)
			pylab.clf()
	
	def get_emmitance(self):
		"return emittance h and v in m rad"
		r_yt = numpy.sqrt(self.coords['Y']**2 + self.coords['T']**2)
		theta_yt = numpy.arctan2(self.coords['Y'], self.coords['T'])
		major_angle = theta_yt[r_yt.argmax()]
		theta_yt = theta_yt - major_angle
		rot_y = r_yt * numpy.sin(theta_yt)
		rot_t = r_yt * numpy.cos(theta_yt)
		emmitance_h = rot_y.max()*rot_t.max()

		r_zp = numpy.sqrt(self.coords['Z']**2 + self.coords['P']**2)
		theta_zp = numpy.arctan2(self.coords['Z'], self.coords['P'])
		major_angle = theta_zp[r_zp.argmax()]
		theta_zp = theta_zp - major_angle
		rot_z = r_zp * numpy.sin(theta_zp)
		rot_p = r_zp * numpy.cos(theta_zp)
		emmitance_v = rot_z.max() * rot_p.max()
		#print "Emmitance (h, v):", emmitance_h, emmitance_v
		return (emmitance_h, emmitance_v)


	def get_twiss(self, emittance):
		"Returns the twiss valuse Beta_h, Alpha_h, Beta_v, Alpha_v, calculated from bunch width extent"
		# emittance may be a single number, or tuple (emittance_h, emittance_v)
		try:
			emittance_h, emittance_v = emittance
		except TypeError:
			emittance_h = emittance_v = emittance
		widths = self.get_widths()
		beta_h = (widths[0] / 2)**2 / emittance_h
		beta_v = (widths[2] / 2)**2 / emittance_v
		#print "beta", beta_h, beta_v
		#gamma_h = (widths[1] / 2)**2 / emittance_h
		#gamma_v = (widths[3] / 2)**2 / emittance_v
		# get T of particle with bigest Y
		# FIXME, maybe can also look at T of particle with smallest Y, and average
		y_p = self.coords['T'][self.coords['Y'].argmax()]
		alpha_h = - y_p / sqrt(emittance_h / beta_h)
		z_p = self.coords['P'][self.coords['Z'].argmax()]
		alpha_v = - z_p / sqrt(emittance_v / beta_v)
		#print "gamma", gamma_h, gamma_v	

		# following is dangerous, due do stat fluctuations causing roots of negative numbers
		#	alpha_h = sqrt((gamma_h*beta_h)-1 )
		#	alpha_v = sqrt((gamma_v*beta_v)-1 )

		#print "twiss (bh,ah,bv,av):", beta_h, alpha_h, beta_v, alpha_v
		return beta_h, alpha_h, beta_v, alpha_v


