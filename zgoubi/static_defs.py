# -*- coding: utf-8 -*-

from __future__ import division, print_function
import os.path
from zgoubi.elements import zgoubi_element, PARTICUL
from zgoubi.constants import *
import copy

__all__ = ["OBJET1", "OBJET2", "OBJET3", "OBJET_bunch", "OBJET5",
"MCOBJET3", "zgoubi_particul", "ELECTRON", "PROTON", "MUON",
"IMMORTAL_MUON", "IMMORTAL_PION", "CHANGREF_NEW", "SPNTRK", "FAKE_ELEM"]

nl="\n"

class OBJET1(zgoubi_element):
	"""Beam made of grid of particles, default params give just 1 reference particle.
	Equivilent to OBJET with a KOBJ=1
	"""
	def __init__(self ,**settings):
		self._zgoubi_name = "OBJET"
		self._class_name = "OBJET1"
		self.label1 = ""
		self.label2 = ""
		self._params = {}
		self._params['BORO'] = 0
		for p in ['Y','T','Z','P','X','D']:
			self._params[str('I'+p)] = 1
			self._params[str('P'+p)] = 0
			self._params[str(p+'R')] = 0 # ref should be 0
		self._params[str('DR')] = 1 # appart from dp/p
		object.__setattr__(self, "ready", True)
		self.set(settings)

	def output(self):
		assert(self.IY * self.IT * self.IZ * self.IP * self.IX * self.ID < 10000)
		assert(max(self.IY, self.IT, self.IZ, self.IP, self.IX,self.ID) < 41 )
		
		# local short cuts for long function names
		f = self.f2s
		i = self.i2s
		
		out = "'OBJET'" +nl
		out += f(self.BORO) +nl
		out += "1" + nl
		out += i(self.IY) +' '+ i(self.IT) +' '+ i(self.IZ) +' '+ i(self.IP) +' '+ i(self.IX) +' '+ i(self.ID) +nl
		out += f(self.PY) +' '+ f(self.PT) +' '+ f(self.PZ) +' '+ f(self.PP) +' '+ f(self.PX) +' '+ f(self.PD) +nl
		out += f(self.YR) +' '+ f(self.TR) +' '+ f(self.ZR) +' '+ f(self.PR) +' '+ f(self.XR) +' '+ f(self.DR) +nl
		return out


class OBJET2(zgoubi_element):
	"""Beam made of particles, with coords give explicitly
	Equivilent to OBJET with a KOBJ=2
	"""
	def __init__(self ,**settings):
		self._zgoubi_name = "OBJET"
		self._class_name = "OBJET2"
		self.label1 = ""
		self.label2 = ""
		self._params = {}
		self._params['BORO'] = 0
		self._params['particles'] = []
		self._params['IMAX'] = 0
		self._params['IDMAX'] = 0
		object.__setattr__(self, "ready", True)
		self.set(settings)
		self.sorted = False
		
	def add(self, **settings):
		"add a particle"
		particle_coords = dict(Y=0, T=0, Z=0, P=0, X=0, D=1, LET=' ')
		for k, v in settings.items():
			particle_coords[k] = v
			
		self._params['particles'].append(particle_coords)
		# keep list sorted by D, needed for output()
		if self.sorted:
			self._params['particles'].sort(key=itemgetter('D'))
		
		
	def clear(self):
		"remove all particles"
		self._params['particles'] = []
	
	def output(self):
		out=''
		self.IMAX = len(self.particles)
		assert(self.IMAX <= 10000)
		assert(self.IMAX >= 1)
		
		#count unique 'D' values
		if self.sorted:
			self.IDMAX = len(set([x['D'] for x in self.particles]))
		else:
			self.IDMAX = self.IMAX
		
		# local short cuts for long function names
		f = self.f2s
		i = self.i2s
		
		out += "'OBJET'" +nl
		out += f(self.BORO) +nl
		out += "2" + nl
		out += i(self.IMAX) +' '+ i(self.IDMAX) + nl
		
		for part in self.particles:
			out += f(part['Y']) +' '+ f(part['T']) +' '+ f(part['Z']) +' '
			out += f(part['P']) +' '+ f(part['X']) +' '+ f(part['D'])
			#if (part['LET'] != ''):
			out += " '"+ part['LET'] + "'"
			out += nl
		
		# assume that we want to track all particles
		for x in range(len(self.particles)):
			out += i(1) + ' '
			if ((x+1)%10 == 0):
				out += nl # add a new line after 10 values
		
		#print out
		return out


class OBJET3(zgoubi_element):
	"""Read beam coordinates from file. KOBJ=3 read from formatted file.
		KOBJ=3.01 read from simple text file
	"""
	def __init__(self ,**settings):
		self._zgoubi_name = "OBJET"
		self._class_name = "OBJET3"
		self._params = {}
		self.label1 = ""
		self.label2 = ""
		self._params['BORO'] = 0
		self._params['IT1'] = 1
		self._params['IT2'] = 1
		self._params['ITStep'] = 1
		self._params['IP1'] = 1
		self._params['IP2'] = 1
		self._params['IPStep'] = 1
		self._params['TAG'] = '*'
		self._params['FNAME'] = 'nofilename'
		self._params['FTYPE'] = 'unspecified'
		for p in ['Y','T','Z','P','X','D','Ti']:
			self._params[str(p+'F')] = 1.0
			self._params[str(p+'R')] = 0.0
		self._params['InitC'] = 0
		object.__setattr__(self, "ready", True)
		self.set(settings)

	def output(self):
		# local short cuts for long function names
		f = self.f2s
		i = self.i2s
		
		out = "'OBJET'" +nl
		out += f(self.BORO) +nl
		if self.FTYPE == 'unformatted':
			out += "3.01" + " HEADER_4"+ nl
		elif self.FTYPE == 'formatted':
			out += "3" + nl
		else:
			print("Error - specify FTYPE, formatted or unformatted")
			sys.exit(1)

		out += i(self.IT1) +' '+ i(self.IT2) +' '+ i(self.ITStep) +nl
		out += i(self.IP1) +' '+ i(self.IP2) +' '+ i(self.IPStep) +nl
		out += f(self.YF)+' '+f(self.TF)+' '+f(self.ZF)+' '+ f(self.PF) +' '+ f(self.XF)+' '+f(self.DF)+' '+ f(self.TiF)+' '+ self.TAG +nl
		out += f(self.YR) +' '+ f(self.TR) +' '+ f(self.ZR) +' '+ f(self.PR) +' '+ f(self.XR) +' '+ f(self.DR)+' '+f(self.TiR)+nl
		out += i(self.InitC) +nl
		out += self.FNAME +nl
		return out

class OBJET_bunch(zgoubi_element):
	def __init__(self,bunch=None, binary=False,**settings):
		self._zgoubi_name = "OBJET"
		self._class_name = "OBJET_bunch"
		self._params = {}
		self.label1 = ""
		self.label2 = ""
		self.bunch = bunch
		self.binary = binary

	def setup(self, rundir):
		if self.binary:
			self.bunch.write_YTZPSD(os.path.join(rundir, "b_coords.dat"), binary=True)
		else:
			self.bunch.write_YTZPSD(os.path.join(rundir, "coords.dat"), binary=False)

	def output(self):
		if self.bunch is None:
			raise BadLineError("OBJET_bunch has no bunch set")

		f = self.f2s
		i = self.i2s
		out = "'OBJET'" +nl
		out += f(self.bunch.get_bunch_rigidity() *1000) + nl # convert from T.m to kG.cm
		out += "3.01" + " HEADER_4"+ nl
		out += "1 " + i(len(self.bunch)) + " 1" + nl
		out += "1 1 1" + nl
		out += "100 1000 100 1000 100 1 1 *" + nl
		out += "0 0 0 0 0 0 0" + nl
		out += "0" + nl
		if self.binary:
			out += "b_coords.dat" + nl
		else:
			out += "coords.dat" + nl

		return out



class OBJET5(zgoubi_element):
	"""Beam made of particles, for use with matrix
	Equivilent to OBJET with a KOBJ=5
	"""
	def __init__(self ,**settings):
		self._zgoubi_name = "OBJET"
		self._class_name = "OBJET5"
		self._params = {}
		self.label1 = ""
		self.label2 = ""
		self._params['BORO'] = 0
		self._params['ellipses'] = None
		for p in ['Y','T','Z','P','X','D']:
			self._params[str('P'+p)] = 1
			self._params[str(p+'R')] = 0
		object.__setattr__(self, "ready", True)
		self.set(settings)

	def add_ellipse(self, **settings):
		"Add a beam ellipse"
		ellipse_twiss = dict(alpha_y=0, beta_y=0, alpha_z=0 ,beta_z=0, alpha_s=0, beta_s=0,
		                     disp_y=0, disp_py=0, disp_z=0, disp_pz=0)
		for k, v in settings.items():
			ellipse_twiss[k] = v
			
		self._params['ellipses'] = ellipse_twiss
		
	def clear_ellipse(self):
		"remove all particles"
		self._params['ellipses'] = None

	def output(self):
		# local short cuts for long function names
		f = self.f2s
		i = self.i2s
		kobj = 5

		if self._params['ellipses']:
			kobj = 5.01
		
		out = "'OBJET'" +nl
		out += f(self.BORO) +nl
		out += str(kobj) + nl
		out += f(self.PY) +' '+ f(self.PT) +' '+ f(self.PZ) +' '+ f(self.PP) +' '+ f(self.PX) +' '+ f(self.PD) +nl
		out += f(self.YR) +' '+ f(self.TR) +' '+ f(self.ZR) +' '+ f(self.PR) +' '+ f(self.XR) +' '+ f(self.DR) +nl
		if self._params['ellipses']:
			e = self._params['ellipses']
			out += f(e["alpha_y"]) +' '+ f(e["beta_y"]) +' '+f(e["alpha_z"]) +' '+f(e["beta_z"]) +' '+f(e["alpha_s"]) +' '+f(e["beta_s"]) +' '+f(e["disp_y"]) +' '+f(e["disp_py"]) +' '+f(e["disp_z"]) +' '+f(e["disp_pz"])
		return out


class MCOBJET3(zgoubi_element):
	"""Monte Carlo generation of a 6-D object on a phase space ellipse
	Equivilent to MCOBJET with a KOBJ=3
	"""
	def __init__(self ,**settings):
		self._zgoubi_name = "MCOBJET"
		self._class_name = "MCOBJET3"
		self._params = {}
		self.label1 = ""
		self.label2 = ""
		self._params['BORO'] = 0
		self._params['IMAX'] = 0
		self._params['particles'] = []
		for p in ['Y','T','Z','P','X','D']:
			self._params[str('K'+p)] = 1
			self._params[str(p+'0')] = 0
		for p in ['y','z','x']:
			self._params[str('alpha_'+p)] = 0
			self._params[str('beta_'+p)] = 0
			self._params[str('emit_'+p)] = 0
			self._params[str('n_cutoff_'+p)] = 0
			self._params[str('n_cutoff2_'+p)] = 0
		self._params[str('I1')] = 139339
		self._params[str('I2')] = 139397
		self._params[str('I3')] = 178393
		object.__setattr__(self, "ready", True)
		self.set(settings)

	def output(self):
		# local short cuts for long function names
		f = self.f2s
		i = self.i2s
		
		out = "'MCOBJET'" +nl
		out += f(self.BORO) +nl
		out += "3" + nl
		out += i(self.IMAX) + nl
		out += i(self.KY) +' '+ i(self.KT) +' '+ i(self.KZ) +' '+ i(self.KP) +' '+ i(self.KX) +' '+ i(self.KD) +nl
		out += f(self.Y0) +' '+ f(self.T0) +' '+ f(self.Z0) +' '+ f(self.P0) +' '+ f(self.X0) +' '+ f(self.D0) +nl

		if self.n_cutoff_y < 0:
			out += f(self.alpha_y) +' '+ f(self.beta_y) +' '+ f(self.emit_y) +' '+ i(self.n_cutoff_y) +' '+i(self.n_cutoff2_y) +nl
		else:
			out += f(self.alpha_y) +' '+ f(self.beta_y) +' '+ f(self.emit_y) +' '+ i(self.n_cutoff_y) +nl

		if self.n_cutoff_z < 0:
			out += f(self.alpha_z) +' '+ f(self.beta_z) +' '+ f(self.emit_z) +' '+ i(self.n_cutoff_z) +' '+i(self.n_cutoff2_z) +nl
		else:
			out += f(self.alpha_z) +' '+ f(self.beta_z) +' '+ f(self.emit_z) +' '+ i(self.n_cutoff_z) +nl

		if self.n_cutoff_x < 0:
			out += f(self.alpha_x) +' '+ f(self.beta_x) +' '+ f(self.emit_x) +' '+ i(self.n_cutoff_x) +' '+i(self.n_cutoff2_x) +nl
		else:
			out += f(self.alpha_x) +' '+ f(self.beta_x) +' '+ f(self.emit_x) +' '+ i(self.n_cutoff_x) +nl

		out += i(self.I1) +' '+ i(self.I2) +' '+ i(self.I3) +nl
		return out


# define some useful particles
# constants defined in zgoubi_constants.py

class zgoubi_particul(PARTICUL):
	def __neg__(self):
		"Return an anti-particle, by inverting charge"
		new_particul = copy.copy(self)
		self._params["Q"] *= -1
			
		return new_particul

class ELECTRON(zgoubi_particul):
	def __init__(self):
		super(self.__class__, self).__init__()
		self._class_name = "ELECTRON"
		self.set(M=ELECTRON_MASS/1e6)
		self.set(Q=ELECTRON_CHARGE)
		self.set(G=ELECTRON_ANOM_MAG_MOM)
		self.set(tau=ELECTRON_MEAN_LIFE)

class PROTON(zgoubi_particul):
	def __init__(self):
		super(self.__class__, self).__init__()
		self._class_name = "PROTON"
		self.set(M=PROTON_MASS/1e6)
		self.set(Q=PROTON_CHARGE)
		self.set(G=PROTON_ANOM_MAG_MOM)
		self.set(tau=PROTON_MEAN_LIFE)

class MUON(zgoubi_particul):
	def __init__(self):
		super(self.__class__, self).__init__()
		self._class_name = "MUON"
		self.set(M=MUON_MASS/1e6)
		self.set(Q=MUON_CHARGE)
		self.set(G=MUON_ANOM_MAG_MOM)
		self.set(tau=MUON_MEAN_LIFE)

class IMMORTAL_MUON(zgoubi_particul):
	def __init__(self):
		super(self.__class__, self).__init__()
		self._class_name = "IMMORTAL_MUON"
		self.set(M=MUON_MASS/1e6)
		self.set(Q=MUON_CHARGE)
		self.set(G=MUON_ANOM_MAG_MOM)
		self.set(tau=0)

class IMMORTAL_PION(zgoubi_particul):
	def __init__(self):
		super(self.__class__, self).__init__()
		self._class_name = "IMMORTAL_PION"
		self.set(M=PION_MASS/1e6)
		self.set(Q=PION_CHARGE)
		self.set(G=0)
		self.set(tau=0)


class CHANGREF_NEW(zgoubi_element):
	""" updated CHANGREF works with format XS,YS,ZS for longitudinal, horizontal and vertical shifts
		and XR,YR,ZR for rotations about the longitudinal, horizontal and vertical axes
		order can be specified, for example order=["ZR","YS"].
	"""
	def __init__(self ,label1=None,**settings):
		self._zgoubi_name = "CHANGREF"
		self._class_name = "CHANGREF_NEW"
		self._params = {}
		self.label1 = ""
		self.label2 = ""
		if label1 is not None:
			self.label1 = label1
		self._params['XS'] = None
		self._params['YS'] = None
		self._params['ZS'] = None
		self._params['XR'] = None
		self._params['YR'] = None
		self._params['ZR'] = None
		#default order of transformations
		self._params['order'] = ["XS","YS","ZS","XR","YR","ZR"]
		object.__setattr__(self, "ready", True)
		self.set(settings)

	def output(self):
		# local short cuts for long function names
		f = self.f2s
		i = self.i2s

		out = "'CHANGREF'" +' '+self.label1+nl
		for elem in self.order:
			if elem == "XS" and self.XS is not None:
				out += "XS " + f(self.XS)+' '
			if elem == "YS" and self.YS is not None:
				out += "YS " + f(self.YS)+' '
			if elem == "ZS" and self.ZS is not None:
				out += "ZS " + f(self.ZS)+' '
			if elem == "XR" and self.XR is not None:
				out += "XR " + f(self.XR)+' '
			if elem == "YR" and self.YR is not None:
				out += "YR " + f(self.YR)+' '
			if elem == "ZR" and self.ZR is not None:
				out += "ZR " + f(self.ZR)+' '
		out += '\n'

		return out


class SPNTRK(zgoubi_element):
	""" SPNTRK assigns the initial spin. Include KSO options 1-3 (longitudinal, horizontal and vertical spin respectively) and KSO=4  which list the spin of each particle, where the number of particles should be equal to that defined in OBJET. If KSO=4, the initial spin should be assigned using the "spin_vector" list

		spin_vector = [[SX,SY,SZ],...] and so on for each particle in OBJET

    KSO=5 is not suppoerted by current Zgoubi and so is not included here.
	"""
	def __init__(self ,**settings):
		self._zgoubi_name = "SPNTRK"
		self._class_name = "SPNTRK"
		self._params = {}
		self.label1 = ""
		self.label2 = ""
		self._params['KSO'] = None
		self._params['spin_vector'] = []
		object.__setattr__(self, "ready", True)
		self.set(settings)

	def output(self):
		# local short cuts for long function names
		f = self.f2s
		i = self.i2s

		out = "'SPNTRK'" +nl
		out += i(self.KSO) + nl
		if self.KSO == 4:
			for spin in self.spin_vector:
				out += f(spin[0])+' '+f(spin[1])+' '+f(spin[2])+nl

		out += '\n'

		return out

		
class FAKE_ELEM(zgoubi_element):
	def __init__(self, data=""):
		self.data = data
		self.label1 = ""
		self.label2 = ""
		self._params = {}
		self._zgoubi_name = "####"
		self._class_name = "FAKE_ELEM"

	def output(self):
		return self.data+'\n'


