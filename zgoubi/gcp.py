from __future__ import division, print_function
import numpy
from zgoubi.core import *
from zgoubi.constants import *
from zgoubi.common import *
from zgoubi.utils import *
from zgoubi.exceptions import *

"""This module contains high level functions for analysing periodic cells and rings

Primary functions are
get_cell_properties()
get_cell_tracks()

"""

data_def = [
('KE',numpy.float64),
('stable',numpy.bool), # has closed orbit and stable tunes around it
('found_co',numpy.bool), # has closed orbit (may be unstable vertically)
('stable_tm_YT',numpy.bool), # has a stable transfer matrix in YT
('stable_tm_ZP',numpy.bool), # has a stable transfer matrix in ZP
('Y',numpy.float64),
('T',numpy.float64),
('Z',numpy.float64),
('P',numpy.float64),
('BETA_Y',numpy.float64),
('BETA_Z',numpy.float64),
('ALPHA_Y',numpy.float64),
('ALPHA_Z',numpy.float64),
('GAMMA_Y',numpy.float64),
('GAMMA_Z',numpy.float64),
('DISP_Y',numpy.float64),
('DISP_Z',numpy.float64),
('DISP_PY',numpy.float64),
('DISP_PZ',numpy.float64),
('NU_Y',numpy.float64),
('NU_Z',numpy.float64),
('matrix', numpy.float64, (6,6)),
('matrix_trace_YT', numpy.float64),
('matrix_trace_ZP', numpy.float64),
('tof',numpy.float64),
('S',numpy.float64),
('MAX_BY', numpy.float64),
('MIN_BY', numpy.float64),
('MAX_BZ', numpy.float64),
('MIN_BZ', numpy.float64),
('twiss_profile', numpy.object),
('full_twiss_profile', numpy.object),
('DA', numpy.object),
('DA_angles', numpy.object),
('ftrack', numpy.object),
('ptrack', numpy.object),
('phase_space', numpy.object),
]

data_def_nonperiodic = [
('KE',numpy.float64),
('stable',numpy.bool), # some functions skip over "unstable" orbits
('Y',numpy.float64),
('T',numpy.float64),
('Z',numpy.float64),
('P',numpy.float64),
('BETA_Y',numpy.float64),
('BETA_Z',numpy.float64),
('ALPHA_Y',numpy.float64),
('ALPHA_Z',numpy.float64),
('GAMMA_Y',numpy.float64),
('GAMMA_Z',numpy.float64),
('DISP_Y',numpy.float64),
('DISP_Z',numpy.float64),
('DISP_PY',numpy.float64),
('DISP_PZ',numpy.float64),
('Y0',numpy.float64),
('T0',numpy.float64),
('Z0',numpy.float64),
('P0',numpy.float64),
('BETA_Y0',numpy.float64),
('BETA_Z0',numpy.float64),
('ALPHA_Y0',numpy.float64),
('ALPHA_Z0',numpy.float64),
('GAMMA_Y0',numpy.float64),
('GAMMA_Z0',numpy.float64),
('DISP_Y0',numpy.float64),
('DISP_Z0',numpy.float64),
('DISP_PY0',numpy.float64),
('DISP_PZ0',numpy.float64),
('NU_Y',numpy.float64),
('NU_Z',numpy.float64),
('matrix', numpy.float64, (6,6)),
('tof',numpy.float64),
('S',numpy.float64),
('MAX_BY', numpy.float64),
('MIN_BY', numpy.float64),
('MAX_BZ', numpy.float64),
('MIN_BZ', numpy.float64),
('twiss_profile', numpy.object),
('full_twiss_profile', numpy.object),
('ftrack', numpy.object),
('ptrack', numpy.object),
]


class GCPData(numpy.ndarray):
	"Subclass of numpy.ndarray to add an info field"
	# based on example in https://docs.scipy.org/doc/numpy/user/basics.subclassing.html
	def __new__(subtype, shape, dtype=None, buffer=None, offset=0, strides=None, order=None, info=None):
		info2 = dict(periodic=True, particle="")
		info2.update(info)
		if dtype is None:
			if info2["periodic"]:
				dtype = data_def
			else:
				dtype = data_def_nonperiodic

		obj = numpy.ndarray.__new__(subtype, shape, dtype, buffer, offset, strides,order)
		obj.info = info2
		return obj

	def __array_finalize__(self, obj):
		if obj is None: return
		self.info = getattr(obj, 'info', None)

	@classmethod
	def from_ndarray(cls, data):
		if "Y0" in data.dtype.names:
			nc = cls(data.shape, info={"periodic":False})
		else:
			nc = cls(data.shape, info={"periodic":True})
		nc[:] = data[:]
		return nc



def part_info(particle):
	"Look up a particle by name and return the zgoubi PARTICUL, mass and charge sign"
	if isinstance(particle, PARTICUL):
		part_ob = particle
	elif particle == "p":
		part_ob = PROTON()
	elif particle == "e":
		part_ob = ELECTRON()
	elif particle == "mu-":
		part_ob = IMMORTAL_MUON()
	elif particle == "mu+":
		part_ob = IMMORTAL_MUON()
	elif particle == "pi+":
		part_ob = IMMORTAL_PION()
	elif particle == "pi-":
		part_ob = IMMORTAL_PION()
	else:
		raise ValueError("Unknown particle")
	mass = part_ob.M *1e6
	charge_sign = part_ob.Q / -ELECTRON_CHARGE
	return part_ob, mass, charge_sign


def get_cell_properties(cell, min_ke, max_ke=None, ke_steps=1, particle=None, tol=1e-6, stop_at_first_unstable=False, closed_orbit_range=None, closed_orbit_range_count=None, closed_orbit_init_YTZP=None, reuse_co_coords=True, closed_orbit_debug=False, full_tracking=False, smart_co_search=False):
	"""Get the closed orbits and basic properties of a periodic cell.

	cell: A PyZgoubi Line object containing the beamline elements
	min_ke, max_ke, ke_steps: kinetic energy in eV. For a single step just set min_ke
	particle: "p", "e", "mu-", "mu+" or a PARTICUL() instance
	tol: tolerance when finding closed orbit
	stop_at_first_unstable: If an energy is found to be unstable, give up
	closed_orbit_range, closed_orbit_range_count, closed_orbit_init_YTZ: see zgoubi.utils.find_closed_orbit_range()
	reuse_co_coords: use closed orbit from previous energy to start search for next energy
	closed_orbit_debug: output debuging information
	full_tracking=True is required in order get minimum and maximum magnetic fields along orbit

	returns orbit_data, an array with ke_steps elements, with the following data
	KE : particle KE in eV
	stable, found_co, stable_tm_YT, stable_tm_ZP : boolean stability flags
	Y, T, Z, P : closed orbit at start of cell, in cm and mrad
	BETA_Y, BETA_Z, ALPHA_Y, ALPHA_Z, GAMMA_Y, GAMMA_Z, DISP_Y, DISP_Z, DISP_PY, DISP_PZ, NU_Y, NU_Z: periodic twiss functions
	tof, S: time of flight and path length along closed orbit
	matrix, matrix_trace_YT, matrix_trace_ZP: transfer matrix and traces
	MAX_BY, MIN_BY, MAX_BZ, MIN_BZ: minimum and maximum fields seen along closed orbit (requires full_tracking=True)
	"""

	if max_ke is None: max_ke = min_ke
	ke_list = numpy.linspace(min_ke, max_ke, ke_steps)

	orbit_data = GCPData(ke_steps, info=dict(periodic=True, particle=particle))
	for n, particle_ke in enumerate(ke_list):
		orbit_data['KE'][n] = particle_ke

	# get closed orbits
	tline = Line('test_line')
	tline.add_input_files(cell.input_files)
	ob = OBJET2()
	tline.add(ob)
	part_ob, mass, charge_sign = part_info(particle)
	tline.add(part_ob)

	tline.add(DRIFT("fco", XL=0* cm_))

	tline.add(cell)

	tline.add(DRIFT("end", XL=0* cm_))
	tline.add(FAISCNL("end", FNAME='zgoubi.fai'))
	tline.add(REBELOTE(NPASS=20, K=99))
	tline.add(END())

	search_coords = [0,0,0,0]
	if closed_orbit_init_YTZP is not None:
		search_coords = closed_orbit_init_YTZP
	if closed_orbit_range is not None and closed_orbit_range_count is None:
		# if range but no count, then do 10 steps in dimensions that are used	
		closed_orbit_range_count = []
		for r in closed_orbit_range:
			if r == 0: closed_orbit_range_count.append(0)
			else: closed_orbit_range_count.append(10)

	for n, particle_ke in enumerate(ke_list):
		orbit_data['found_co'][n] = False
		orbit_data['stable_tm_YT'][n] = False
		orbit_data['stable_tm_ZP'][n] = False
		orbit_data['stable'][n] = False
		print("closed orbit, energy = ", particle_ke)
		rigidity = ke_to_rigidity(particle_ke,mass) / charge_sign
		ob.set(BORO=rigidity)
		good_cos = (orbit_data['found_co']*1).sum()
		if smart_co_search and good_cos >= 2:
			for coordn, coord in enumerate("Y"):
				good_data = orbit_data[orbit_data['found_co']]
				good_data = good_data[-5:]
				co_poly = numpy.polyfit(good_data['KE'], good_data[coord], min(good_cos,4)-1) # get linear or quad fit
				search_coords[coordn] = numpy.poly1d(co_poly)(particle_ke)

		if closed_orbit_debug:
			record_fname = "closedorbit_%s.log"%n
		else:
			record_fname = None

		if closed_orbit_range is None:
			closed_orbit =  find_closed_orbit(tline, init_YTZP=search_coords, tol=tol, max_iterations=50, record_fname=record_fname)
		else:
			closed_orbit =  find_closed_orbit_range(tline, range_YTZP=closed_orbit_range, count_YTZP=closed_orbit_range_count, init_YTZP=search_coords, tol=tol, max_iterations=50, record_fname=record_fname)
	
		if closed_orbit is not None:
			orbit_data['Y'][n],orbit_data['T'][n],orbit_data['Z'][n], orbit_data['P'][n]= closed_orbit
			orbit_data['found_co'][n] = True
			search_coords = closed_orbit
		else:
			zlog.warn("No closed orbit at: %s"% particle_ke)
			if closed_orbit_debug:
				plot_find_closed_orbit(data_fname=record_fname, outfile=record_fname+".pdf")
				print("Search plots writen to ", record_fname+".pdf")
			if stop_at_first_unstable: break

	# get tunes and twiss
	tline = Line('test_line')
	tline.add_input_files(cell.input_files)
	ob = OBJET5()
	tline.add(ob)
	tline.add(part_ob)
	tline.add(DRIFT("fco", XL=0* cm_))
	tline.add(FAISCNL(FNAME='zgoubi.fai',))

	tline.add(cell)

	tline.add(DRIFT("end", XL=0* cm_))
	tline.add(FAISCNL("end",FNAME='zgoubi.fai',))
	tline.add(MATRIX(IORD=1, IFOC=11))
	tline.add(END())

	for n, particle_ke in enumerate(ke_list):
		print("twiss, energy = ", particle_ke)
		
		rigidity = ke_to_rigidity(particle_ke,mass) / charge_sign
		ob.set(BORO=rigidity)
		
		ref_Y,ref_T,ref_Z,ref_P = orbit_data['Y'][n],orbit_data['T'][n],orbit_data['Z'][n], orbit_data['P'][n]
		#print "ref_Y,ref_T,ref_Z,ref_P", ref_Y,ref_T,ref_Z,ref_P
		ob.set(YR=ref_Y, TR=ref_T, ZR=ref_Z, PR=ref_P, XR=0, DR=1)
		step_disp = 0.001 #cm
		step_ang =  0.01 #mrad
		ob.set(PY=step_disp, PT=step_ang, PZ=step_disp, PP=step_ang, PX=step_disp, PD=0.001)
		tline.full_tracking(full_tracking)

		res =tline.run()
		try:
			orbit_data['matrix'][n] = res.get_transfer_matrix()
		except BadLineError:
			continue
		
		orbit_data['matrix_trace_YT'][n] = orbit_data['matrix'][n][0,0] + orbit_data['matrix'][n][1,1]
		orbit_data['matrix_trace_ZP'][n] = orbit_data['matrix'][n][2,2] + orbit_data['matrix'][n][3,3]
		
		orbit_data['stable_tm_YT'][n] = abs(orbit_data['matrix_trace_YT'][n]) < 2
		orbit_data['stable_tm_ZP'][n] = abs(orbit_data['matrix_trace_ZP'][n]) < 2

		if orbit_data['stable_tm_YT'][n] and orbit_data['stable_tm_ZP'][n] and orbit_data['found_co'][n]:
			orbit_data['stable'][n] = True
		
		tune = res.get_tune()
		orbit_data['NU_Y'][n],orbit_data['NU_Z'][n] = tune
		
		if tune[0] == 0 or tune[1] == 0 or isnan(tune[0]) or isnan(tune[1]):
			zlog.warn("tune is zero or nan")
			orbit_data['stable'][n] = False
			#tline.run(xterm=1)
		
		twiss = res.get_twiss_parameters()
		orbit_data['BETA_Y'][n] = twiss['beta_y']
		orbit_data['ALPHA_Y'][n] = twiss['alpha_y']
		orbit_data['GAMMA_Y'][n] = twiss['gamma_y']
		orbit_data['DISP_Y'][n] = twiss['disp_y']
		orbit_data['DISP_PY'][n] = twiss['disp_py']
		orbit_data['BETA_Z'][n] = twiss['beta_z']
		orbit_data['ALPHA_Z'][n] = twiss['alpha_z']
		orbit_data['GAMMA_Z'][n] = twiss['gamma_z']
		orbit_data['DISP_Z'][n] = twiss['disp_z']
		orbit_data['DISP_PZ'][n] = twiss['disp_pz']

		if not orbit_data['found_co'][n]: continue
		if full_tracking:
			ptrack = res.get_all('plt')
			by = ptrack['BY']
			bz = ptrack['BZ']
			orbit_data['MAX_BY'][n] = by.max()
			orbit_data['MIN_BY'][n] = by.min()
			orbit_data['MAX_BZ'][n] = bz.max()
			orbit_data['MIN_BZ'][n] = bz.min()

		ftrack = res.get_all('fai')
		ftrack = ftrack[ftrack['ID']==1]
		orbit_data['tof'][n] = ftrack['tof'][-1]
		orbit_data['S'][n] = ftrack['S'][-1]
		res.clean()

		if orbit_data['NU_Y'][n] == -1 or orbit_data['NU_Z'][n] == -1:
			orbit_data['stable'][n] = False
			if stop_at_first_unstable: break
		
	return orbit_data


def get_cell_properties_nonperiodic(cell, min_ke, max_ke=None, ke_steps=1, particle=None, init_YTZP=None, init_twiss=None, full_tracking=False):
	"""Get the basic properties of a non-periodic cell. 

	Works similarly to get_cell_properties(), but rather than finding a periodic solution for closed orbit and twiss parameters, takes them as input:
	init_YTZP: list of the Y, T, Z and P at the start of the cell
	init_twiss: initial twiss params, use twiss_param_array() to create

	cell: A PyZgoubi Line object containing the beamline elements
	min_ke, max_ke, ke_steps: kinetic energy in eV. For a single step just set min_ke
	particle: "p", "e", "mu-", "mu+", or a PARTICUL() instance
	full_tracking=True is required in order get minimum and maximum magnetic fields along orbit

	returns orbit_data, an array with ke_steps elements, with the following data
	KE : particle KE in eV
	stable : boolean stability flags, always true, but useful for controlling which orbits are processed by other functions
	Y0, T0, Z0, P0 : position and angles at start of cell, in cm and mrad
	Y, T, Z, P : position and angles at end of cell, in cm and mrad
	BETA_Y0, BETA_Z0, ALPHA_Y0, ALPHA_Z0, GAMMA_Y0, GAMMA_Z0, DISP_Y0, DISP_Z0, DISP_PY0, DISP_PZ0, NU_Y0, NU_Z0: twiss functions at start of cell
	BETA_Y, BETA_Z, ALPHA_Y, ALPHA_Z, GAMMA_Y, GAMMA_Z, DISP_Y, DISP_Z, DISP_PY, DISP_PZ, NU_Y, NU_Z: twiss functions at end of cell
	tof, S: time of flight and path length along closed orbit
	matrix, matrix_trace_YT, matrix_trace_ZP: transfer matrix and traces
	MAX_BY, MIN_BY, MAX_BZ, MIN_BZ: minimum and maximum fields seen along closed orbit (requires full_tracking=True)
	"""

	if max_ke is None: max_ke = min_ke
	if init_YTZP is None: init_YTZP=[0,0,0,0]
	ke_list = numpy.linspace(min_ke, max_ke, ke_steps)

	orbit_data = GCPData(ke_steps, info=dict(periodic=False, particle=particle))
	for n, particle_ke in enumerate(ke_list):
		orbit_data['KE'][n] = particle_ke
		orbit_data['Y0'][n],orbit_data['T0'][n],orbit_data['Z0'][n], orbit_data['P0'][n] = init_YTZP
	orbit_data["BETA_Y0"] = init_twiss["beta_y"]
	orbit_data["ALPHA_Y0"] = init_twiss["alpha_y"]
	orbit_data["GAMMA_Y0"] = init_twiss["gamma_y"]
	orbit_data["DISP_Y0"] = init_twiss["disp_y"]
	orbit_data["BETA_Z0"] = init_twiss["beta_z"]
	orbit_data["ALPHA_Z0"] = init_twiss["alpha_z"]
	orbit_data["GAMMA_Z0"] = init_twiss["gamma_z"]
	orbit_data["DISP_Z0"] = init_twiss["disp_z"]
	orbit_data["stable"] = True

	part_ob, mass, charge_sign = part_info(particle)

	# get tunes and twiss
	tline = Line('test_line')
	tline.add_input_files(cell.input_files)
	ob = OBJET5()
	tline.add(ob)
	tline.add(part_ob)
	tline.add(DRIFT("fco", XL=0* cm_))
	tline.add(FAISCNL(FNAME='zgoubi.fai',))

	tline.add(cell)

	tline.add(DRIFT("end", XL=0* cm_))
	tline.add(FAISCNL("end",FNAME='zgoubi.fai',))
	tline.add(MATRIX(IORD=1, IFOC=11))
	tline.add(END())

	for n, particle_ke in enumerate(ke_list):
		print("twiss, energy = ", particle_ke)
		
		rigidity = ke_to_rigidity(particle_ke,mass) / charge_sign
		ob.set(BORO=rigidity)
		
		ref_Y,ref_T,ref_Z,ref_P = orbit_data['Y0'][n],orbit_data['T0'][n],orbit_data['Z0'][n], orbit_data['P0'][n]
		#print "ref_Y,ref_T,ref_Z,ref_P", ref_Y,ref_T,ref_Z,ref_P
		ob.set(YR=ref_Y, TR=ref_T, ZR=ref_Z, PR=ref_P, XR=0, DR=1)
		step_disp = 0.001 #cm
		step_ang =  0.01 #mrad
		ob.set(PY=step_disp, PT=step_ang, PZ=step_disp, PP=step_ang, PX=step_disp, PD=0.001)
		tline.full_tracking(full_tracking)

		res =tline.run()
		try:
			orbit_data['matrix'][n] = res.get_transfer_matrix()
		except BadLineError:
			continue
		
		tune = res.get_tune()
		orbit_data['NU_Y'][n],orbit_data['NU_Z'][n] = tune

		if full_tracking:
			ptrack = res.get_all('plt')
			by = ptrack['BY']
			bz = ptrack['BZ']
			orbit_data['MAX_BY'][n] = by.max()
			orbit_data['MIN_BY'][n] = by.min()
			orbit_data['MAX_BZ'][n] = bz.max()
			orbit_data['MIN_BZ'][n] = bz.min()

		orbit_data['BETA_Y0'][n] = init_twiss['beta_y']
		orbit_data['ALPHA_Y0'][n] = init_twiss['alpha_y']
		orbit_data['GAMMA_Y0'][n] = init_twiss['gamma_y']
		orbit_data['DISP_Y0'][n] = init_twiss['disp_y']
		orbit_data['DISP_PY0'][n] = init_twiss['disp_py']
		orbit_data['BETA_Z0'][n] = init_twiss['beta_z']
		orbit_data['ALPHA_Z0'][n] = init_twiss['alpha_z']
		orbit_data['GAMMA_Z0'][n] = init_twiss['gamma_z']
		orbit_data['DISP_Z0'][n] = init_twiss['disp_z']
		orbit_data['DISP_PZ0'][n] = init_twiss['disp_pz']

		# FIXME: can probably just propergate these thought the matrix
		# see SY Lee pg 48
		twiss_profiles = get_twiss_profiles(tline,None, input_twiss_parameters=init_twiss, calc_dispersion=0)
		orbit_data['BETA_Y'][n] = twiss_profiles['beta_y'][-1]
		orbit_data['ALPHA_Y'][n] = twiss_profiles['alpha_y'][-1]
		orbit_data['GAMMA_Y'][n] = twiss_profiles['gamma_y'][-1]
		orbit_data['DISP_Y'][n] = twiss_profiles['disp_y'][-1]
		orbit_data['DISP_PY'][n] = twiss_profiles['disp_py'][-1]
		orbit_data['BETA_Z'][n] = twiss_profiles['beta_z'][-1]
		orbit_data['ALPHA_Z'][n] = twiss_profiles['alpha_z'][-1]
		orbit_data['GAMMA_Z'][n] = twiss_profiles['gamma_z'][-1]
		orbit_data['DISP_Z'][n] = twiss_profiles['disp_z'][-1]
		orbit_data['DISP_PZ'][n] = twiss_profiles['disp_pz'][-1]


		ftrack = res.get_all('fai')
		ftrack = ftrack[ftrack['ID']==1]
		for key in "tof S Y T Z P".split():
			orbit_data[key][n] = ftrack[key][-1]
		res.clean()
		
	return orbit_data


def get_cell_tracks(cell, data, particle, full_tracking=False, xterm=False, add_faiscnl=True):
	"""Get tracks along the closed orbit for values in data
	cell: periodic cell
	data: the data structure return from get_cell_properties()
	particle: see get_cell_properties()
	full_tracking: record all steps in the magnet
	add_faiscnl: insert a faiscnl (beam store) after each element

	tracks are added in the data structures ftrack and ptrack fields
	
	"""
	if not isinstance(data, GCPData):
		data = GCPData.from_ndarray(data)
	for n, particle_ke in enumerate(data['KE']):
		if data.info["periodic"]:
			ref_Y,ref_T,ref_Z,ref_P = data[n]['Y'], data[n]['T'], data[n]['Z'], data[n]['P']
		else:
			ref_Y,ref_T,ref_Z,ref_P = data[n]['Y0'], data[n]['T0'], data[n]['Z0'], data[n]['P0']

		tracks = get_tracks(cell=cell, start_YTZP=[ref_Y,ref_T,ref_Z,ref_P],
		                    particle=particle, ke=particle_ke, full_tracking=full_tracking, xterm=xterm, add_faiscnl=add_faiscnl)
		data[n]['ftrack'] = tracks['ftrack']
		data[n]['ptrack'] = tracks['ptrack']


def get_tracks(cell, start_YTZP, particle, ke, full_tracking=False, return_zgoubi_files=False, xterm=False, add_faiscnl=True):
	"""Run a particle through a cell from a give starting point and return track from fai and plt files.
	
	This is mostly used by other functions in this module, but can be useful for debugging a lattice
	
	add_faiscnl: insert a faiscnl (beam store) after each element
	"""
	split = 1
	tline = Line('test_line')
	tline.add_input_files(cell.input_files)
	ob = OBJET2()
	tline.add(ob)
	part_ob, mass, charge_sign = part_info(particle)
	tline.add(part_ob)
	
	tline.add(DRIFT('start', XL=0))

	for e in cell.elements():
		if split == 1:
			tline.add(e)
			if add_faiscnl and (hasattr(e, "XL") or hasattr(e, "AT")):
				tline.add(FAISCNL(FNAME='zgoubi.fai',))
		else:
			for es in split_element(e, split):
				tline.add(es)
				if add_faiscnl and (hasattr(e, "XL") or hasattr(e, "AT")):
					tline.add(FAISCNL(FNAME='zgoubi.fai',))

	tline.add(DRIFT('end', XL=0))
	tline.add(END())

	rigidity = ke_to_rigidity(ke,mass) / charge_sign
	ob.set(BORO=rigidity)

	ref_Y,ref_T,ref_Z,ref_P = start_YTZP
	ob.clear()
	ob.add(Y=ref_Y, T=ref_T, Z=ref_Z, P=ref_P, X=0, D=1)
	if full_tracking:
		tline.full_tracking(True, drift_to_multi=True)

	res = tline.run(xterm=xterm)
	try:
		ftrack = res.get_all('fai')
	except IOError:
		ftrack = None
	if full_tracking:
		try:
			ptrack = res.get_all('plt')
		except EmptyFileError:
			zlog.warn("Empty plt file")
			ptrack = None
	else:
		ptrack = None
	ret_data = dict(ftrack=ftrack, ptrack=ptrack)

	if return_zgoubi_files:
		ret_data['dat_file'] = res.dat()
		ret_data['res_file'] = res.res()

	return ret_data


def plot_cell_properties(data, output_prefix="results/cell_", file_fmt=".pdf", ncells=1):
	"""Produce a set of standard plots from the data structure
	
	data: the data structure returned by get_cell_properties()

	"""
	import os
	from matplotlib import pyplot
	if not hasattr(pyplot, "tight_layout"): pyplot.tight_layout = lambda :None 
	pyplot.clf()

	mkdir_p(os.path.dirname(output_prefix))

	#data = numpy.load(data_file_path)
	n_vars = len(data.dtype.names)

	# drop rows that dont start with zero
	stable_data =  data[data['stable']]

	n_energies =  stable_data.shape[0]
	if n_energies == 0:
		zlog.warn("No stable properties to plot")
		return

	off_x, off_y = 0.1, -3
	# horizontal
	for n in range(n_energies):
		pyplot.plot([stable_data['Y'][n]], [stable_data['T'][n]], 'bx')
		if (n_energies < 5) or ((n%int(n_energies/5))==0):
			pyplot.annotate("%.3f"%(stable_data['KE'][n]/1e6), xy=(stable_data['Y'][n], stable_data['T'][n]),
			               xycoords='data', xytext=(10, 10), textcoords='offset points',
						   arrowprops=dict(arrowstyle="->", color="grey"))
	pyplot.title("Horizontal closed orbit")
	pyplot.xlabel("Y (cm)")
	pyplot.ylabel("T (mrad)")
	if (stable_data['T'].max() - stable_data['T'].min()) < 1e-3:
		y_min = stable_data['T'].min()
		y_max = stable_data['T'].max()
		pyplot.ylim([ y_min-(y_max-y_min)*1 , y_max+(y_max-y_min)*1 ])
	pyplot.tight_layout()
	pyplot.savefig("%s_hor"%(output_prefix)+file_fmt)
	label = "Min = %s, Max = %s" % (stable_data['Y'].min(), stable_data['Y'].max())

	# vert
	pyplot.clf()
	off_x, off_y = 0, 0
	for n in range(n_energies):
		pyplot.plot([stable_data['Z'][n]], [stable_data['P'][n]], 'bx')
		if (n_energies < 5) or ((n%int(n_energies/5))==0):
			pyplot.annotate("%.3f"%(stable_data['KE'][n]/1e6), xy=(stable_data['Z'][n], stable_data['P'][n]),
			               xycoords='data', xytext=(10, 10), textcoords='offset points',
						   arrowprops=dict(arrowstyle="->", color="grey"))
	pyplot.title("Vertical closed orbit")
	pyplot.xlabel("Z (cm)")
	pyplot.ylabel("P (mrad)")
	pyplot.tight_layout()
	pyplot.savefig("%s_ver"%(output_prefix)+file_fmt)
	label = "Min = %s, Max = %s" % (stable_data['Z'].min(), stable_data['Z'].max())

	#tune
	pyplot.clf()
	#pyplot.title("Cell Tune")
	pyplot.xlabel("KE (MeV)")
	pyplot.ylabel("Cell Tune")

	pyplot.plot(stable_data['KE']/1e6,stable_data['NU_Y'], "b-x", label='Horizontal')
	pyplot.plot(stable_data['KE']/1e6,stable_data['NU_Z'], "g-x", label='Vertical')
	pyplot.legend(loc='best')

	pyplot.tight_layout()
	pyplot.savefig("%s_tune"%(output_prefix)+file_fmt)

	if(ncells>1):
		pyplot.clf()
		#pyplot.title("Ring Tune")
		pyplot.xlabel("KE (MeV)")
		pyplot.ylabel("Ring Tune")

		pyplot.plot(stable_data['KE']/1e6,stable_data['NU_Y']*ncells, "b-x", label='Horizontal')
		pyplot.plot(stable_data['KE']/1e6,stable_data['NU_Z']*ncells, "g-x", label='Vertical')
		pyplot.legend(loc='best')

		pyplot.savefig("%s_ring_tune"%(output_prefix)+file_fmt)

	#tof
	pyplot.clf()
	#pyplot.title("Cell time of flight")
	pyplot.xlabel("KE (MeV)")
	pyplot.ylabel("Cell Time of flight ($\mu$s)")
	pyplot.plot(stable_data['KE']/1e6,stable_data['tof'], "bx")
	pyplot.tight_layout()
	pyplot.savefig("%s_tof"%(output_prefix)+file_fmt)

	if(ncells>1):
		pyplot.clf()
		#pyplot.title("Cell time of flight")
		pyplot.xlabel("KE (MeV)")
		pyplot.ylabel("Ring Time of flight ($\mu$s)")
		pyplot.plot(stable_data['KE']/1e6, stable_data['tof']*ncells, "bx")
		pyplot.tight_layout()
		pyplot.savefig("%s_ring_tof"%(output_prefix)+file_fmt)
		label = "Mean time of flight = %s microseconds" % (stable_data['tof']*ncells).mean()

	#freq
	pyplot.clf()
	#pyplot.title("Revolution Frequency")
	pyplot.xlabel("KE (MeV)")
	pyplot.ylabel("Cell revolution frequence (MHz)")
	pyplot.plot(stable_data['KE']/1e6,1/(stable_data['tof']), "bx")
	pyplot.tight_layout()
	pyplot.savefig("%s_freq"%(output_prefix)+file_fmt)
	
	if(ncells>1):
		pyplot.clf()
		#pyplot.title("Revolution Frequency")
		pyplot.xlabel("KE (MeV)")
		pyplot.ylabel("Revolution frequence (MHz)")
		pyplot.plot(stable_data['KE']/1e6,1/(stable_data['tof'])/ncells, "bx")
		pyplot.tight_layout()
		pyplot.savefig("%s_ring_freq"%(output_prefix)+file_fmt)
		label = "Min = %s, Max = %s" % ((1/(stable_data['tof'])/ncells).min(), (1/(stable_data['tof'])/ncells).max())

	#path length
	pyplot.clf()
	#pyplot.title("Path length")
	pyplot.xlabel("KE (MeV)")
	pyplot.ylabel("Cell path length")
	pyplot.plot(stable_data['KE']/1e6,(stable_data['S']), "bx")
	pyplot.tight_layout()
	pyplot.savefig("%s_pathlen"%(output_prefix)+file_fmt)
	
	if(ncells>1):
		pyplot.clf()
		#pyplot.title("Path length")
		pyplot.xlabel("KE (MeV)")
		pyplot.ylabel("Ring path length")
		pyplot.plot(stable_data['KE']/1e6,(stable_data['S'])*ncells, "bx")
		pyplot.tight_layout()
		pyplot.savefig("%s_ring_pathlen"%(output_prefix)+file_fmt)
		label = "Min = %s, Max = %s" % (((stable_data['S'])*ncells).min(), ((stable_data['S'])*ncells).max())
	page_images = []


def cell_properties_table(data, keys, sep="\t"):
	"""Quick function for printing parameter tables, eg:
	
	    print gcp.cell_properties_table(data_arc, ["KE", "stable", "Y", "T", "NU_Y", "NU_Z"])

	"""
	out = ""
	out += sep.join(keys) + "\n"
	for d in data:
		row = [str(d[k]) for k in keys]
		out += sep.join(row) + "\n"
	return out


def get_twiss_params(cell, data, particle, output_prefix="results/twiss_profiles_", full_tracking=False, calc_dispersion=False):
	"""Get the periodic Twiss (Courant and Snyder) parameters for the cell. Tracks a bunch of particles containing pairs offset in each plane.

	Must pass in data returned from get_cell_properties() to give initial conditions. The profile is added to the 'twiss_profile' key in the data structure. If full_tracking is enabled, the Twiss parameters are every integration step are recorded, otherwise parameters are recorded at the end of each magnetic element.


	"""
	if not isinstance(data, GCPData):
		data = GCPData.from_ndarray(data)
	mkdir_p(os.path.dirname(output_prefix))

	tline = Line('test_line')
	tline.add_input_files(cell.input_files)
	ob = OBJET5(DR=1, PY=1e-4,PT=1e-3,PZ=1e-4,PP=1e-3,PX=1e-3,PD=1e-3)
	tline.add(ob)
	part_ob, mass, charge_sign = part_info(particle)
	tline.add(part_ob)
	
	tline.add(DRIFT('start', XL=0))
	tline.add(FAISCNL(FNAME='zgoubi.fai',))

	for e in cell.elements():
		isphysical = False
		try:
			dummy = e.XL
			isphysical = True
		except AttributeError:pass
		try:
			dummy = e.AT
			isphysical = True
		except AttributeError:pass
		tline.add(e)
		if isphysical:
			tline.add(FAISCNL(FNAME='zgoubi.fai',))

	tline.add(DRIFT('end', XL=0))
	tline.add(FAISCNL(FNAME='zgoubi.fai',))
	tline.add(MATRIX(IORD=1,IFOC=11))	
	tline.add(END())

	for n, particle_ke in enumerate(data['KE']):
		if not data[n]['stable']: continue
		print("energy = ", particle_ke)
		rigidity = ke_to_rigidity(particle_ke,mass) / charge_sign
		ob.set(BORO=rigidity)

		ref_Y,ref_T,ref_Z,ref_P = data[n]['Y'], data[n]['T'], data[n]['Z'], data[n]['P']
		ob.set(YR=ref_Y, TR=ref_T, ZR=0, PR=0, XR=0)
		init_twiss = twiss_param_array()
		pn = ""
		if not data.info["periodic"]:
			pn="0"
		init_twiss["beta_y"] = data[n]["BETA_Y"+pn]
		init_twiss['alpha_y'] = data[n]['ALPHA_Y'+pn]
		init_twiss['gamma_y'] = data[n]['GAMMA_Y'+pn]
		init_twiss['disp_y'] = data[n]['DISP_Y'+pn]
		init_twiss['disp_py'] = data[n]['DISP_PY'+pn]
		init_twiss['beta_z'] = data[n]['BETA_Z'+pn]
		init_twiss['alpha_z'] = data[n]['ALPHA_Z'+pn]
		init_twiss['gamma_z'] = data[n]['GAMMA_Z'+pn]
		init_twiss['disp_z'] = data[n]['DISP_Z'+pn]
		init_twiss['disp_pz'] = data[n]['DISP_PZ'+pn]
		#data[n][["BETA_Y", "ALPHA_Y", "GAMMA_Y", "DISP_Y", "DISP_PY", "BETA_Z", "ALPHA_Z", "GAMMA_Z", "DISP_Z", "DISP_PZ"]]
	
		tline.full_tracking(False)
		twiss_profiles = get_twiss_profiles(tline,'%s%s.txt'%(output_prefix, particle_ke), input_twiss_parameters=init_twiss, calc_dispersion=calc_dispersion)
		data[n]['twiss_profile'] = twiss_profiles
		if full_tracking:
			tline.full_tracking(True)
			twiss_profiles = get_twiss_profiles(tline, '%s%s_full.txt'%(output_prefix, particle_ke), calc_dispersion=calc_dispersion, input_twiss_parameters=init_twiss)
			data[n]['full_twiss_profile'] = twiss_profiles


def plot_twiss_params(data, output_prefix="results/twiss_profiles_", fields=None):
	"""Plot the Twiss profiles found with get_twiss_params()

	fields: list of fields to plot, beta_y, beta_z, alpha_y, alpha_z, gamma_y, gamma_z, disp_y, disp_z, disp_py, disp_pz
	"""
	from matplotlib import pyplot
	if fields is None:
		fields = ["beta_y", "beta_z", "alpha_y", "alpha_z"]

	stable_data =  data[data['stable']].copy()

	field_cats = []
	for f in fields:
		if f not in stable_data[0]['twiss_profile'].dtype.names:
			raise ValueError("Not a valid twiss key: %s"%f)
		cat = f[:-1]
		if cat not in field_cats: field_cats.append(cat)
	
	if len(field_cats) > 2:
		raise ValueError("More than 2 types of field not currently supported")
		

	full_tracking_message = 0
	for n, particle_ke in enumerate(stable_data['KE']):
		if stable_data[n]['full_twiss_profile'] is not None:
			print(stable_data[n]['full_twiss_profile'])
			twiss_profiles = stable_data[n]['full_twiss_profile']
		elif stable_data[n]['twiss_profile'] is not None:
			if full_tracking_message == 0:
				zlog.warn("get_twiss_params() called without full_tracking=True, so only plotting twiss at element ends")
				full_tracking_message = 1
			twiss_profiles = stable_data[n]['twiss_profile']
		else:
			zlog.error("Call get_twiss_params() before plot_twiss_params()")
			raise ValueError


		pyplot.clf()
		lns = []
		ax1 = pyplot.axes()
		for f in fields:
			if f.startswith(field_cats[0]):
				l = ax1.plot(twiss_profiles['s'],twiss_profiles[f],"-", label=f)
				lns.append(l)
		ax1.set_xlabel("Path length (m)")
		ax1.set_ylabel(field_cats[0].partition("_")[0], color='k')

		if len(field_cats) > 1:
			ax2 = ax1.twinx()
			for f in fields:
				if f.startswith(field_cats[1]):
					l = ax2.plot(twiss_profiles['s'],twiss_profiles[f],"--", label=f)
					lns.append(l)
			ax2.set_ylabel(field_cats[1].partition("_")[0], color='k')

		lns = sum(lns, [])
		pyplot.legend(lns, [l.get_label() for l in lns],loc="best") # stackoverflow.com/questions/5484922
		pyplot.savefig('%s%s.pdf'%(output_prefix, particle_ke))


def plot_cell_tracks(cell, data, particle, output_file="results/track.pdf", show=False, plot_unstable=False, draw_field_midplane=False, sector_width=None, aspect="equal", style=None, draw_tracks=True, min_y=None, max_y=None, y_steps=None, angle=0, plot_extents=None):
	"""Plot particle track through cell, starting from in the closed orbits stored in data.

	output_file: file to save plot
	show: display plot to screen
	plot_unstable: show tracks for initial coordinates marked as unstable in data structure
	draw_field_midplane: show the magnet field on the magnet midplane
	sector_width: width to draw sector magnets
	aspect: "equal", use same scale for x and y, "auto" stretch to fit
	style: style for plotting, see lab_plot.LabPlot
	draw_tracks: draw particle tracks
	min_y, max_y, y_steps, angle: starting coordinates for test particles for sampling the midplane field
	plot_extents: list of extents [left, right, bottom, top]

	"""
	from zgoubi.lab_plot import LabPlot
	import pylab
	pylab.close()
	n_vars = len(data.dtype.names)

	if plot_unstable==False:
		# drop rows that dont start with zero
		stable_data =  data[data['stable']].copy()
	else:
		stable_data = data.copy()
	
	cell = copy.deepcopy(cell)
	cell.full_tracking(True, drift_to_multi=True)
	#cell = uniquify_labels(cell)
	
	cell2 = Line("test_line")
	#cell2.add(DRIFT("gcpstar", XL=0* cm_)) # not sure if needed, but they probably show up a bug in labplot when using an FFAG
	cell2.add(cell)
	#cell2.add(DRIFT("gcpend", XL=0* cm_))

	part_ob, mass, charge_sign = part_info(particle)

	# if line contains BENDS, then line will be drawn with first energy, but zgoubi will adjust angles for other particles
	if len(data['KE']>0):
		boro = ke_to_rigidity(data['KE'][0],mass) / charge_sign
	else:
		boro = None
	
	lp = LabPlot(cell2, boro=boro, sector_width=sector_width, aspect=aspect)
	lp.set_noel_offset(3) # get_cell_tracks adds 3 elements to start of line
	if style:
		lp.set_style(style)
	

	if draw_tracks:
		get_cell_tracks(cell2, stable_data, particle, full_tracking=True, add_faiscnl=False)
		for d in stable_data:
			ftrack = d["ftrack"]
			ptrack = d["ptrack"]
			if ftrack is None and ptrack is None:
				zlog.error("Failed to read fai or plt files")
			lp.add_tracks(ftrack=ftrack, ptrack=ptrack, draw=1)
	
	if draw_field_midplane:
		if min_y is None or max_y is None or y_steps is None:
			raise ValueError("When using draw_field_midplane, you must set min_y, max_y and y_steps")

		for Y in numpy.linspace(min_y, max_y, y_steps):
			track = get_tracks(cell2, [Y, angle, 0, 0], particle, 1e15, full_tracking=True, add_faiscnl=False)
			lp.add_tracks(ftrack=track["ftrack"], ptrack=track["ptrack"], draw=0, field=1)


	if draw_field_midplane:
		lp.draw(draw_field_midplane=draw_field_midplane, draw_tracks=draw_tracks, field_steps=y_steps, plot_extents=plot_extents)
	else:
		lp.draw(plot_extents=plot_extents)

	lp.save(output_file)
	if show:
		lp.show()


def profile_get_tracks(magnet, min_y, max_y, y_steps, angle=0):
	"""Internal function used by profile1d() and profile2d()

	"""

	tline = Line('t')
	try:
		tline.add_input_files(magnet.input_files)
	except AttributeError:
		pass
	tline.add(magnet)

	# test bunch of high rigidity particles
	testparticles = GCPData(y_steps, info=dict(periodic=True, particle="p"))
	testparticles['Y'] = numpy.linspace(min_y, max_y, y_steps)
	testparticles['T'] = angle
	testparticles['KE'] = 1e20
	testparticles['stable'] =True

	get_cell_tracks(tline, testparticles, 'p', full_tracking=True)

	fieldpoints_list = []
	for pt in testparticles['ptrack']:
		fieldpoints = []
		if pt is None:
			zlog.warn('Empty track')
			continue
		for pp in pt:
			y, z, x = pp['Y'], pp['Z'], pp['X']
			by, bz, bx = pp['BY'], pp['BZ'], pp['BX']
			if isnan(x) or isnan(y) or isnan(z):
				continue
			fieldpoints.append([y,z,x,by, bz, bx])
	
		field_map_data = numpy.array(fieldpoints)
		fieldpoints_list.append(field_map_data)
	
	return fieldpoints_list


def profile2d(magnet, min_y, max_y, y_steps, angle=0):
	"""Get a 2d interpolated field profile. Used by plot_element_fields

	returns a y_steps by y_steps grid of interpolated field values, and the extents of that grid:

	    returns int_field, (xmin,xmax,ymin,ymax)

	"""
	import scipy.interpolate
	if not hasattr(scipy.interpolate, "griddata"):
		print("profile2d() requires scipy > 0.9")
		raise

	field_map_data = profile_get_tracks(magnet, min_y, max_y, y_steps, angle)
	
	# due to scales being different we need to rescale
	# see https://github.com/scipy/scipy/issues/2975
	# get typecal steps, in x and y. should be close to xpas, and gap between orbits
	x_diff_means = []
	y_pos_means = []
	for fmp in field_map_data:
		fmp_xs = fmp[:,2]
		fmp_ys = fmp[:,0]
		x_difs = fmp_xs[1:] - fmp_xs[:-1]
		x_diff_means.append(x_difs.mean())
		y_pos_means.append(fmp_ys.mean())

	x_step = numpy.array(x_diff_means).mean()
	y_pos_means = numpy.array(y_pos_means)
	y_step =  (y_pos_means[1:] - y_pos_means[:-1]).mean()

	field_map_data = numpy.vstack(field_map_data)
	if len(field_map_data) == 0:
		zlog.error("No tracks, can't make field profile")
		raise NoTrackError

	points = field_map_data[:,0:3:2]
	values = field_map_data[:,4] #bz

	ymin,ymax,xmin,xmax = points[:,0].min(), points[:,0].max(), points[:,1].min(), points[:,1].max()

	# scale positions, we return orginal extents, so this does not effect results
	points[:,0] /= y_step
	points[:,1] /= x_step
	mean_sy = points[:,0].mean()
	points[:,0] -= mean_sy
	symin,symax,sxmin,sxmax = points[:,0].min(), points[:,0].max(), points[:,1].min(), points[:,1].max()

	nsteps = 1j*y_steps*2 # j because of how mgrid works
	grid_y, grid_x = numpy.mgrid[symin:symax:nsteps, sxmin:sxmax:nsteps]

	if numpy.abs(field_map_data[:,1]).max() > 1e-6:
		zlog.warn("Some field points are not at z=0. Plot will be of projection onto z=0")

	int_field = scipy.interpolate.griddata(points, values, (grid_y, grid_x), method="linear")

	return int_field, (xmin,xmax,ymin,ymax)


def profile1d(magnet, min_y, max_y, y_steps, angle=0):
	"""Returns a 1d transverse horizontal (radial) profile of the magnet. At mid point, or ACN for DIPOLES. Used by plot_element_fields()
	
	returns field, ys
	"""
	import scipy.interpolate
	if not hasattr(scipy.interpolate, "griddata"):
		print("profile1d() requires scipy > 0.9")
		raise

	field_map_data = profile_get_tracks(magnet, min_y, max_y, y_steps, angle)
	
	# due to scales being different we need to rescale
	# see https://github.com/scipy/scipy/issues/2975
	# get typecal steps, in x and y. should be close to xpas, and gap between orbits
	x_diff_means = []
	y_pos_means = []
	for fmp in field_map_data:
		fmp_xs = fmp[:,2]
		fmp_ys = fmp[:,0]
		x_difs = fmp_xs[1:] - fmp_xs[:-1]
		x_diff_means.append(x_difs.mean())
		y_pos_means.append(fmp_ys.mean())

	x_step = numpy.array(x_diff_means).mean()
	y_pos_means = numpy.array(y_pos_means)
	y_step =  (y_pos_means[1:] - y_pos_means[:-1]).mean()

	field_map_data = numpy.vstack(field_map_data)
	if len(field_map_data) == 0:
		zlog.error("No tracks, can't make field profile")
		raise NoTrackError

	points = field_map_data[:,0:3:2]
	values = field_map_data[:,4] #bz

	ymin,ymax,xmin,xmax = points[:,0].min(), points[:,0].max(), points[:,1].min(), points[:,1].max()

	xmid = (xmax + xmin)/2
	if hasattr(magnet, "elements"):
		magnet = magnet.elements().next()
	if magnet._zgoubi_name in ["DIPOLES"]:
		xmid = radians(magnet._looped_data[0]['ACN'])

	# scale positions, we return original extents, so this does not effect results
	points[:,0] /= y_step
	points[:,1] /= x_step
	points[:,0] -= points[:,0].mean()
	sxmid = xmid / x_step
	symin,symax,sxmin,sxmax = points[:,0].min(), points[:,0].max(), points[:,1].min(), points[:,1].max()

	nsteps = 1j * y_steps # j because of how mgrid works
	grid_y, grid_x = numpy.mgrid[symin:symax:nsteps, sxmid:sxmid:1j]
	if numpy.abs(field_map_data[:,1]).max() > 1e-6:
		zlog.warn("Some field points are not at z=0. Plot will be of projection onto z=0")

	int_field = scipy.interpolate.griddata(points, values, (grid_y, grid_x), method="linear")
	int_field = int_field.reshape([-1])

	return int_field, numpy.linspace(ymin,ymax,y_steps)


def plot_element_fields(cell, min_y, max_y, y_steps, angle=None, output_prefix_rad="results/mag_field_rad_%s.pdf", output_prefix_long="", output_prefix_2d="", rad_method="interp", extra_test_range=0.5):
	"""Plots radial, longitudinal and 2d mid-plane magnetic field for each element, using high rigidity rays to sample field.
	
	min_y, max_y, y_steps, angle set the starting positions and angle for the test rays
	If angle=None, then 0 will be used for rectangular magnets, and AT/2 for sector magnets

	output_prefix_rad, output_prefix_long, output_prefix_long should contain a '%s' that is replaced with the element label, eg
	"results/mag_field_rad_%s.pdf"
	output_prefix_long and output_prefix_2d will be skipped if they are left blank
	The longitudinal profile is along the straight lines of the test rays, not the closed orbit.

	For the radial plot method can be:
	interp: points along a radial line are interpolated from a 2d profile
	max: use the max field seen by the ray

	extra_test_range factor to extend the range (min_y, max_y) for test particles, to improve interpolation edges
	
	In 2d plots the trick of making rectangular magnets by using large radii is recognised for radii over 1e6cm

	"""
	from matplotlib import pyplot
	if not hasattr(pyplot, "tight_layout"): pyplot.tight_layout = lambda :None

	estart = min_y - (max_y-min_y)*extra_test_range
	estop = max_y + (max_y-min_y)*extra_test_range
	ecount = int(y_steps *(1+extra_test_range))

	for e in cell.elements():
		ismagnet = False
		if angle:
			this_angle = angle
		else:
			this_angle = 0
		try:
			dummy = e.XL
			ismagnet = True
			mag_shape = "rect"
		except AttributeError:pass
		try:
			dummy = e.AT
			ismagnet = True
			mag_shape = "sect"
			if angle is None:
				this_angle = radians(-e.AT/2)*1000 # AT in deg, Y in mrad
		except AttributeError:pass
		if e._zgoubi_name == "DRIFT":
			ismagnet = False

		if e._zgoubi_name == "POLARMES":
			ismagnet = True
			mag_shape = "sect"


		if ismagnet:
			mag_rm = 0
			if hasattr(e, "AT") and hasattr(e, "RM"):
				mag_rm = e.RM
			tline = Line("m")
			tline.add_input_files(cell.input_files)
			tline.add(e)

			if rad_method == "interp":
				bz, ys = profile1d(tline, estart, estop, ecount, this_angle)
				ys -= mag_rm
				pyplot.clf()
				mask = numpy.logical_and(ys>=min_y, ys<=max_y)
				pyplot.plot(ys[mask], bz[mask], '-b')
			elif rad_method == "max":
				testparticles = GCPData(ecount, info=dict(periodic=True, particle="p"))
				testparticles['Y'] = numpy.linspace(min_y, max_y, y_steps)
				testparticles['T'] = this_angle
				testparticles['KE'] = 1e15
				testparticles['stable'] =True
				get_cell_tracks(tline, testparticles, 'p', full_tracking=True)

				max_bz = [ pt["BZ"].max() for pt in testparticles['ptrack'] ]
				pyplot.plot(testparticles['Y'], max_bz, "-b")
			else:
				raise ValueError('rad_method should be "interp" or "max"')

			pyplot.xlabel("Y (cm)")
			pyplot.ylabel("$B_z$ (kGauss)")
			pyplot.grid()
			pyplot.tight_layout()
			pyplot.savefig(output_prefix_rad%e.label1)
			pyplot.clf()

			if output_prefix_long:
				testparticles = GCPData(y_steps, info=dict(periodic=True, particle="p"))
				testparticles['Y'] = numpy.linspace(min_y, max_y, y_steps)
				testparticles['T'] = this_angle
				testparticles['KE'] = 1e15
				testparticles['stable'] =True

				get_cell_tracks(tline, testparticles, 'p', full_tracking=True)
				pyplot.clf()
				for pt in testparticles['ptrack']:
					pyplot.plot(pt["S"], pt["BZ"], '-b')
				pyplot.xlabel("S (cm)")
				pyplot.ylabel("$B_z$ (kGauss)")
				pyplot.grid()

				pyplot.tight_layout()
				pyplot.savefig(output_prefix_long%e.label1)
				pyplot.clf()

			if output_prefix_2d:
				bz2, extents = profile2d(tline, estart, estop, ecount, this_angle)
				#check for rectangular magnets made with DIPOLES.
				if e.RM > 1e6:
					# looks like using DIPOLES with large RM to fake rectangular magnet, so convert radians to cm
					extents = extents[0]*mag_rm, extents[1]*mag_rm, extents[2], extents[3]

				extents = [extents[0], extents[1], extents[2]-mag_rm, extents[3]-mag_rm]

				gridn = bz2.shape[0]
				grid_base = extents[2]
				grid_size = extents[3] - extents[2]
				nstart = int(floor((min_y - extents[2]) / grid_size * gridn))
				nstop = int(ceil((max_y - extents[2]) / grid_size * gridn))
				bz2 = bz2[nstart:nstop+1,:]
				extents[2] = nstart * grid_size/gridn + grid_base
				extents[3] = nstop * grid_size/gridn + grid_base

				bz_max = numpy.nanmax(numpy.abs(bz2))
				pyplot.clf()
				pyplot.imshow(bz2, extent=extents, origin="lower", vmin=-bz_max, vmax=bz_max, aspect="auto", interpolation="bilinear", cmap="bwr")

				cbar = pyplot.colorbar()
				cbar.ax.set_ylabel("$B_z$ (kGauss)")
				if mag_shape == "rect" or mag_rm > 1e6:
					pyplot.xlabel("S (cm)")
				else:
					pyplot.xlabel("S (radians)")
				pyplot.ylabel("Y (cm)")
				pyplot.tight_layout()
				pyplot.savefig(output_prefix_2d%e.label1)
				pyplot.clf()


def get_dynamic_aperture(cell, data, particle, npass, nangles=3, tol=0.01, quick_mode=False, debug_log=None, island_avoid=0.01, start = 1e-6):
	"""Get Dynamic Aperture.
	
	cell: the cell to run
	data: data structure containing the closed orbit to be used initial conditions, see get_cell_properties()

	tol: tolerance 0.01 give result to 1% precision
	npass: number of passes through the cell
	quick_mode: "+y+z" 1 particle at [dY,0,dZ,0]
	            "+t+p" 1 particle at [0,dT,0,dP]
	            False:16 combinations of plus and minus offsets are used
	nangles: number of realspace angles to use. 1 will give the 45degree DA, i.e. with equal horizontal and vertical excitation, 2 will give pure horizontal (0 deg) and pure vertical (90deg), for other numbers will use numpy.linspace(0,pi/2,nangles).
	island_avoid: when an unstable amplitude is found take a small step up, to see if it just a small island. set to zero to disable
	debug_log: file name to write debug information to
	start: starting emittance 

	From the starting emittance a search for the stability boundary is made.

	"""

	if debug_log:
		debug_log = open(debug_log, "w")
		print("get_dynamic_aperture, npass=", npass, " nangles=", nangles, " quick_mode=", quick_mode, file=debug_log)
		if data['stable'].sum() == 0:
			print("no stable orbits", file=debug_log)

	min_val = 1e-10 # unstable if below
	step = 5 # step factor when unbounded

	if nangles == 1:
		angles = numpy.array([pi/4])
	else:
		angles = numpy.linspace(0,pi/2,nangles)

	def get_cur(cur, bound_min, bound_max):
		if bound_min is None and bound_max is None:
			# first lap
			cur = start
		elif bound_min is None:
			# not found bottom bound yet
			cur /= step
		elif bound_max is None:
			# not found top bound yet
			cur *= step
		else:
			# bisect
			cur = (bound_min + bound_max) /2
			if cur == bound_min or cur == bound_max:
				print("WARN: Bounds equal, machine precision reached. Reduce tolerance")
		return cur

	def update_bounds(cur, stab, bound_min, bound_max):
		if stab and (bound_min is None or bound_min < cur):
			bound_min = cur
		elif not stab and (bound_max is None or bound_max > cur):
			bound_max = cur
		else:
			print("WARN:  bounds unchanged", cur, stab, bound_min, bound_max, " Reduce tolerance")

		return bound_min, bound_max
	
	def is_stable(cur_emit, tline, angle, it):
		emit_h = sin(pi/2 - angle) * cur_emit # like 'cos(angle)' but goes to zero better
		emit_v = sin(angle) * cur_emit


		# get offsets from closed orbit
		current_YTZP = emittance_to_coords(emit_h, emit_v, [alpha_y,alpha_z], [beta_y, beta_z])
		dY, dT, dZ, dP = current_YTZP[0][0], current_YTZP[1][1], current_YTZP[0][2], current_YTZP[1][3]
		if debug_log:
			print("cur_emit=", cur_emit, "  emit_h=", emit_h, " emit_v=", emit_v, file=debug_log)
			print("dY, dT, dZ, dP", dY, dT, dZ, dP, file=debug_log)

		#print dY, dT, dZ, dP
		ob = tline.get_objet()
		ob.clear()
		pn = 0
		if quick_mode:
			if quick_mode == "+y+z":
				ob.add(Y=Yc+dY, T=Tc, Z=Zc+dZ, P=Pc, LET="A", D=1)
			elif quick_mode == "+t+p":
				ob.add(Y=Yc, T=Tc+dT, Z=Zc, P=Pc+dP, LET="B", D=1)
			else:
				raise ValueError('quick mode must be "+y+z" or "+t+p"')
			pn=1
		else:
			for yt_coords in [[dY,0],[-dY,0],[0,dT],[0,-dT]]:
				for zp_coords in [[dZ,0],[-dZ,0],[0,dP],[0,-dP]]:
					#print pn,  yt_coords + zp_coords
					Ye1, Te1, Ze1, Pe1 = yt_coords + zp_coords
					ob.add(Y=Yc+Ye1, T=Tc+Te1, Z=Zc+Ze1, P=Pc+Pe1, LET=chr(ord('A')+(pn%26)), D=1)
					pn+=1
		
		print("DA: angle %.2f iteration %s"%(degrees(angle), it))
		res = tline.run(xterm =0)
		
		if res.test_rebelote(): # atlease one particle survived
			iexs = res.get_all(file="fai")['IEX'] # get losses
			stab_c = sum(iexs==1)
			if stab_c == pn:
				stab = True
			else:
				stab = False

		else: # all particles lost
			stab = False
			stab_c = 0
		if debug_log:
			print("stab=", stab, "  stab_c=", stab_c, file=debug_log)

		res.clean()

		return stab, pn, stab_c

	def is_stable_test(cur_emit, tline, angle, it):
		da     = 0.12345
		island = 0.00312
		#island_size = island * 0.02 # big island
		island_size = island * 0.01 # big island
		#island_size = island * 0.001 # tiny island

		in_island = (  island-island_size < cur_emit < island+island_size  )
		stab = (cur_emit < da and not in_island)
		if debug_log:
			print("In island", in_island, file=debug_log)
		return stab, 16, 0

	tline = Line('test_line')
	tline.add_input_files(cell.input_files)
	ob = OBJET2()
	tline.add(ob)
	part_ob, mass, charge_sign = part_info(particle)
	tline.add(part_ob)
	tline.add(cell)
	

	tline.add(REBELOTE(NPASS=npass, K=99))
	tline.add(DRIFT("end",XL=1e-12))
	tline.add(FAISCNL("end", FNAME='zgoubi.fai'))
	tline.add(END())
	tline.full_tracking(False)

	for n, particle_ke in enumerate(data['KE']):
		if not data[n]['stable']: continue
		print("energy = ", particle_ke)
		rigidity = ke_to_rigidity(particle_ke,mass) / charge_sign
		ob.set(BORO=rigidity)
		# closed orbit
		Yc, Tc, Zc, Pc = [data[n]['Y'], data[n]['T'],data[n]['Z'],data[n]['P'] ]
		alpha_y, beta_y, alpha_z, beta_z = data[n]['ALPHA_Y'], data[n]['BETA_Y'],data[n]['ALPHA_Z'],data[n]['BETA_Z']
		
		data[n]['DA'] = numpy.zeros([nangles])
		data[n]['DA_angles'] = angles

		for an, angle in enumerate(angles):

			cur_emit = 0 # just need to create it, get set by get_cur
			it=0
			bound_min = None
			bound_max = None
			if debug_log:
				print("Start search KE=", particle_ke, " angle=", angle, file=debug_log)
			
			while it<100:
				it += 1
				cur_emit = get_cur(cur_emit, bound_min, bound_max)
				if debug_log:
					print("\ncur_emit=", cur_emit, file=debug_log)

				stab, count_p, count_s  = is_stable(cur_emit, tline, angle, it)

				if not stab and island_avoid != 0 and not (count_p == 16 and count_s <= 8):
					# check if this is just a small unstable island
					stab, count_p, count_s = is_stable(cur_emit*(1+island_avoid), tline, angle, it)
					if stab:
						print("Stepped over small unstable island at", cur_emit, file=debug_log)

				if debug_log:
					print("stab=", stab, file=debug_log)
					debug_log.flush()

				bound_min, bound_max = update_bounds(cur_emit, stab, bound_min, bound_max)

				if bound_min is not None and bound_max is not None and ((bound_max-bound_min)/bound_min < tol):
					data[n]['DA'][an] = bound_min
					break

				if not stab and cur_emit < min_val:
					zlog.warn("Not stable above min_val")
					data[n]['DA'][an] = 0
					break
			else:
				zlog.warn("Maximum iterations reached")
				data[n]['DA'][an] = 0

		print("DA", data[n]['DA'])
		

def get_phase_space(cell, data, particle, npass, emits=None):
	"""Track particle for npass turns

	"""

	tline = Line('test_line')
	tline.add_input_files(cell.input_files)
	ob = OBJET2()
	tline.add(ob)
	part_ob, mass, charge_sign = part_info(particle)
	tline.add(part_ob)
	tline.add(cell)
	tline.add(FAISCNL("end", FNAME='zgoubi.fai'))
	tline.add(REBELOTE(NPASS=npass, K=99))
	tline.add(END())
	tline.full_tracking(False)

	for n, particle_ke in enumerate(data['KE']):
	#	if not data[n]['stable']: continue
		if not data[n]['found_co']: continue
		print("get_phase_space ke", particle_ke, "n")
		rigidity = ke_to_rigidity(particle_ke,mass) / charge_sign
		ob = tline.get_objet()
		ob.set(BORO=rigidity)
		Yc, Tc, Zc, Pc = [data[n]['Y'], data[n]['T'],data[n]['Z'],data[n]['P'] ]
		alpha_y, beta_y, alpha_z, beta_z = data[n]['ALPHA_Y'], data[n]['BETA_Y'],data[n]['ALPHA_Z'],data[n]['BETA_Z']
		
		ob.clear()
		for m, emit in enumerate(emits):
			# horizontal particle
			current_YTZP = emittance_to_coords(emit/sqrt(2), emit/sqrt(2), [alpha_y,alpha_z], [beta_y, beta_z])
			Ye1, Te1, Ze1, Pe1 = current_YTZP[0][0], current_YTZP[0][1], current_YTZP[0][2], current_YTZP[0][3]
			ob.add(Y=Yc+Ye1, T=Tc+Te1, Z=Zc+0, P=Pc+0, LET='A', D=1)

			# vertical particle
			current_YTZP = emittance_to_coords(emit*emit/sqrt(2), emit/sqrt(2), [alpha_y,alpha_z], [beta_y, beta_z])
			Ye1, Te1, Ze1, Pe1 = current_YTZP[0][0], current_YTZP[0][1], current_YTZP[0][2], current_YTZP[0][3]
			ob.add(Y=Yc+0, T=Tc+0, Z=Zc+Ze1, P=Pc+Pe1, LET='B', D=1)
			
		res = tline.run()
		if res.test_rebelote():
			stab = True
			fai_data = res.get_all('fai')
		else:
			stab = False

			try:
				fai_data = res.get_all('fai')
			except IOError:
				print("No fai file")
				continue
		data[n]['phase_space'] = fai_data



def plot_dynamic_aperture(cell, data, particle, npass, output_prefix="results/da_"):
	"""Make phase space plots showing dynamic aperture

	"""

	from matplotlib import pyplot
	mkdir_p(os.path.dirname(output_prefix))

	tline = Line('test_line')
	tline.add_input_files(cell.input_files)
	ob = OBJET2()
	tline.add(ob)
	part_ob, mass, charge_sign = part_info(particle)
	tline.add(part_ob)
	tline.add(cell)
	tline.add(FAISCNL("end", FNAME='zgoubi.fai'))
	tline.add(REBELOTE(NPASS=npass, K=99))
	tline.add(END())
	tline.full_tracking(False)

	for n, particle_ke in enumerate(data['KE']):
		if not data[n]['stable']: continue
		pyplot.clf()
		print("energy = ", particle_ke)
		rigidity = ke_to_rigidity(particle_ke,mass) / charge_sign
		ob = tline.get_objet()
		ob.set(BORO=rigidity)
		Yc, Tc, Zc, Pc = [data[n]['Y'], data[n]['T'],data[n]['Z'],data[n]['P'] ]
		alpha_y, beta_y, alpha_z, beta_z = data[n]['ALPHA_Y'], data[n]['BETA_Y'],data[n]['ALPHA_Z'],data[n]['BETA_Z']
		

		emit_h_max_stable = data[n]['DA'][0]
		emit_v_max_stable = data[n]['DA'][-1]
		
		tracks = []
		# then test a range from small to unstable
		#for m, emit in enumerate(numpy.linspace(1e-12, 2, 20)): # linearly spaced in emittance
		for m, emit in enumerate((numpy.linspace(1e-12, (2)**0.5, 20))**2): # linearly spaced on plot
			print("plot_da ke", particle_ke, "n", n, "m", m)
			ob.clear()
			# horizontal particle
			current_YTZP = emittance_to_coords(emit*emit_h_max_stable/sqrt(2), emit*emit_h_max_stable/sqrt(2), [alpha_y,alpha_z], [beta_y, beta_z])
			Ye1, Te1, Ze1, Pe1 = current_YTZP[0][0], current_YTZP[0][1], current_YTZP[0][2], current_YTZP[0][3]
			ob.add(Y=Yc+Ye1, T=Tc+Te1, Z=Zc+0, P=Pc+0, LET='A', D=1)

			# vertical particle
			current_YTZP = emittance_to_coords(emit*emit_v_max_stable/sqrt(2), emit*emit_v_max_stable/sqrt(2), [alpha_y,alpha_z], [beta_y, beta_z])
			Ye1, Te1, Ze1, Pe1 = current_YTZP[0][0], current_YTZP[0][1], current_YTZP[0][2], current_YTZP[0][3]
			ob.add(Y=Yc+0, T=Tc+0, Z=Zc+Ze1, P=Pc+Pe1, LET='B', D=1)
			
			res = tline.run()
			if res.test_rebelote():
				stab = True
				fai_data = res.get_all('fai')
			else:
				stab = False

				try:
					fai_data = res.get_all('fai')
				except IOError:
					print("No fai file")
					continue
			tracks.append([fai_data,emit,stab])

			#pyplot.plot(fai_data['Y'], fai_data['T'], ','+color)

		import pickle
		pickle.dump(tracks, open('%s_%s.pickle'%(output_prefix, particle_ke), "w"))

		print("dynamic aperture", numpy.median(data[n]['DA']))
		for fai_data,emit,stab in tracks:
			fai_data = fai_data[fai_data["LET"]=="A"]
			print(emit, stab, len(fai_data))
			stab = fai_data["PASS"].max() == npass+1
			if stab:
				pyplot.plot(fai_data['Y'], fai_data['T'], ',k')
				plot_width_Y = fai_data['Y'].max() - fai_data['Y'].min()
				plot_width_T = fai_data['T'].max() - fai_data['T'].min()
				#pyplot.annotate("%4g"%(emit*1e6), [fai_data['Y'][0], fai_data['T'][0]])
			else:
				pyplot.plot(fai_data['Y'], fai_data['T'], ',r')

		pyplot.xlabel("Y (cm)")
		pyplot.ylabel("T (mrad)")
		pyplot.xlim(Yc-0.6*plot_width_Y, Yc+0.6*plot_width_Y)
		pyplot.ylim(Tc-0.6*plot_width_T, Tc+0.6*plot_width_T)
		pyplot.savefig('%sh_%s.pdf'%(output_prefix, particle_ke))

		for fai_data,emit,stab in tracks:
			fai_data = fai_data[fai_data["LET"]=="B"]
			stab = fai_data["PASS"].max() == npass+1
			if stab:
				pyplot.plot(fai_data['Z'], fai_data['P'], ',k')
				plot_width_Z = fai_data['Z'].max() - fai_data['Z'].min()
				plot_width_P = fai_data['P'].max() - fai_data['P'].min()
				#pyplot.annotate("%4g"%(emit*1e6), [fai_data['Z'][0], fai_data['P'][0]])
			else:
				pyplot.plot(fai_data['Z'], fai_data['P'], ',r')
		pyplot.xlabel("Z (cm)")
		pyplot.ylabel("P (mrad)")
		pyplot.xlim(Zc-0.6*plot_width_Z, Zc+0.6*plot_width_Z)
		pyplot.ylim(Pc-0.6*plot_width_P, Pc+0.6*plot_width_P)
		pyplot.savefig('%sv_%s.pdf'%(output_prefix, particle_ke))
