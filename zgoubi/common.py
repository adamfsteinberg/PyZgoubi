#!/usr/bin/env python
"""Common functions that can be used anywhere in pyzgoubi

"""

from __future__ import division, print_function
import errno
import os

def mkdir_p(path):
	"""Make a directory. Make parents if needed. Don't raise error if it already exists

	From: http://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python code by tzotzioy http://stackoverflow.com/users/6899/
	"""
	try:
		os.makedirs(path)
	except OSError as exc:
		if exc.errno == errno.EEXIST:
			pass
		else:
			raise

def open_file_or_name(forn, mode="r", mkdir=False):
	"""Pass either a filename or file handle like object. Returns a file like object

	If mode is 'w' or 'a' then can set mkdir=True to make sure the directory is created if needed.
	"""
	if hasattr(forn, 'readline'):
		return forn
	else:
		if mkdir and mode in ["w", "a"] and os.path.dirname(forn):
			mkdir_p(os.path.dirname(forn))
		return open(forn, mode)


def twiss_param_array(**kwargs):
	"""Helper function to make the twiss param array used as input to get_twiss_profiles() or get_cell_properties_nonperiodic()

	tp = twiss_param_array(beta_y=2, alpha_y=-1 ...)
	all unspecified will be left at 0

	"""
	import numpy
	tp = numpy.zeros(1, dtype=[('beta_y', 'f8'), ('alpha_y', 'f8'), ('gamma_y', 'f8'), ('disp_y', 'f8'), ('disp_py', 'f8'),
	                        ('beta_z', 'f8'), ('alpha_z', 'f8'), ('gamma_z', 'f8'), ('disp_z', 'f8'), ('disp_pz', 'f8')])
	for k,v in kwargs.items():
		tp[k.lower()] = v
	
	if "beta_y" in kwargs and 'alpha_y' in kwargs and not 'gamma_y' in kwargs:
		tp['gamma_y'] = (1+tp['alpha_y']**2) / tp['beta_y']
	if "beta_z" in kwargs and 'alpha_z' in kwargs and not 'gamma_z' in kwargs:
		tp['gamma_z'] = (1+tp['alpha_z']**2) / tp['beta_z']

	return tp

