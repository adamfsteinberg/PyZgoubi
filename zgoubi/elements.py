#!/usr/bin/env python
# -*- coding: utf-8 -*-

#       pyzgoubi - python interface to zgoubi
#       Copyright 2008-2015 Sam Tygier <Sam.Tygier@hep.manchester.ac.uk>
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.

# a base class for all the beam line objects

from __future__ import division, print_function
import zgoubi_metadata.elements
from numbers import Number

__all__ = ["zgoubi_element", "tXPAS", "zgoubi_element_def"]

class zgoubi_element(object):
	"A base class for zgoubi elements"
	def __init__(self):
		pass

	def set_param(self, key, val):
		if key in self._params.keys():
			self._params[key] = val
		else:
			raise ValueError("no such param: '" + str(key) + "' In element " + self._zgoubi_name)

	def set(self, *dsettings, **settings):
		"""Set a parameter value::
			my_element.set(XL=5)

		can also use a dictionary to set values::
			s = {'XL':5, B_0=0.2}
			my_element.set(s)
		"""
		#try to merge the two dicts
		try:
			settings.update(dsettings[0])
		except IndexError:
			pass

		for key, val in settings.items():
			self.set_param(key, val)

	def get(self, key):
		"Get a parameter"
		return self._params[key]

	#FIXME investigae if these would be faster as staticmethods
	def f2s(self, f):
		"format float for printing"
		#out = "%e" % float(f)
		#out = "%s" % float(f)
		out = "%.12e" % float(f)
		return out

	def i2s(self, i):
		"format integer for printing"
		out = str(int(i))
		return out

	def l2s(self, l):
		"format label for printing"
		out = l[:max_label_size]
		return out

	def x2s(self, i):
		"format xpas for printing"
		try:
			out = self.f2s(i)
		except TypeError:
			out = '#'+self.i2s(i[0])+ "|"+self.i2s(i[1])+ "|"+self.i2s(i[2])
		return out

	def __getattr__(self, name):
		"allow dot access to parameters"
		if name != '_params':
			try:
				return self._params[name]
			except KeyError:
				pass
		return object.__getattribute__(self, name)

	def list_params(self):
		"Return a list of parameters"
		return list(self._params)

	def reverse(self):
		"Flip the element along the beam line direction, i.e. the entrance and exit properties are swapped"
		if self._zgoubi_name in  ["DIPOLES", "FFAG"]:
			sub_swap_pairs = "G0_E,G0_S KAPPA_E,KAPPA_S NCE,NCS CE_0,CS_0 CE_1,CS_1 CE_2,CS_2 CE_3,CS_3 CE_4,CS_4 CE_5,CS_5 SHIFT_E,SHIFT_S OMEGA_E,OMEGA_S THETA_E,THETA_S R1_E,R1_S U1_E,U1_S U2_E,U2_S R2_E,R2_S"
			for sub_element in self.subelements:
				for swap_pair in sub_swap_pairs.split():
					p1, p2 = swap_pair.split(",")
					sub_element[p1], sub_element[p2] = sub_element[p2], sub_element[p1]
				sub_element["ACN"] = self._params["AT"] - sub_element["ACN"]
				sub_element["OMEGA_E"] *= -1
				sub_element["OMEGA_S"] *= -1
				sub_element["THETA_E"] *= -1
				sub_element["THETA_S"] *= -1
		elif self._zgoubi_name == "CHANGREF":
			self._params["YCE"] *= -1

	def __neg__(self):
		new_e = copy.deepcopy(self)
		new_e.reverse()
		return new_e

	def set_plot_hint(self, **hints):
		"Add hints to help lab_plot"
		if not hasattr(self, "plot_hints"): self.plot_hints = {}
		self.plot_hints.update(hints)

class tXPAS(object):
	"Type for XPAS to handle output"
	def __init__(self, val):
		if isinstance(val, tXPAS):
			val = val.val
		self.val = val

	def __str__(self):
		if isinstance(self.val, Number):
			return "%.12e" % float(self.val)
		else:
			return "#{:d}|{:d}|{:d}".format(*self.val)

param_type_classes = {"E":float, "I":int, "X":tXPAS, "A80": str}

class zgoubi_element_def(zgoubi_element):
	def __init__(self, cdefs, label1='', label2='', **settings):
		self.cdefs = cdefs
		self.label1 = label1
		self.label2 = label2
		self.template = self.cdefs["template_s"]
		self._params = self.cdefs["init_params"].copy()
		self._params_types = self.cdefs["params_types"]
		self.set(settings)

		self.has_subelements = "subelements" in self.cdefs
		if self.has_subelements:
			self.subelements = []
			self.subelement_template = self.cdefs["subelements"][0]["template_s"]
			self.subelement_params = self.cdefs["subelements"][0]["init_params"]
			self.subelement_params_types = self.cdefs["subelements"][0]["params_types"]

	def set_param(self, key, val):
		if key in self._params.keys():
			param_type = param_type_classes[self._params_types[key]]
			self._params[key] = param_type(val)
		else:
			raise ValueError("no such param: '" + str(key) + "' In element " + self._zgoubi_name)

	def add(self, **kwargs):
		if not self.has_subelements:
			raise NotImplemented("Element %s does not have sub-elements"%self._zgoubi_name)
		new_sub_params = self.subelement_params.copy()
		for key, val in  kwargs.items():
			if key in self.subelement_params.keys():
				new_sub_params[key] = val
			else:
				raise ValueError("no such param: '" + str(key) + "' in sub element of " + self._zgoubi_name)
		self.subelements.append(new_sub_params)

	def output(self):
		"Output the element in Zgoubi.dat format"

		# render the conditional section if needed
		cond_out = None
		if "cond_section" in self.cdefs:
			sc_matched = None
			for cs in self.cdefs["cond_section"]:
				if "equals" in cs:
					if self._params[cs["key"]] == cs["equals"]:
						sc_matched = cs
						break
			if sc_matched is None:
				raise ValueError("No matching conditional section")
			cond_out = sc_matched["template_s"].format(**self._params)

		# render sub-element section
		subelement1 = ""
		if self.has_subelements:
			for subelement in self.subelements:
				subelement1 += self.subelement_template.format(**subelement) + "\n"
			subelement1 = subelement1[:-1] # trim  extra newline

		out = "'{zn}' {l1} {l2}\n".format(zn=self._zgoubi_name, l1=self.label1, l2=self.label2)
		out += self.template.format(cond_section=cond_out, subelement1=subelement1, **self._params)
		out += "\n"
		return out

def output_types(s):
	try:
		return dict(I=":d", E=":.12e")[s]
	except KeyError:
		pass
	if s.startswith("A"):
		return ""
	if s == "X":
		#XPAS outputs uings __str__ method of tXPAS class
		return ""
	raise ValueError("Unknown type in definition: %s"%s)

def render_template_string(template, params):
	"Create templates that can be used with string.format"
	template_s = ""
	param_db = params

	for row in template:
		template_row = []
		for item in row:
			if item == "cond_section":
				template_row.append("{%s}"%item)
			elif item.startswith("subelement"):
				template_row.append("{%s}"%item)
			else:
				if not item in param_db:
					raise ValueError("Key used in template, but not defined: %s"%item)
				template_row.append( "{%s%s}"%(item, output_types(param_db[item]["type"])) )
		template_s += " ".join(template_row) +"\n"
	if template_s[-1] == "\n":
		template_s = template_s[:-1]
	return template_s

def make_zgoubi_element(cname, defs):
	def init_func(self,*args, **kargs):
		zgoubi_element_def.__init__(self, defs, *args, **kargs)

	new_class = type(cname, (zgoubi_element_def,),
	            {"__init__":init_func,
	            "_class_name": defs['zgoubi_name'],
	            "_zgoubi_name": defs['zgoubi_name'],
	            "_class_defs": defs})

	globals()[cname] = new_class
	__all__.append(cname)

for e_name, e_def in zgoubi_metadata.elements.get_parsed().items():

	# simple elements just have this template
	if "template" in e_def:
		template_s = render_template_string(e_def["template"], e_def["params"])
		e_def["template_s"] = template_s
		e_def["init_params"] = { _k:_v["default"] for _k, _v in e_def["params"].items() }
		e_def["params_types"] = { _k:_v["type"] for _k, _v in e_def["params"].items() }
	else: # elements like MARKER and END
		e_def["template_s"] = ""
		e_def["init_params"] = {}
		e_def["params_types"] = {}

	# if the element has a conditional section, then each version needs a template
	if "cond_section" in e_def:
		for cs in e_def["cond_section"]:
			template_s = render_template_string(cs["template"], e_def["params"])
			cs["template_s"] = template_s

	if "subelements" in e_def:
		for se in e_def["subelements"]:
			template_s = render_template_string(se["template"], se["params"])
			se["template_s"] = template_s
			se["init_params"] = { _k:_v["default"] for _k,_v in se["params"].items() }
			se["params_types"] = { _k:_v["type"] for _k,_v in se["params"].items() }

	make_zgoubi_element(e_name, e_def)


import zgoubi.static_defs
from zgoubi.static_defs import *
__all__ += zgoubi.static_defs.__all__
