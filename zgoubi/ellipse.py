#!/usr/bin/env python
"Algorithm by J. Scott Berg to find the smallest circle that enclose a set of ellipses that lie on the midplane."

from __future__ import division, print_function
from math import sqrt

def ivl_add(l, ivl):
	"""Add the interval ivl=(min, max) to the interval list l
	"""
	i = 0
	while i < len(l):
		if l[i][0] < ivl[0]:
			i = i + 1
		else:
			if ivl[1] == l[i][0]:
				l[i][0] = ivl[0]
			elif ivl[0] == l[i][1]:
				if i+1 < len(l) and ivl[1] == l[i+1][0]:
					l[i][1] = l[i+1][1]
					l.pop(i+1)
				else:
					l[i][1] = ivl[1]
			else:
				l.insert(i, list(ivl))
			return
	l.append(list(ivl))
	
def solee(c1, a1, c2, a2):
	return 0.5*(c1-a1+c2+a2)

def solemm(c1, a1, c2, a2, b2):
	b2a2 = (b2 - a2) * (b2 + a2)
	x = c2 - c1 - a1
	if x >= 0:
		return c2 + (x * b2a2 + b2 * sqrt(b2a2 * (x - a2) * (x + a2))) / \
			(a2 * a2)
	else:
		return c2 + (x - b2) * (x + b2) / \
			(b2 * sqrt((x - a2) * (x + a2) / b2a2) - x)

def solemp(c1, a1, c2, a2, b2):
	b2a2 = (b2 - a2) * (b2 + a2)
	x = c2 - c1 + a1
	if x <= 0:
		return c2 + (x * b2a2 - b2 * sqrt(b2a2 * (x - a2) * (x + a2))) / \
			(a2 * a2)
	else:
		return c2 + (b2 - x) * (b2 + x) / \
			(b2 * sqrt((x - a2) * (x + a2) / b2a2) + x)

def solC4(c1, a1, b1, c2, a2, b2, sqrtplus):
	c21 = c2 - c1
	c212 = c21 * c21
	b12 = b1 * b1
	b22 = b2 * b2
	b1a1 = (b1 - a1) * (b1 + a1)
	b2a2 = (b2 - a2) * (b2 + a2)
	b2b1 = (b2 - b1) * (b2 + b1)
	badiff = (b2 * a1 - b1 * a2)*(b2 * a1 + b1 * a2)
	r = sqrt(b1a1 * b2a2 * (c212 * b12 * b22 + badiff * b2b1))
	if sqrtplus:
		sign = 1
	else:
		sign = -1
	if (c2 < c1) == sqrtplus:
		return c1 + (-b1a1 * b22 * c21 + sign * r) / badiff
	else:
		return c1 + b1a1 * (b22 * c212 + b2a2 * b2b1) / \
			(b22 * b1a1 * c21 + sign * r)

def ielg(c1, a1, c2, a2, b2):
	b22oa2 = b2 * b2 / a2
	c21 = c2 - c1
	z = []
	if a1 < b22oa2:
		if c21 <= a1 + a2 - 2 * b22oa2:
			z.append(solee(c2, a2, c1, a1))
		elif c21 >= -a1 - a2 + 2 * b22oa2:
			z.append(solee(c1, a1, c2, a2))
		else:
			if a1 <= b2:
				if c21 < a1 - a2:
					z.append(solemm(c1, a1, c2, a2, b2))
				if c21 > a2 - a1:
					z.append(solemp(c1, a1, c2, a2, b2))
			else:
				rt = sqrt((b2 * b2 - a2 * a2) * (a1 * a1 - b2 * b2)) / b2
				if c21 < -rt:
					z.append(solemm(c1, a1, c2, a2, b2))
					if c21 > a2 - a1:
						z.append(solemp(c1, a1, c2, a2, b2))
				elif c21 > rt:
					z.append(solemp(c1, a1, c2, a2, b2))
					if c21 < a1 - a2:
						z.append(solemm(c1, a1, c2, a2, b2))
	elif c21 < a2 - a1:
		z.append(solee(c2, a2, c1, a1))
	elif c21 > a1 - a2:
		z.append(solee(c1, a1, c2, a2))
	return z

def ielpp(c1, a1, b1, c2, a2, b2):
	c21 = c2 - c1
	a12 = a1 * a1
	b12 = b1 * b1
	b22 = b2 * b2
	b22oa2 = b22 / a2
	d1 = (b1 - a1) * (b1 + a1) / a1
	z = []
	if b12 > b2 * a1:
		e1 = sqrt((b12 * b12 - b22 * a12) * (b2 - a2) * (b2 + a2)) / (b2 * a1)
	else:
		e1 = 0
	if b1 > b2:
		f = sqrt((b2 * a1 - b1 * a2) * (b1 - b2) * (b2 * a1 + b1 * a2) * \
					 (b1 + b2)) / (b1 * b2)
	else:
		f = 0
	if c21 <= a1 + a2 - 2 * b22oa2:
		z.append(solee(c2, a2, c1, a1))
	elif c21 < a1 - a2:
		if b12 <= b2 * a1 or c21 <= -d1 - e1 or c21 >= -d1 + e1:
			z.append(solemm(c1, a1, c2, a2, b2))
		elif b1 < b2 or abs(c21) > f:
			z.append(solC4(c1, a1, b1, c2, a2, b2, False))
	if c21 >= -a1 - a2 + 2 * b22oa2:
		z.append(solee(c1, a1, c2, a2))
	elif c21 > a2 - a1:
		if b12 <= b2 * a1 or c21 <= d1 - e1 or c21 >= d1 + e1:
			z.append(solemp(c1, a1, c2, a2, b2))
		elif b1 < b2 or abs(c21) > f:
			z.append(solC4(c1, a1, b1, c2, a2, b2, True))
	return z

def intersect_ellipses(e1, e2):
	a1 = e1[0]
	b1 = e1[1]
	c1 = e1[2]
	a2 = e2[0]
	b2 = e2[1]
	c2 = e2[2]
	z = []
	if b1 <= a1 and b2 <= a2:
		c21 = c2 - c1
		if c21 > abs(a2 - a1):
			z.append(solee(c1, a1, c2, a2))
		elif -c21 > abs(a2 - a1):
			z.append(solee(c2, a2, c1, a1))
		return z
	elif b1 <= a1 and b2 > a2:
		return ielg(c1, a1, c2, a2, b2)
	elif b1 > a1 and b2 <= a2:
		return ielg(c2, a2, c1, a1, b1)
	elif b2 * b2 * a1 >= b1 * b1 * a2:
		return ielpp(c1, a1, b1, c2, a2, b2)
	else:
		return ielpp(c2, a2, b2, c1, a1, b1)

def ellipse_radius2(e, z):
	a = e[0]
	b = e[1]
	c = e[2]
	if b > a:
		zm = c - (b - a) * (b + a) / a
		zp = c + (b - a) * (b + a) / a
	else:
		zm = zp = c
	if z <= zm:
		return (c + a - z) * (c + a - z)
	elif z >= zp:
		return (c - a - z) * (c - a - z)
	else:
		return b * b * ((c - z) * (c - z) / ((b - a) * (b + a)) + 1.0)

class BestCircle:
	"""Class to allow one to find a circle enclosing a set of ellipses.
	append((a, b, c)) adds an ellipse to the set of ellipses, where a is the
	half width, b is the half height, and c is the horizontal center.

	radius(z) returns the radius of a circle with horizontal center at z
	which encloses all the ellipses.
	
	get_circle() returns a sequence (z, r) giving the horizontal center z and
	radius r of the smallest circle enclosing all of the ellipses
	"""
	def __init__(self):
		self.l = []
	def append(self, e):
		"""Add an ellipse e to the set of ellipses
		"""
		if len(self.l)==0:
			self.l.append((e, []))
			ivl_add(self.l[0][1], (-float('inf'), float('inf')))
		else:
			l1 = []
			i = 0
			while i < len(self.l):
				z = intersect_ellipses(e, self.l[i][0])
				a1 = e[0]
				c1 = e[2]
				a2 = self.l[i][0][0]
				c2 = self.l[i][0][2]
				li = self.l[i][1]
				if len(z) == 0:
					if c1+a1 >= c2+a2 and c1-a1 <= c2-a2:
						for ivl in li:
							ivl_add(l1, tuple(ivl))
						self.l.pop(i)
					else:
						i = i + 1
				elif len(z) == 1:
					if c1+a1 >= c2+a2 and c1-a1 <= c2-a2:
						for ivl in li:
							ivl_add(l1, tuple(ivl))
						self.l.pop(i)
					elif c1+a1 <= c2+a2 and c1-a1 >= c2-a2:
						i = i + 1
					else:
						j = 0
						while j < len(li):
							if z[0] <= li[j][0] and c1-a1 <= c2-a2 or \
									li[j][1] <= z[0] and c1+a1 >= c2+a2:
								ivl_add(l1, tuple(li[j]))
								li.pop(j)
							elif li[j][0] < z[0] and z[0] < li[j][1]:
								if c1+a1 >= c2+a2:
									ivl_add(l1, (li[j][0], z[0]))
									li[j][0] = z[0]
								else:
									ivl_add(l1, (z[0], li[j][1]))
									li[j][1] = z[0]
								j = j + 1
							else:
								j = j + 1
						if len(li) == 0:
							self.l.pop(i)
						else:
							i = i + 1
				else:
					z1 = z[0]
					z2 = z[1]
					if (z1 > z2):
						z1, z2 = z2, z1
					j = 0
					while j < len(li):
						if z2 <= li[j][0] and c1-a1 <= c2-a2 or \
								li[j][1] <= z1 and c1+a1 >= c2+a2:
							ivl_add(l1, tuple(li[j]))
							li.pop(j)
						elif li[j][0] < z1 and z1 < li[j][1] <= z2:
							if c1+a1 >= c2+a2:
								ivl_add(l1, (li[j][0], z1))
								li[j][0] = z1
							else:
								ivl_add(l1, (z1, li[j][1]))
								li[j][1] = z1
							j = j + 1
						elif li[j][0] < z1 and z2 < li[j][1]:
							if c1+a1 >= c2+a2 and c1-a1 <= c2-a2:
								if z1 < z2:
									ivl_add(l1, (li[j][0], z1))
									ivl_add(l1, (z2, li[j][1]))
									li[j][0] = z1
									li[j][1] = z2
									j = j + 1
								else:
									ivl_add(l1, tuple(li[j]))
									li.pop(j)
							elif z1 < z2:
								ivl_add(l1, (z1, z2))
								li.insert(j, (li[j][0], z1))
								j = j + 1
								li[j][0] = z2
								j = j + 1
						elif z1 <= li[j][0] and li[j][1] <= z2 and \
								c1+a1 < c2+a2:
							ivl_add(l1, tuple(li[j]))
							li.pop(j)
						elif z1 <= li[j][0] < z2 and z2 < li[j][1]:
							if c1-a1 <= c2-a2:
								ivl_add(l1, (z2, li[j][1]))
								li[j][1] = z2
							else:
								ivl_add(l1, (li[j][0], z2))
								li[j][0] = z2
							j = j + 1
						else:
							j = j + 1
					if len(li) == 0:
						self.l.pop(i)
					else:
						i = i + 1
			if len(l1) > 0:
				self.l.append((e, l1))

	def radius(self, z):
		"""Get the radius of a circle, centered at z, enclosing all of the
		ellipses
		"""
		for ei in self.l:
			for ivl in ei[1]:
				if ivl[0] <= z <= ivl[1]:
					return sqrt(ellipse_radius2(ei[0], z))
		return 0

	def get_circle(self):
		"""Return (z,r), where z is the center and r is the radius of the
		circle with the smallest radius that encloses all the ellipses
		"""
		r = float('inf')
		z = 0.0
		for ei in self.l:
			for ivl in ei[1]:
				if ivl[0] != float('inf') and ellipse_radius2(ei[0], ivl[0]) < r:
					z = ivl[0]
					r = ellipse_radius2(ei[0], z)
			a = ei[0][0]
			b = ei[0][1]
			if a >= b and a*a < r or a <= b and b*b < r:
				j = 0
				c = ei[0][2]
				while j < len(ei[1]) and c >= ei[1][j][1]:
					j = j + 1
				if j < len(ei[1]) and ei[1][j][0] < c:
					z = c
					if a >= b:
						r = a*a
					else:
						r = b*b
		return (z, sqrt(r))


