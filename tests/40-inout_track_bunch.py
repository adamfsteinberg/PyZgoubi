
mass = PROTON_MASS
energy = 1e6


b_orig = Bunch.gen_halo_x_xp_y_yp(1e2, 1e-3, 1e-3, 4, 5, 1e-3, 2e-2, ke=energy, mass=mass, charge=1)
b_orig.particles()[0]['D'] = 1.1

b_orig_4d = numpy.column_stack([b_orig.particles()[col] for col in 'YTZP'])



line_seg = Line("lineseg")

line_seg.add(PROTON())

t_bunch = b_orig
for x in range(10):
	t_bunch = line_seg.track_bunch(t_bunch)


#print t_bunch.particles()[0]

try:
	# Numpy 1.16 changes how a view() interacts with padding, so use safer structured_to_unstructured
	from numpy.lib.recfunctions import structured_to_unstructured
	start = structured_to_unstructured(b_orig.particles()[['Y', 'P', 'Z', 'T', 'D']])
	end = structured_to_unstructured(t_bunch.particles()[['Y', 'P', 'Z', 'T', 'D']])
except ImportError:
	start = b_orig.particles()[['Y', 'P', 'Z', 'T', 'D']].view(float).reshape([-1, 5])
	end = t_bunch.particles()[['Y', 'P', 'Z', 'T', 'D']].view(float).reshape([-1, 5])


for n, c in enumerate(zip(start[:10], end[:10])):
	print("%s %s -> %s"%(n, list(c[0]), list(c[1])))
	print(abs((c[0] - c[1]) / numpy.maximum(c[0], c[1])))

errors = abs((start - end) / numpy.maximum(start, end))
print("%r" % errors[0][0])
print("mean errors in YTZP")
print(errors.mean(0))
assert(numpy.all(errors.mean(0) < [1e-16, 2e-16, 1e-16, 2e-16, 1e-16])), "error to big"



