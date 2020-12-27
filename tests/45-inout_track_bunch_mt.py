from numpy.lib.recfunctions import structured_to_unstructured

mass = PROTON_MASS
energy = 1e6


b_orig = Bunch.gen_halo_x_xp_y_yp(1e4, 1e-3, 1e-3, 4, 5, 1e-3, 2e-2, ke=energy, mass=mass, charge=1)
b_orig.particles()[0]['D'] = 1.1

b_orig_4d =  structured_to_unstructured(b_orig.particles()[['Y', 'T', 'Z', 'P']])


line_seg = Line("lineseg")

line_seg.add(PROTON())
line_seg.add(QUADRUPO(XL=1, R_0=5, B_0=0.1, XPAS=(50,500,50)))


import time
t0 = time.time()
st_bunch = line_seg.track_bunch_mt(b_orig, n_threads=1, max_particles=1e3)
t1 = time.time()
st_time = t1-t0

t0 = time.time()
mt2_bunch = line_seg.track_bunch_mt(b_orig, n_threads=2, max_particles=1e3)
t1 = time.time()
mt2_time = t1-t0


t0 = time.time()
mt4_bunch = line_seg.track_bunch_mt(b_orig, n_threads=4, max_particles=1e3)
t1 = time.time()
mt4_time = t1-t0

#print t_bunch.particles()[0]

st_end = structured_to_unstructured(st_bunch.particles()[['Y', 'P', 'Z', 'T', 'D']])
mt2_end = structured_to_unstructured(mt2_bunch.particles()[['Y', 'P', 'Z', 'T', 'D']])
mt4_end = structured_to_unstructured(mt4_bunch.particles()[['Y', 'P', 'Z', 'T', 'D']])


print("%r" % st_end[0])
print("%r" % mt2_end[0])
print("%r" % mt4_end[0])
print()
print("%r" % st_end[1])
print("%r" % mt2_end[1])
print("%r" % mt4_end[1])

errors = abs((mt2_end - st_end) / numpy.maximum(mt2_end, st_end))
print("%r" % errors[0][0])
print("mean errors in YTZPD: single vs 2 thread")
print(errors.mean(0))
assert(numpy.all(errors.mean(0) < [1e-16, 2e-16, 1e-16, 2e-16, 1e-16])), "error to big"


errors = abs((mt4_end - st_end) / numpy.maximum(mt4_end, st_end))
print("%r" % errors[0][0])
print("mean errors in YTZPD: signle vs 4 thread")
print(errors.mean(0))
assert(numpy.all(errors.mean(0) < [1e-16, 2e-16, 1e-16, 2e-16, 1e-16])), "error to big"

print("single thread:", st_time, "s")
print("2 threads    :", mt2_time, "s")
print("4 threads    :", mt4_time, "s")

