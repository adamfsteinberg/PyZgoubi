
my_bunch = Bunch.gen_halo_x_xp_y_yp(10, 1e-3, 1e-3, 4, 5, 1e-3, 2e-2, ke=50e6, mass=PROTON_MASS, charge=1)

my_bunch.write_YTZPSD("mybunch.dat")

my_bunch2 = Bunch.read_YTZPSD("mybunch.dat")


print(my_bunch.particles())
print(my_bunch2.particles())

parts1 = my_bunch.particles()
parts2 = my_bunch2.particles()

for name in my_bunch.particles().dtype.names:
	errors = parts2[name] - parts1[name]
	print(errors)
	assert(errors.max()<1e-15)


