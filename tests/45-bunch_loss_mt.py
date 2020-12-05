
inbunch = Bunch(nparticles=100, ke=100e6, mass=PROTON_MASS, charge=1)

inbunch.particles()['Y'] = numpy.linspace(0,1,len(inbunch))


test_line = Line("t")

test_line.add(CHAMBR(IA=1, IFORM=1, ZL=1, YL=50))
test_line.add(DRIFT(XL=1))

print("These should raise some warnings")
outbunch = test_line.track_bunch_mt(inbunch)
assert (len(outbunch) == 50)
