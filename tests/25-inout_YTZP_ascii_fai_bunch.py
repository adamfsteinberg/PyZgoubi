binary = False
#binary = True

make_line('line')
mass = PROTON_MASS
energy = 1e6


b_orig = Bunch.gen_halo_x_xp_y_yp(1e2, 1e-3, 1e-3, 4,5,1e-3, 2e-2, ke=energy, mass=mass, charge=1)

b_orig_4d = numpy.column_stack([b_orig.particles()[col] for col in 'YTZP'])


#ob = OBJET2()
#ob.set(BORO=-ke_to_rigidity(energy, mass))
#for p in b_orig.particles():
#	ob.add(D=1, Y=p['Y']*100, T=p['T']*1000, Z=p['Z']*100, P=p['P']*1000)

ob = OBJET_bunch(b_orig, binary=False)

add(ob)

add(PROTON())

#length = 1e-15*m
#d1 = DRIFT(XL=length*cm_, label1="d1")
#add(d1)


if binary:
	add(FAISCNL(FNAME='b_zgoubi.fai'))
else:
	add(FAISCNL(FNAME='zgoubi.fai'))


add(END())

#print output()

res = run(xterm=False)
print(res.res())

if binary:
	fai_data =  numpy.array(res.get_track('bfai', ['Y','T','Z','P'])) / [100, 1000, 100, 1000]
	print(type(fai_data))
else:
	fai_data =  numpy.array(res.get_track('fai', ['Y','T','Z','P'])) / [100, 1000, 100, 1000]



print("%r"%b_orig_4d[0][0])
print("%r"%fai_data[0][0])

errors = abs((b_orig_4d - fai_data) / numpy.maximum(b_orig_4d, fai_data))
print("%r"%errors[0][0])
print("mean errors in YTZP")
print(errors.mean(0))
assert(numpy.all(errors.mean(0) < [1e-13, 2e-13, 1e-13, 2e-13])), "error to big"



