binary = False

make_line('line')
mass = PROTON_MASS
energy = 1e6


b_orig = Bunch.gen_halo_x_xp_y_yp(1e3, 1e-3, 1e-3, 4,5,1e-3, 2e-2)

b_orig_4d = b_orig.particles()[["Y", "T", "Z", "P"]]


ob = OBJET2()
ob.set(BORO=-ke_to_rigidity(energy, mass))
for p in b_orig.particles():
	ob.add(D=1, Y=p['Y']*100, T=p['T']*1000, Z=p['Z']*100, P=p['P']*1000)
add(ob)
add(PROTON())


#length = 1e-6*m
#d1 = DRIFT(XL=length*cm_, label1="d1")
#add(d1)

if binary:
	add(FAISCNL(FNAME='b_zgoubi.fai'))
else:
	add(FAISCNL(FNAME='zgoubi.fai'))

add(END())

#print output()

res = run(xterm=False)

if binary:
	fai_data =  res.get_all('bfai')[['Y','T','Z','P']]
else:
	fai_data =  res.get_all('fai')[['Y','T','Z','P']]
	fai_data["Y"] /= 100
	fai_data["T"] /= 1000
	fai_data["Z"] /= 100
	fai_data["P"] /= 1000

for key in ['Y','T','Z','P']:
	errors = abs((b_orig_4d[key] - fai_data[key]) / numpy.maximum(b_orig_4d[key], fai_data[key]))
	print("mean errors in YTZP")
	print(errors.mean(0))
	assert(errors.mean(0) < [1e-12]), "error to big for %s"%key
