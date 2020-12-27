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

q1 = QUADRUPO(XL=1e-6, XPAS=(1,1,1), IL=2, B_0=1e-50, R_0=1)
add(q1)


if binary:
	add(FAISCNL(FNAME='b_zgoubi.fai'))
else:
	add(FAISCNL(FNAME='zgoubi.fai'))

add(END())

#print output()

res = run(xterm=False)

if binary:
	plt_data =  res.get_all('bplt')[['Y','T','Z','P']]
else:
	plt_data =  res.get_all('plt')[['Y','T','Z','P']]
	plt_data["Y"] /= 100
	plt_data["T"] /= 1000
	plt_data["Z"] /= 100
	plt_data["P"] /= 1000


#select the points from entrance of the magnet
plt_data = plt_data[::3]


for key in ['Y','T','Z','P']:
	errors = abs((b_orig_4d[key] - plt_data[key]) / numpy.maximum(b_orig_4d[key], plt_data[key]))
	print("mean errors in YTZP")
	print(errors.mean(0))
	assert(errors.mean(0) < [1e-12]), "error to big for %s"%key

