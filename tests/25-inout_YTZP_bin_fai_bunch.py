binary = False
#binary = True

make_line('line')
mass = PROTON_MASS
energy = 1e6


b_orig = Bunch.gen_halo_x_xp_y_yp(1e2, 1e-3, 1e-3, 4,5,1e-3, 2e-2, ke=energy, mass=mass, charge=1)

b_orig_4d = b_orig.particles()[["Y", "T", "Z", "P"]]


ob = OBJET_bunch(b_orig, binary=True)
add(ob)
add(PROTON())

if binary:
	add(FAISCNL(FNAME='b_zgoubi.fai'))
else:
	add(FAISCNL(FNAME='zgoubi.fai'))

add(END())

#print output()

res = run(xterm=False)
print(res.res())


try:
	if binary:
		fai_data =  res.get_all('bfai')[['Y','T','Z','P']]
	else:
		fai_data =  res.get_all('fai')[['Y','T','Z','P']]
except:
	print("This will fail if you do not have patch from https://sourceforge.net/tracker/?func=detail&aid=3041984&group_id=205776&atid=995005")
	raise

fai_data["Y"] /= 100
fai_data["T"] /= 1000
fai_data["Z"] /= 100
fai_data["P"] /= 1000

print("%r"%b_orig_4d[0][0])
print("%r"%fai_data[0][0])

for key in ['Y','T','Z','P']:
	errors = abs((b_orig_4d[key] - fai_data[key]) / numpy.maximum(b_orig_4d[key], fai_data[key]))
	print("mean errors in YTZP")
	print(errors.mean(0))
	assert(errors.mean(0) < [1e-12]), "error to big for %s"%key
