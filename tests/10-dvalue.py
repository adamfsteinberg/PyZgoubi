"""Check that D values are being correctly passed through OBJET and MC OBJET
Requires new style fai headers, in Zgoubi SVN from around r251
"""

def are_close(a,b,tol):
	if abs(a) < tol or abs(b) < tol:
		if abs(a-b) < tol:
			return True
	else:
		if (abs(a-b)/a) < tol:
			return True 
	return False



for dval in [1, 2, 2.3, 5, 0.5]:
	myline = Line('line')

	ob = OBJET2()
	ob.set(BORO=ke_to_rigidity(10e6, ELECTRON_MASS))
	ob.add(Y=0, T=0.1, D=dval)
	myline.add(ob)
	myline.add(FAISCNL(FNAME='zgoubi.fai'))
	myline.add(END())


	myresults = myline.run()

	data = myresults.get_all('fai')

	print(data['D0-1'][0], dval - 1)
	assert are_close(data['D0-1'][0],  dval - 1, 1e-15 )
	assert are_close(data['D-1'][0],  dval - 1, 1e-15 )
	myresults.clean()


	myline2 = Line('line')

	ob2 = MCOBJET3()
	ob2.set(BORO=ke_to_rigidity(10e6, ELECTRON_MASS), IMAX=10)
	ob2.set(KY=3, KT=3, KZ=3, KP=3, KX=2, KD=2)
	ob2.set(Y0=0, T0=0, Z0=0, P0=0, X0=0, D0=dval)
	ob2.set(alpha_x=0.906562, alpha_y=-0.925411, alpha_z=-4)
	ob2.set(beta_x=0.303851, beta_y=0.326992, beta_z=3)
	ob2.set(emit_x=0, emit_y=0, emit_z=2E-4)
	ob2.set(n_cutoff_x=2, n_cutoff_y=2, n_cutoff_z=-0.4, n_cutoff2_z=0.40001)

	myline2.add(ob2)
	myline2.add(FAISCNL(FNAME='zgoubi.fai'))
	myline2.add(END())

	myresults2 = myline2.run()

	data2 = myresults2.get_all('fai')

	print(data2['D0-1'])
	assert data2['D0-1'].size == 10
	assert are_close(data2['D0-1'][0],  dval - 1, 1e-15 )
	assert are_close(data2['D0-1'][8],  dval - 1, 1e-15 )
	assert are_close(data2['D-1'][0],  dval - 1, 1e-15 )
	myresults2.clean()
