binary = False
for particle in ['e', 'p']: 
	for energy in [1e3, 5e5, 1e6, 1e9, 2e9]:
		laps=200
		if energy == 1e3 and particle == 'p':
			laps=20000 # do a long test for 1 set up
		print(particle, energy, "eV")
		make_line('line')
		if particle == 'e':
			mass = ELECTRON_MASS
		if particle =='p':
			mass = PROTON_MASS
		gamma = (energy + mass) / mass
		beta = sqrt(1-1/(gamma**2))

		print(gamma, beta)

		ob = OBJET2()
		ob.set(BORO=-ke_to_rigidity(energy, mass))
		ob.add(D=1)
		add(ob)
		if particle == 'e':
			add(ELECTRON())
		if particle =='p':
			add(PROTON())


		length = 1*m
		d1 = DRIFT(XL=length*cm_, label1="d1")
		add(d1)
		if binary:
			add(FAISCNL(FNAME='b_zgoubi.fai'))
		else:
			add(FAISCNL(FNAME='zgoubi.fai'))
		add(REBELOTE(K=99, NPASS=laps))
		add(END())

		res = run(xterm=False)

		if binary:
			fai_data =  res.get_all('bfai')
		else:
			fai_data =  res.get_all('fai')

		for n,p in enumerate(fai_data):
			#print p['PASS'], p['Y'], p['T'], p['tof']
			#print abs(( abs(p['tof'] /1e6) - abs((n+1)/SPEED_OF_LIGHT/beta) ) / abs(p['tof'] /1e6) ), (p['tof'] /1e6), ((n+1)/SPEED_OF_LIGHT/beta)
			assert p['PASS'] == n+1
			assert p['Y'] == p['T'] == p['Z'] == p['P'] == 0
			assert abs(( abs(p['tof'] /1e6) - abs((n+1)/SPEED_OF_LIGHT/beta) ) / abs(p['tof'] /1e6) ) < 1e-6, "error in time of flight %s"%abs(( (p['tof'] /1e6) - ((n+1)/SPEED_OF_LIGHT/beta) ) / (p['tof'] /1e6) )

