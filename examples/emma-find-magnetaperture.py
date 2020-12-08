# -*- coding: utf-8 -*-
#Find magnet aperture assuming a circular beam pipe. The aperture is given by the smallest circle that encloses the beam ellipses in the magnets at various momenta
print("running emma example")

emma = Line('emma')

xpas = (20,20,20)

cells = 42
angle = 360/cells
d_offset = -34.048 * mm
f_offset = -7.514 * mm

#lengths
ld = 210 * mm
sd = 50 * mm

fq = 58.782 * mm
dq = 75.699 * mm

# quad radius
fr = 37 * mm
dr = 53 * mm

#quad field at tip
#fb = -0.403 / 0.055 * 0.037 * T
#db = 0.367 / 0.065 * 0.053 * T

fb = -6.695 * fr * T
db = 4.704 * dr * T

#set particle mass and energy. Calculate reference momentum
mass_ev = ELECTRON_MASS
kinenergy_ev = 15e6
mom_ev = sqrt((kinenergy_ev+mass_ev)**2-mass_ev**2)


qd = QUADRUPO('defoc', XL=dq*cm_, R_0=dr*cm_, B_0=db*kgauss_, XPAS=xpas, KPOS=1)
qf = QUADRUPO('foc', XL=fq*cm_, R_0=fr*cm_, B_0=fb*kgauss_, XPAS=xpas, KPOS=1)

ob = OBJET2()
emma.add(ob)

emma.add(ELECTRON())

emma.add(DRIFT('ld', XL=ld*cm_/2))
emma.add(CHANGREF(ALE=angle))

emma.add(CHANGREF(YCE=d_offset*cm_))
emma.add(qd)
emma.add(CHANGREF(YCE=-d_offset*cm_))

emma.add(DRIFT('sd', XL=sd*cm_))

emma.add(CHANGREF(YCE=f_offset*cm_))
emma.add(qf)
emma.add(CHANGREF(YCE=-f_offset*cm_))

emma.add(DRIFT('ld', XL=ld*cm_/2))

emma.add(FAISCNL(FNAME='zgoubi.fai'))

#add REBELOTE like this to allow a later modification 
reb=REBELOTE(K=99, NPASS=5)
emma.add(reb)

emma.add(END())

rigidity = ke_to_rigidity(kinenergy_ev, mass_ev)
ob.set(BORO=-rigidity)
ob.add(Y=0, T=0, D=1)

#print emma.output()
#emma.run(xterm = False)
#print emma.res()
#emma.clean()

#first step is to find closed orbit at injection, reference and at extraction and find closed orbit trajectories in magnets
closedorbs =  []
Y_co_QD = []
Y_co_QF = []
D_list = [0.666,1.0,1.333]
for D_sel in D_list:
	closedorb_YTZP = find_closed_orbit(emma, init_YTZP=[0.1,0,0,0], tol=1e-6, D=D_sel)
	closedorbs.append(closedorb_YTZP)
	#find YZ everywhere in magnets along closed orbit
	qd.set(IL=2)
	qf.set(IL=2)
	reb.set(NPASS=1)
	ob.clear()	# remove existing particles
	ob.add(Y=closedorb_YTZP[0], T=closedorb_YTZP[1], Z=closedorb_YTZP[2], P=closedorb_YTZP[3], LET='A', D=D_sel)
	r = emma.run(xterm = False)
	label =flatten(r.get_track('plt', ['element_label1']))
	label = [x.strip() for x in label]
	#find first point in F magnet
	QF_start_index = min(find_indices(label,'foc'))
	Y_co = flatten(r.get_track('plt', ['Y']))
	#assign Y data to Y_co_QD and Y_co_QF
	Y_co_QD.append(Y_co[0:QF_start_index])
	Y_co_QF.append(Y_co[QF_start_index:len(Y_co)])
	qd.set(IL=0)
	qf.set(IL=0)
	reb.set(NPASS=5)


#set IL=2 in quadrupoles so that coordinates are written to zgoubi.plt
qd.set(IL=2)
qf.set(IL=2)

#Change from OBJET2 to OBJET5 so that MATRIX can compute transfer matrix etc, necessary for get_twiss_profiles
ob5 = OBJET5()
emma.replace(ob, ob5)
ob5.set(BORO=-rigidity)
ob5.set(PY=1e-4,PT=1e-3,PZ=1e-4,PP=1e-3,PX=1e-3,PD=1e-3)

matrix=MATRIX(IORD=1,IFOC=11)
emma.replace(reb,matrix)

#set normalised emittance of beam in pi m rad. Assume a round beam (equal emittance in horizontal and vertical planes)
emittance = 3e-3

Y_rad_QD = []
Z_rad_QD = []
Y_rad_QF = []
Z_rad_QF = []
#calculate beta profile at each momentum 
for index, closedorb_YTZP in enumerate(closedorbs):
	#calculate relativistic beta*gamma
	ke_mom = mom_to_ke(D_list[index]*mom_ev,mass_ev)
	beta_gamma_rel = ke_to_relativistic_beta_gamma(ke_mom, mass_ev)
	ob5.set(YR=closedorb_YTZP[0],TR=closedorb_YTZP[1],ZR=closedorb_YTZP[2],PR=closedorb_YTZP[3],DR=D_list[index])	
	twiss_profiles = get_twiss_profiles(emma)
	beta_y = twiss_profiles['beta_y']
	beta_z = twiss_profiles['beta_z']
	#convert from beta profile to beam size using assumed emittance. Units converted to cm
	y_rad = [cm_*sqrt(beta*emittance/beta_gamma_rel) for beta in beta_y]
	z_rad = [cm_*sqrt(beta*emittance/beta_gamma_rel) for beta in beta_z]
	Y_rad_QD.append(y_rad[0:QF_start_index])
	Y_rad_QF.append(y_rad[QF_start_index:len(beta_y)])
	Z_rad_QD.append(z_rad[0:QF_start_index])
	Z_rad_QF.append(z_rad[QF_start_index:len(beta_y)])

#set of [Y_rad_QD, Z_rad_QD, Y_co_QD] define ellipses in QD. Use algorithm to find smallest enclosing circle
QD_ellipse_data = numpy.transpose([flatten(Y_rad_QD), flatten(Z_rad_QD),flatten(Y_co_QD)])
QF_ellipse_data = numpy.transpose([flatten(Y_rad_QF), flatten(Z_rad_QF),flatten(Y_co_QF)])

#Use Scott's Ellipse algorithm to find centre and radius of enclosing circle
print("find magnet aperture in D magnet")
get_enclosing_circle(QD_ellipse_data)

print("find magnet aperture in F magnet")
get_enclosing_circle(QF_ellipse_data)
