# -*- coding: utf-8 -*-
print("running emma-tune-diagram example")

#Calculate cell tune of the EMMA FFAG at various momenta and plot on tune diagram including resonance lines up to third order


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

fb = -6.695 * fr * T
db = 4.704 * dr * T

ob2 = OBJET2()
emma.add(ob2)

emma.add(ELECTRON())


emma.add(DRIFT('ld', XL=ld*cm_/2))
emma.add(CHANGREF(ALE=angle))

emma.add(CHANGREF(YCE=d_offset*cm_))
emma.add(QUADRUPO('defoc', XL=dq*cm_, R_0=dr*cm_, B_0=db*kgauss_, XPAS=xpas, IL=0, KPOS=1))
emma.add(CHANGREF(YCE=-d_offset*cm_))

emma.add(DRIFT('sd', XL=sd*cm_))

emma.add(CHANGREF(YCE=f_offset*cm_))
emma.add(QUADRUPO('foc', XL=fq*cm_, R_0=fr*cm_, B_0=fb*kgauss_, XPAS=xpas, IL=0, KPOS=1))
emma.add(CHANGREF(YCE=-f_offset*cm_))

emma.add(DRIFT('ld', XL=ld*cm_/2))

emma.add(FAISCNL(FNAME='zgoubi.fai'))

#add REBELOTE like this to allow a later modification 
reb=REBELOTE(K=99, NPASS=2)
emma.add(reb)


emma.add(END())

rigidity_ref = ke_to_rigidity(10e6, 0.51099892e6)
ob2.set(BORO=-rigidity_ref)
ob2.add(Y=0, T=0, D=1)

r = emma.run(xterm = False)


nmom = 10
D_mom_lower = ke_to_rigidity(10e6, ELECTRON_MASS)/rigidity_ref
D_mom_upper = ke_to_rigidity(20e6, ELECTRON_MASS)/rigidity_ref

#D is p/p0. D_list is the list of D values at which tunes will be calculated
D_list = [D_mom_lower+i*(D_mom_upper-D_mom_lower)/(nmom-1) for i in range(nmom)]


tune_list = []
for D_sel in D_list:
	#first step is to find closed orbit
	closedorb_YTZP = find_closed_orbit(emma, init_YTZP=[0,0,0,0], tol=1e-10, D=D_sel)

	#Change from OBJET2 to OBJET5 so that MATRIX can compute transfer matrix etc
	ob5 = OBJET5()
	emma.replace(ob2, ob5)
	ob5.set(BORO=-rigidity_ref)
	ob5.set(PY=1e-4,PT=1e-3,PZ=1e-4,PP=1e-3,PX=1e-3,PD=1e-3)
	ob5.set(YR=closedorb_YTZP[0],TR=closedorb_YTZP[1],ZR=closedorb_YTZP[2],PR=closedorb_YTZP[3],DR=D_sel)
	matrix=MATRIX(IORD=1,IFOC=11)
	emma.replace(reb,matrix)
	#run zgoubi to find tune etc.
	r = emma.run(xterm = False)

	#find tune calculated by MATRIX over this periodic cell
	tune = r.get_tune()
	tune_list.append(tune)

	#revert to object2
	emma.replace(ob5, ob2)
	emma.replace(matrix,reb)

#required tune_list in the format [[list of horizontal tunes],[list of vertical tunes]], need to take transpose
tune_list_transpose = numpy.transpose(tune_list)
tune_diagram(tune_list_transpose, order=3, xlim=[0,0.5],ylim=[0,0.5])


