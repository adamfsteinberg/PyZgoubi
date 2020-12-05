from zgoubi.lab_plot import LabPlot

print("Drifts and bends")

l = Line("lp test")
for x in range(2):
	l.add(DRIFT("start", XL=0))
	l.add(FAISCNL("start", FNAME='zgoubi.fai'))
	l.add(DRIFT("d1", XL=5))
	l.add(BEND("b1", XPAS=(10,10,10), XL=5, B1=-4, KPOS=3, W_E = radians(0), W_S = radians(0)))
	l.add(DRIFT("d2", XL=5))
	l.add(DRIFT("d2", XL=5))
	l.add(BEND("b2", XPAS=(10,10,10), XL=5, B1=4, KPOS=1, W_E = radians(0), W_S = radians(-20)))
	l.add(DRIFT("d3", XL=5))
	l.add(DRIFT("d4", XL=5))
	l.add(BEND("b3", XPAS=(10,10,10), XL=5, B1=-4, KPOS=1, W_E = radians(0), W_S = radians(0)))
	l.add(DRIFT("d5", XL=5))
	l.add(DRIFT("end",XL=0))
	l.add(FAISCNL("end", FNAME='zgoubi.fai'))

l = uniquify_labels(l)

tline = Line('test_line')
ob = OBJET2()
rigidity = - ke_to_rigidity(10e6,ELECTRON_MASS)
ob.set(BORO=rigidity)
tline.add(ob)
tline.add(ELECTRON())
tline.add(l)
tline.add(REBELOTE(NPASS=0, K=99))
tline.add(END())

for y in numpy.linspace(-10,10,10):
	ob.add(Y=y, T=0, Z=0, P=0, D=1)

tline.full_tracking(drift_to_multi=True)
print(tline)
res= tline.run(xterm=0)
ftrack = res.get_all('fai')
ptrack = res.get_all('plt')

lp = LabPlot(tline, boro=-ke_to_rigidity(10e6,ELECTRON_MASS), style={"track":{"linewidth":1, "color":"g"}})
#lp.add_tracks(ftrack=ftrack)
#lp.add_tracks(ptrack=ptrack)
lp.add_tracks(ftrack, ptrack)
lp.draw()
mkdir_p("plots")
lp.save("plots/41-labplot_1.pdf")


print("Drifts and bends 2")

l = Line("lp test")
for x in range(2):
	l.add(DRIFT("start", XL=0))
	l.add(DRIFT("d1", XL=5))
	l.add(BEND("b1", XPAS=(10,10,10), XL=5, B1=-4, KPOS=3, W_E = radians(0), W_S = radians(0)))
	l.add(DRIFT("d2", XL=5))
	l.add(DRIFT("d2", XL=5))
	l.add(BEND("b2", XPAS=(10,10,10), XL=5, B1=4, KPOS=1, W_E = radians(0), W_S = radians(-20)))
	l.add(DRIFT("d3", XL=5))
	l.add(DRIFT("d4", XL=5))
	l.add(BEND("b3", XPAS=(10,10,10), XL=5, B1=-4, KPOS=1, W_E = radians(0), W_S = radians(0)))
	l.add(DRIFT("d5", XL=5))
	l.add(DRIFT("end",XL=0))

l = uniquify_labels(l)

tline = Line('test_line')
ob = OBJET2()
rigidity = - ke_to_rigidity(10e6,ELECTRON_MASS)
ob.set(BORO=rigidity)
tline.add(ob)
tline.add(ELECTRON())
for e in l.elements():
	tline.add(e)
	tline.add(FAISCNL("f", FNAME='zgoubi.fai'))
tline.add(REBELOTE(NPASS=0, K=99))
tline.add(END())

for y in numpy.linspace(-10,10,10):
	ob.add(Y=y, T=0, Z=0, P=0, D=1)

tline.full_tracking(drift_to_multi=True)
print(tline)
res= tline.run(xterm=0)
ftrack = res.get_all('fai')
ptrack = res.get_all('plt')

lp = LabPlot(tline, boro=-ke_to_rigidity(10e6,ELECTRON_MASS), style={"track":{"linewidth":1, "color":"g"}})
#lp.add_tracks(ftrack=ftrack)
#lp.add_tracks(ptrack=ptrack)
lp.add_tracks(ftrack, ptrack)
lp.draw()
mkdir_p("plots")
lp.save("plots/41-labplot_2.pdf")


print("DIPOLE")

l = Line("lp test")
ref_rid = -ke_to_rigidity(10e6, ELECTRON_MASS)
dipole_angle = 2*pi/36
dipole_field = -0.2*T
dipole_bend_radius  = abs(ref_rid / (dipole_field *kgauss_ ) *cm)
ad = abs(degrees(dipole_angle))
for x in range(2):
	l.add(DRIFT("start", XL=0))
	l.add(FAISCNL("start", FNAME='zgoubi.fai'))
	l.add(DRIFT("d1", XL=5))
	l.add(DIPOLE("b1", AT=ad, RM=dipole_bend_radius*cm_, ACN=ad/2, B_0=dipole_field*kgauss_, N=2,
		OMEGA_E=ad/2, OMEGA_S=-ad/2,
		R1_E=1e9, R2_E=1e9, R1_S=1e9, R2_S=1e9, R1_L=1e9, R2_L=1e9,R3=1e9,
		U1_E=-1e9, U2_E=1e9, U1_S=-1e9, U2_S=1e9, U1_L=-1e9, U2_L=1e9,
		IORDRE=2, Resol=10,
		LAM_E=1,NCE=2,CE_0=0.1,CE_1=2, XI_E=1,SHIFT_E=1,
		LAM_S=0.1,
		KPOS=2, RE=dipole_bend_radius*cm_, RS=dipole_bend_radius*cm_,
		XPAS=1
		))

	if x==1: l.add(FAISCNL("f1", FNAME='zgoubi.fai'))
	l.add(DRIFT("d2", XL=5))
	if x==1: l.add(FAISCNL("f2", FNAME='zgoubi.fai'))
	l.add(DRIFT("d3", XL=5))
	if x==1: l.add(FAISCNL("f3", FNAME='zgoubi.fai'))
	l.add(DIPOLE("b2", AT=ad, RM=dipole_bend_radius*cm_, ACN=ad/2, B_0=-dipole_field*kgauss_, N=-2,
		OMEGA_E=ad/2, OMEGA_S=-ad/2,
		R1_E=1e9, R2_E=1e9, R1_S=1e9, R2_S=1e9, R1_L=1e9, R2_L=1e9,R3=1e9,
		U1_E=-1e9, U2_E=1e9, U1_S=-1e9, U2_S=1e9, U1_L=-1e9, U2_L=1e9,
		IORDRE=2, Resol=10,
		LAM_E=1,NCE=2,CE_0=0.1,CE_1=2, XI_E=1,SHIFT_E=1,
		LAM_S=0.1,
		KPOS=2, RE=dipole_bend_radius*cm_, RS=dipole_bend_radius*cm_,
		XPAS=1
		))

	l.add(DRIFT("d4", XL=5))
	l.add(FAISCNL("end", FNAME='zgoubi.fai'))

l = uniquify_labels(l)

tline = Line('test_line')
ob = OBJET2()
rigidity = - ke_to_rigidity(10e6,ELECTRON_MASS)
ob.set(BORO=rigidity)
tline.add(ob)
tline.add(ELECTRON())
tline.add(l)
tline.add(REBELOTE(NPASS=0, K=99))
tline.add(END())

for y in numpy.linspace(-10,10,10):
	ob.add(Y=y, T=0, Z=0, P=0, D=1)

tline.full_tracking(drift_to_multi=True)
print(tline)
res= tline.run(xterm=0)
ftrack = res.get_all('fai')
ptrack = res.get_all('plt')

lp = LabPlot(tline, boro=-ke_to_rigidity(10e6,ELECTRON_MASS), style={"track":{"linewidth":1, "color":"g"}})
#lp.add_tracks(ftrack=ftrack)
#lp.add_tracks(ptrack=ptrack)
lp.add_tracks(ftrack, ptrack)
lp.draw()
mkdir_p("plots")
lp.save("plots/41-labplot_3.pdf")


print("FFAG")

l = Line("lp test")
ref_rid = -ke_to_rigidity(10e6, ELECTRON_MASS)
dipole_angle = 2*pi/36
dipole_field = -0.2*T
dipole_bend_radius  = abs(ref_rid / (dipole_field *kgauss_ ) *cm)
ad = abs(degrees(dipole_angle))
for x in range(2):
	l.add(DRIFT("start", XL=0))
	l.add(FAISCNL("start", FNAME='zgoubi.fai'))
	l.add(DRIFT("d1", XL=5))
	
	big = 1e6
	c0, c1 = .1455, 2.2670
	rm = 100
	
	il = 2
	
	ffag = FFAG('ff', N =1, AT=20, RM=rm, IL=il,
	             KIRD=0, RESOL=2, XPAS=0.5, KPOS=2, RE=rm, RS=rm)
	ffag.add(ACN = 10, BZ_0 = 0.2, K = 10,
				  G0_E = 2, KAPPA_E = 0, 
				  NCE = 2, CE_0 = c0, CE_1 = c1,
				  OMEGA_E = 5, R1_E=big, U1_E=-big, U2_E=big, R2_E=big,
				  G0_S=2, KAPPA_S=0, 
				  NCS = 2, CS_0 = c0, CS_1 = c1,
				  OMEGA_S = -5, R1_S=big, U1_S=-big, U2_S=big, R2_S=big,
				  KAPPA_L = -1)
	l.add(ffag)

	if x==1: l.add(FAISCNL("f1", FNAME='zgoubi.fai'))
	l.add(DRIFT("d2", XL=5))
	if x==1: l.add(FAISCNL("f2", FNAME='zgoubi.fai'))
	l.add(DRIFT("d3", XL=5))
	if x==1: l.add(FAISCNL("f3", FNAME='zgoubi.fai'))

	ffag = FFAG('ff', N =1, AT=20, RM=rm, IL=il,
	             KIRD=0, RESOL=2, XPAS=0.5, KPOS=2, RE=rm, RS=rm)
	ffag.add(ACN = 10, BZ_0 = -0.2, K = 10,
				  G0_E = 2, KAPPA_E = 0, 
				  NCE = 2, CE_0 = c0, CE_1 = c1,
				  OMEGA_E = 5, R1_E=big, U1_E=-big, U2_E=big, R2_E=big,
				  G0_S=2, KAPPA_S=0, 
				  NCS = 2, CS_0 = c0, CS_1 = c1,
				  OMEGA_S = -5, R1_S=big, U1_S=-big, U2_S=big, R2_S=big,
				  KAPPA_L = -1)
	l.add(ffag)

	l.add(DRIFT("d4", XL=5))
	l.add(FAISCNL("end", FNAME='zgoubi.fai'))

l = uniquify_labels(l)

tline = Line('test_line')
ob = OBJET2()
rigidity = - ke_to_rigidity(10e6,ELECTRON_MASS)
ob.set(BORO=rigidity)
tline.add(ob)
tline.add(ELECTRON())
tline.add(l)
tline.add(REBELOTE(NPASS=0, K=99))
tline.add(END())

for y in numpy.linspace(-10,10,10):
	ob.add(Y=y, T=0, Z=0, P=0, D=1)

#tline.full_tracking(False, drift_to_multi=False)
print(tline)
res= tline.run(xterm=0)
ftrack = res.get_all('fai')
ptrack = res.get_all('plt')

lp = LabPlot(tline, boro=-ke_to_rigidity(10e6,ELECTRON_MASS), style={"track":{"linewidth":1, "color":"g"}})
#lp.add_tracks(ftrack=ftrack)
#lp.add_tracks(ptrack=ptrack)
lp.add_tracks(ftrack, ptrack)
lp.draw()
mkdir_p("plots")
lp.save("plots/41-labplot_4.pdf")


print("FFAG 2")

l = Line("lp test")
ref_rid = -ke_to_rigidity(10e6, ELECTRON_MASS)
dipole_angle = 2*pi/36
dipole_field = -0.2*T
dipole_bend_radius  = abs(ref_rid / (dipole_field *kgauss_ ) *cm)
ad = abs(degrees(dipole_angle))

for x in range(1):
	
	#l.add(MULTIPOL("d", XL=0, B_1=1e-10, R_0=1, KPOS=1))
	#l.add(FAISCNL("", FNAME='zgoubi.fai'))
	big = 1e6
	c0, c1 = .1455, 2.2670
	rm = 100
	
	il = 2
	
	ffag = FFAG('ff', N =1, AT=20, RM=rm, IL=il,
	             KIRD=0, RESOL=2, XPAS=0.5, KPOS=2, RE=rm, RS=rm)
	ffag.add(ACN = 10, BZ_0 = 0.2, K = 10,
				  G0_E = 2, KAPPA_E = 0, 
				  NCE = 2, CE_0 = c0, CE_1 = c1,
				  OMEGA_E = 5, R1_E=big, U1_E=-big, U2_E=big, R2_E=big,
				  G0_S=2, KAPPA_S=0, 
				  NCS = 2, CS_0 = c0, CS_1 = c1,
				  OMEGA_S = -5, R1_S=big, U1_S=-big, U2_S=big, R2_S=big,
				  KAPPA_L = -1)
	l.add(ffag)

	l.add(FAISCNL("", FNAME='zgoubi.fai'))
	l.add(DRIFT("d", XL=5))
	l.add(FAISCNL("", FNAME='zgoubi.fai'))
	#l.add(MULTIPOL("d", XL=0, B_1=1e-10, R_0=1, KPOS=1))
	#l.add(FAISCNL("", FNAME='zgoubi.fai'))

l = uniquify_labels(l)

tline = Line('test_line')
ob = OBJET2()
rigidity =  ke_to_rigidity(10e6,ELECTRON_MASS)
ob.set(BORO=rigidity)
tline.add(ob)
tline.add(ELECTRON())
tline.add(DRIFT("start", XL=0))
tline.add(FAISCNL("", FNAME='zgoubi.fai'))
tline.add(DRIFT("", XL=10))
for e in l.elements():
	tline.add(e)
tline.add(DRIFT("end", XL=10))
tline.add(FAISCNL("", FNAME='zgoubi.fai'))
tline.add(END())

for y in numpy.linspace(-10,10,2):
	ob.add(Y=y, T=0, Z=0, P=0, D=1)

tline.full_tracking(True, drift_to_multi=False)
res= tline.run(xterm=0)
ftrack = res.get_all('fai')
ptrack = res.get_all('plt')

lp = LabPlot(tline, boro=-ke_to_rigidity(10e6,ELECTRON_MASS), style={"track":{"linewidth":1, "color":"g"}})
#lp.add_tracks(ftrack=ftrack)
#lp.add_tracks(ptrack=ptrack)
lp.add_tracks(ftrack, ptrack)
lp.draw()
mkdir_p("plots")
lp.save("plots/41-labplot_5.pdf")
