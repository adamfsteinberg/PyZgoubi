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


ob = OBJET2()
emma.add(ob)

emma.add(ELECTRON())

emma.add(DRIFT('ld', XL=ld*cm_/2))
emma.add(CHANGREF(ALE=angle))

emma.add(CHANGREF(YCE=d_offset*cm_))
emma.add(QUADRUPO('defoc', XL=dq*cm_, R_0=dr*cm_, B_0=db*kgauss_, XPAS=xpas, KPOS=1))
emma.add(CHANGREF(YCE=-d_offset*cm_))

emma.add(DRIFT('sd', XL=sd*cm_))

emma.add(CHANGREF(YCE=f_offset*cm_))
emma.add(QUADRUPO('foc', XL=fq*cm_, R_0=fr*cm_, B_0=fb*kgauss_, XPAS=xpas, KPOS=1))
emma.add(CHANGREF(YCE=-f_offset*cm_))

emma.add(DRIFT('ld', XL=ld*cm_/2))

emma.add(FAISCNL(FNAME='zgoubi.fai'))

emma.add(REBELOTE(K=99, NPASS=10))

emma.add(END())

rigidity = ke_to_rigidity(10e6, 0.51099892e6)
ob.set(BORO=-rigidity)
ob.add(Y=0, T=0, D=1)

#print emma.output()
#emma.run(xterm = False)
#print emma.res()
#emma.clean()

find_closed_orbit(emma, init_YTZP=[0,0,0,0], tol=1e-6)
