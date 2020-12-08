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

fb = -6.695 * fr * T
db = 4.704 * dr * T

ob = OBJET2()
emma.add(ob)

emma.add(ELECTRON())

emma.add(DRIFT('ld', XL=ld*cm_/2))
emma.add(CHANGREF(ALE=angle))

emma.add(CHANGREF(YCE=d_offset*cm_))
emma.add(QUADRUPO('defoc', XL=dq*cm_, R_0=dr*cm_, B_0=db*kgauss_, XPAS=xpas, IL=2, KPOS=1))
emma.add(CHANGREF(YCE=-d_offset*cm_))

emma.add(DRIFT('sd', XL=sd*cm_))

emma.add(CHANGREF(YCE=f_offset*cm_))
emma.add(QUADRUPO('foc', XL=fq*cm_, R_0=fr*cm_, B_0=fb*kgauss_, XPAS=xpas, IL=2, KPOS=1))
emma.add(CHANGREF(YCE=-f_offset*cm_))

emma.add(DRIFT('ld', XL=ld*cm_/2))

emma.add(FAISCNL(FNAME='zgoubi.fai'))

emma.add(REBELOTE(K=99, NPASS=10))

emma.add(END())

rigidity = ke_to_rigidity(10e6, 0.51099892e6)
ob.set(BORO=-rigidity)
ob.add(Y=0, T=0, D=1)

#print emma.output()
res = emma.run()

fai_data =  res.get_all('fai')
print("fai")
print("PASS Y T")
for p in fai_data:
	print(p['PASS'], p['Y'], p['T'])

plt_data = res.get_all('plt')
print("plt")
print("PASS Y T element_label1")
for p in plt_data:
	if p['PASS'] == 3:
		print(p['PASS'], p['Y'], p['T'], p['element_label1'])
		





#print emma.res()
emma.clean()

