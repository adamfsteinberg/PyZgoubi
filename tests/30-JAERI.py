# -*- coding: utf-8 -*-
def are_close(a,b,tol):
	if abs(a) < tol or abs(b) < tol:
		if abs(a-b) < tol:
			return True
	else:
		if (abs(a-b)/a) < tol:
			return True
	return False


head = Line("head")
ob = OBJET2(BORO=ke_to_rigidity(35e6, PROTON_MASS))
head.add(ob)
head.add(PROTON())


xpas = (10,20,10)

cell = Line("cell")
d1 = DRIFT('d1', XL=70)
d2 = DRIFT('d2', XL=275)
d3 = DRIFT('d3', XL=450)
d4 = DRIFT('d4', XL=500)
db = DRIFT('db', XL=150)
dch1 = DRIFT('dch1', XL=20)
dch2 = DRIFT('dch2', XL=10)


bore_rad = 20 * cm
qf1_b0 = 0.475020 * bore_rad * T # field at tip, T
qf1 = QUADRUPO('qf1', XL=25, B_0=qf1_b0 * kgauss_, R_0=bore_rad * cm_, XPAS=xpas, KPOS=1)

qd1_b0 = -0.497895 * bore_rad * T # field at tip, T
qd1 = QUADRUPO('qd1', XL=50, B_0=qd1_b0 * kgauss_, R_0=bore_rad * cm_, XPAS=xpas, KPOS=1)

qf2_b0 = 0.540764 * bore_rad * T # field at tip, T
qf2 = QUADRUPO('qf2', XL=25, B_0=qf2_b0 * kgauss_, R_0=bore_rad * cm_, XPAS=xpas, KPOS=1)

qd2_b0 = -0.642090 * bore_rad * T # field at tip, T
qd2 = QUADRUPO('qd2', XL=50, B_0=qd2_b0 * kgauss_, R_0=bore_rad * cm_, XPAS=xpas, KPOS=1)

qf3_b0 = 0.474652 * bore_rad * T # field at tip, T
qf3 = QUADRUPO('qf3', XL=25, B_0=qf3_b0 * kgauss_, R_0=bore_rad * cm_, XPAS=xpas, KPOS=1)

qd3_b0 = -0.499856 * bore_rad * T # field at tip, T
qd3 = QUADRUPO('qd3', XL=50, B_0=qd3_b0 * kgauss_, R_0=bore_rad * cm_, XPAS=xpas, KPOS=1)

angle = 2*pi/16
hbf_xl = 300 * 2* sin(angle/2) / angle
hbf_b1 =   2 * ke_to_rigidity(35e6, PROTON_MASS) / hbf_xl * sin(angle /2)
hbf = BEND('hbf', XL=hbf_xl, B1= hbf_b1, X_E=0, X_S=0, XPAS=hbf_xl/20.0, KPOS=3)

arc1 = Line('arc1')
arc1.add(qf1,d1,hbf,d1,qd1,d1,hbf,d1,qf2)

arc2 = Line('arc2')
arc2.add(qf2,d1,db,d1,qd2,d1,db,d1,qf2)

arc3 = Line('arc3')
arc3.add(qf2,d1,hbf,d1,qd1,d1,hbf,d1,qf1)

arc4 = Line('arc4')
arc4.add(qf1,d3,qd3,d3,qf3)

quarter_ring = Line('quarter_ring')
quarter_ring.add(d3,qf3,-arc4,arc1,arc2,arc3,qf1,d3,qd3)

print(quarter_ring)


twissline = Line('twissline')
twissline.add(OBJET5(BORO=ke_to_rigidity(35e6, PROTON_MASS),DR=1, PY=1e-4,PT=1e-3,PZ=1e-4,PP=1e-3,PX=1e-3,PD=1e-3))
for e in quarter_ring.elements():
	twissline.add(e)
	twissline.add(FAISCNL(FNAME='zgoubi.fai'))
twissline.add(MATRIX(IORD=1,IFOC=11))
twissline.add(END())
twissline.full_tracking(False)
r = twissline.run(xterm=False)
r.parse_matrix()

tune = r.get_tune()
print("tune", tune)
assert are_close( tune[0], 0.20999960, 1e-6 ) # zgoubi 261 values
assert are_close( tune[1], 0.19500105, 1e-6 ) # zgoubi 261 values

twissparam = r.get_twiss_parameters()
print("twiss params", twissparam)

twiss_profiles = get_twiss_profiles(twissline,'twiss_profiles.txt')


chris_data_raw = """ D3       4.500000  16.320573  -2.255745   2.717455   0.464487  
QF3      4.750000  16.891063  -0.000000   2.602646  -0.000000  
QF3      5.000000  16.320573   2.255745   2.717455  -0.464487  
D3       9.500000   3.573126   0.577021  15.957381  -2.477719  
QD3     10.000000   3.579598  -0.590586  16.107473   2.192167  
D3      14.500000  16.525063  -2.286184   3.676651   0.570238  
QF1     14.750000  17.103256   0.000000   3.535704  -0.000000  
QF1     15.000000  16.525063   2.286185   3.676651  -0.570238  
D1      15.700000  13.509036   2.022425   4.651594  -0.822538  
HBF     18.700000   3.826173   1.037552  12.830683  -1.903825  
D1      19.400000   2.639529   0.657653  15.672648  -2.156125  
QD1     19.900000   2.458573  -0.278502  15.558750   2.372859  
D1      20.600000   3.063237  -0.585303  12.445564   2.074550  
HBF     23.600000   9.522331  -1.455899   3.833673   0.796081  
D1      24.300000  11.721119  -1.685228   2.927977   0.497771  
QF2     24.550000  12.108927   0.154305   2.815727  -0.042924  
QF2     24.800000  11.570812   1.969974   2.972031  -0.590433  
D1      25.500000   9.019541   1.674700   4.020983  -0.908070  
DB      27.000000   4.944535   1.041970   7.766169  -1.588721  
D1      27.700000   3.692468   0.746697  10.212725  -1.906358  
QD2     28.200000   3.692468  -0.746697  10.212725   1.906358  
D1      28.900000   4.944535  -1.041970   7.766169   1.588721  
DB      30.400000   9.019541  -1.674700   4.020983   0.908070  
D1      31.100000  11.570812  -1.969974   2.972031   0.590433  
QF2     31.350000  12.108927  -0.154305   2.815727   0.042924  
QF2     31.600000  11.721119   1.685228   2.927977  -0.497771  
D1      32.300000   9.522331   1.455899   3.833673  -0.796081  
HBF     35.300000   3.063237   0.585303  12.445564  -2.074550  
D1      36.000000   2.458573   0.278502  15.558750  -2.372859  
QD1     36.500000   2.639529  -0.657653  15.672648   2.156125  
D1      37.200000   3.826173  -1.037552  12.830683   1.903825  
HBF     40.200000  13.509036  -2.022425   4.651594   0.822538  
D1      40.900000  16.525063  -2.286185   3.676651   0.570238  
QF1     41.150000  17.103256  -0.000000   3.535704   0.000000  
QF1     41.400000  16.525063   2.286184   3.676651  -0.570238  
D3      45.900000   3.579598   0.590586  16.107473  -2.192167  
QD3     46.400000   3.573126  -0.577021  15.957381   2.477719  """

chris_data = [x.split() for x in  chris_data_raw.split('\n')]

zgoubi_data_t = [twiss_profiles['label'],twiss_profiles['s'],twiss_profiles['beta_y'],twiss_profiles['alpha_y'],twiss_profiles['beta_z'],twiss_profiles['alpha_z']]
zgoubi_data = list(zip(*zgoubi_data_t))

cd = chris_data


print("#element\tdistance\tbetah\talphah\tbetav\talphav")
for n, tp in enumerate (chris_data):
	tp = zgoubi_data[n]
	cd = [chris_data[n][0]] + [float(x) for x in chris_data[n][1:]]

	print("%s\t%s\t%s\t%s\t%s\t%s" %  (tp[0],tp[1], tp[2], tp[3],tp[4],tp[5]))
	print("%s\t%s\t%s\t%s\t%s\t%s" % tuple(cd))

	assert are_close( tp[1], cd[1],1e-8 )
	assert are_close( tp[2], cd[2],1e-4 )
	assert are_close( tp[3], cd[3],1e-4 )
	assert are_close( tp[4], cd[4],1e-4 )
	assert are_close( tp[5], cd[5],1e-4 )
	



