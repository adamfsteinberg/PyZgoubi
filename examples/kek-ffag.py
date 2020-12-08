#KEK 150 DFD triplet FFAG

ffagex = Line('ffagex')

cells = 12
angle = 360/cells

#Number of FFAG magnets per cell
nffag = 3

#Reference radius of FFAG in meters
rm = 5.4

#Field index
kindex = 7.6

#Field at reference in Tesla in F and D
bz_d = -1.21744691
bz_f = 1.69055873

#Azimuthal positions of FFAG magnets w.r.t left hand boundary
acn1 = 6.465
acn2 = 15.0
acn3 = 23.535

#Fringe field extent [m], shape, coefficients
g0 = 6.3e-2
kappa = 3
nc = 4
c0 = 0.1455
c1 = 2.2670  
c2 = -.6395  
c3 = 1.1558

#azimuth of entrance/exit EFB w.r.t ACN
omega_d = 1.715
omega_f = 5.12

#EFB shaping parameters
r1 = 1e6
u1 = -1e6
u2 = 1e6
r2 = 1e6

#Tracking parameters - 
#KIRD=0 chooses analytic computation of field derivatives
#RESOL when KIRD=0 can be 2 or 4 - order of field derivatives used
#XPAS sets the integration step [cm]
#KPOS allows positioning of magnet. Normally set to 2
kird = 0
resol = 2
xpas = 0.25
kpos = 2


ob = OBJET2()
ffagex.add(ob)

ffagex.add(PROTON())

#Add FFAG
ffag_trip = FFAG('ffag_trip', N =3, AT=angle, RM=rm*cm_, 
KIRD=kird, RESOL=resol, XPAS=xpas, KPOS=2)
ffag_trip.add(ACN = acn1, BZ_0 = bz_d*kgauss_, K = kindex,
              G0_E = g0*cm_, KAPPA_E = kappa, 
              NCE = nc, CE_0 = c0, CE_1 = c1, CE_2 = c2, CE_3 = c3,
              OMEGA_E = omega_d, R1_E = r1, U1_E = u1, U2_E = u2, R2_E = r2,
              G0_S=g0*cm_, KAPPA_S=kappa, 
              NCS = nc, CS_0 = c0, CS_1 = c1, CS_2 = c2, CS_3 = c3,
              OMEGA_S = -omega_d, R1_S = r1, U1_S = u1, U2_S = u2, R2_S = r2,
              KAPPA_L = -1)             
ffag_trip.add(ACN  =  acn2, BZ_0 = bz_f*kgauss_, K = kindex,
              G0_E = g0*cm_, KAPPA_E = kappa, 
              NCE = nc, CE_0 = c0, CE_1 = c1, CE_2 = c2, CE_3 = c3,
              OMEGA_E = omega_f, R1_E = r1, U1_E = u1, U2_E = u2, R2_E = r2,
              G0_S = g0*cm_, KAPPA_S = kappa, 
              NCS = nc, CS_0 = c0, CS_1 = c1, CS_2 = c2, CS_3 = c3,
              OMEGA_S = -omega_f, R1_S = r1, U1_S = u1, U2_S = u2, R2_S = r2,
              KAPPA_L = -1)  
ffag_trip.add(ACN  =  acn3, BZ_0 = bz_d*kgauss_, K = kindex,
              G0_E = g0*cm_, KAPPA_E = kappa, 
              NCE = nc, CE_0 = c0, CE_1 = c1, CE_2 = c2, CE_3 = c3,
              OMEGA_E = omega_d, R1_E = r1, U1_E = u1, U2_E = u2, R2_E = r2,
              G0_S = g0*cm_, KAPPA_S = kappa, 
              NCS = nc, CS_0 = c0, CS_1 = c1, CS_2 = c2, CS_3 = c3,
              OMEGA_S = -omega_d, R1_S = r1, U1_S = u1, U2_S = u2, R2_S = r2,
              KAPPA_L = -1)
ffagex.add(ffag_trip)

ffagex.add(FAISCNL(FNAME='zgoubi.fai'))

#add REBELOTE like this to allow a later modification 
reb=REBELOTE(K=99, NPASS=10)
ffagex.add(reb)

ffagex.add(END())

rigidity = ke_to_rigidity(150e6, 938.272013e6)
ob.set(BORO=rigidity)
ob.add(Y=517, T=0, D=1)

print(ffagex.output())
res = ffagex.run(xterm = False)
print(res.res())
res.clean()

#must give reasonable initial guess at closed orbit 
closedorb_YTZP = find_closed_orbit(ffagex, init_YTZP=[517,0,0,0], D=1.0, tol=1e-6)

#track along closed orbit for one turn
ob.clear()
ob.add(Y=closedorb_YTZP[0],T=closedorb_YTZP[1],Z=closedorb_YTZP[2],P=closedorb_YTZP[3],DR=1.0)
reb.set(NPASS = 1)
ffagex.full_tracking(True)
r = ffagex.run(xterm = False)
s_co = list(numpy.transpose(r.get_track('plt', ['S']))[0])
y_co = list(numpy.transpose(r.get_track('plt', ['Y']))[0])

#Change from OBJET2 to OBJET5 so that MATRIX can compute transfer matrix etc
ob5 = OBJET5()
ffagex.replace(ob, ob5)
ob5.set(BORO=rigidity)
ob5.set(PY=1e-4,PT=1e-3,PZ=1e-4,PP=1e-3,PX=1e-3,PD=1e-3)
ob5.set(YR=closedorb_YTZP[0],TR=closedorb_YTZP[1],ZR=closedorb_YTZP[2],PR=closedorb_YTZP[3],DR=1.0)
matrix=MATRIX(IORD=1,IFOC=11)
ffagex.replace(reb,matrix)
#run zgoubi to find tune etc.
r = ffagex.run(xterm = False)


#find tune calculated by MATRIX over this periodic cell
tune = r.get_tune()
print("tune ",tune)

#get twiss parameters at end of cell, returns [beta_y,alpha_y,gamma_y,disp_y,disp_py,beta_z,alpha_z,gamma_z,disp_z,disp_pz]
twissparam = r.get_twiss_parameters()
betayz = [twissparam['beta_y'][0],twissparam['beta_z'][0]]
alphayz = [twissparam['alpha_y'][0],twissparam['alpha_z'][0]]
gammayz = [twissparam['gamma_y'][0],twissparam['gamma_z'][0]]

#get periodic twiss parameters and dispersion
ffagex.full_tracking(True)
twiss_profiles = get_twiss_profiles(ffagex)


#extract s coordinate and beta_y from twiss_profiles 
s = twiss_profiles['s']
beta_y = twiss_profiles['beta_y']
beta_z = twiss_profiles['beta_z']


import pylab as plt
plt.subplot(211)
plt.plot([t*cm for t in s_co], [t*cm for t in y_co],'k-')
plt.xlim(cm*min(s_co),cm*max(s_co))
plt.ylabel('closed orbit [m]')
plt.title('150 MeV in KEK DFD FFAG')
plt.subplot(212)
plt.plot(s, beta_y, 'r-', label='horizontal')
plt.plot(s, beta_z, 'b-', label='vertical')
plt.xlabel('s [m]')
plt.ylabel('beta function [m]')
plt.ylim(ymin=0)
plt.xlim(min(s),max(s))
plt.legend(loc='lower center')
plt.show()