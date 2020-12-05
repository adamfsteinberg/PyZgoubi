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

ffagex.add(REBELOTE(K=99, NPASS=10))

ffagex.add(END())

rigidity = ke_to_rigidity(150e6, 938.272013e6)
ob.set(BORO=rigidity)
ob.add(Y=517, T=0, D=1)

print(ffagex.output())
res = ffagex.run(xterm = False)
print(res.res())
res.clean()

#must give reasonable initial guess at closed orbit
#co1 = find_closed_orbit(ffagex, init_YTZP=[479.7,0,0,0], D=0.52, tol=1e-10)
co1 = find_closed_orbit(ffagex, init_YTZP=[480,0,0,0], D=0.52, tol=1e-10)
assert (co1 is not None)
#co2 = find_closed_orbit(ffagex, init_YTZP=[517,0,0,0], D=1.0, tol=1e-6)
co2 = find_closed_orbit(ffagex, init_YTZP=[510,0,0,0], D=1.0, tol=1e-10)
assert (co2 is not None)


print(co1)
Y1, T1, Z1, P1 = co1
Y0, T0, Z0, P0 = 4.79700585e+02, -6.44703079e-08, 0.0, 0.0 # from a zgoubi-5.0.0 run on 64bit amd Opteron

assert( abs((Y0-Y1)/Y0) < 1e-7  )
#assert( abs((T0-T1)/T0) < 1e-4  ) # value is 10^-7 mrad, don't need much accuracy
assert( abs(T1) < 1e-6  ) # just be small
assert( abs((Z0-Z1)) < 1e-10  )
assert( abs((P0-P1)) < 1e-10  )



print(co2)
Y1, T1, Z1, P1 = co2
Y0, T0, Z0, P0 = 5.17498257e+02, -1.59120900e-07, 0.0, 0.0 # from a zgoubi-5.0.0 run on 64bit amd Opteron

assert( abs((Y0-Y1)/Y0) < 1e-7  )
#assert( abs((T0-T1)/T0) < 1e-4  )
assert( abs(T1) < 1e-6  ) # just be small
assert( abs((Z0-Z1)) < 1e-10  )
assert( abs((P0-P1)) < 1e-10  )


