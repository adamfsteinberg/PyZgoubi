# -*- coding: utf-8 -*-
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
emma.add(QUADRUPO('defoc', XL=dq*cm_, R_0=dr*cm_, B_0=db*kgauss_, XPAS=xpas, IL=2,KPOS=1))
emma.add(CHANGREF(YCE=-d_offset*cm_))

emma.add(DRIFT('sd', XL=sd*cm_))

emma.add(CHANGREF(YCE=f_offset*cm_))
emma.add(QUADRUPO('foc', XL=fq*cm_, R_0=fr*cm_, B_0=fb*kgauss_, XPAS=xpas, IL=2, KPOS=1))
emma.add(CHANGREF(YCE=-f_offset*cm_))

emma.add(DRIFT('ld', XL=ld*cm_/2))

emma.add(FAISCNL(FNAME='zgoubi.fai'))

#add REBELOTE like this to allow a later modification 
reb=REBELOTE(K=99, NPASS=20)
emma.add(reb)


emma.add(END())

ke = 10e6 #Kinetic energy 10MeV
rigidity = ke_to_rigidity(ke, 0.51099892e6)
ob.set(BORO=-rigidity)

#first step is to find closed orbit
closedorb_YTZP = find_closed_orbit(emma, init_YTZP=[0,0,0,0], tol=1e-6)

#find tof, used for phase slip calculation
ob.add(Y=closedorb_YTZP[0],T=closedorb_YTZP[1],Z=closedorb_YTZP[2],P=closedorb_YTZP[3])
r = emma.run(xterm = False)
tof_ref = r.get_track('fai', ['tof'])[0][0] #TOF in ms

gamma_lorentz = ke_to_gamma(ELECTRON_MASS, ke)

#find longitudinal parameters - phase slip, momentum compaction factor and transition gamma
phase_slip = calc_phase_slip(emma, tof_ref)
alpha_c, gamma_t = calc_momentum_compaction(phase_slip, gamma_lorentz)

#Change from OBJET2 to OBJET5 so that MATRIX can compute transfer matrix etc
ob5 = OBJET5()
emma.replace(ob, ob5)
ob5.set(BORO=-rigidity)
ob5.set(PY=1e-4,PT=1e-3,PZ=1e-4,PP=1e-3,PX=1e-3,PD=1e-3)
ob5.set(YR=closedorb_YTZP[0],TR=closedorb_YTZP[1],ZR=closedorb_YTZP[2],PR=closedorb_YTZP[3],DR=1.0)
matrix=MATRIX(IORD=1,IFOC=11)
emma.replace(reb,matrix)
#run zgoubi to find tune etc.
r = emma.run(xterm = False)

#find tune calculated by MATRIX over this periodic cell
tune = r.get_tune()

#get twiss parameters at end of cell, returns [beta_y,alpha_y,gamma_y,disp_y,disp_py,beta_z,alpha_z,gamma_z,disp_z,disp_pz]
twissparam = r.get_twiss_parameters()
betayz = [twissparam['beta_y'][0],twissparam['beta_z'][0]]
alphayz = [twissparam['alpha_y'][0],twissparam['alpha_z'][0]]
gammayz = [twissparam['gamma_y'][0],twissparam['gamma_z'][0]]

#get periodic twiss parameters and dispersion
emma.full_tracking(True)

twiss_profiles = get_twiss_profiles(emma,'twiss_profiles.txt')

#Note - Could alternatively specify twiss parameters at beginning of cell
#twiss_profiles = get_twiss_profiles(emma,'twiss_profiles.txt',input_twiss_parameters = twissparam)

print("phase_slip ",phase_slip)
print("momentum compaction factor ",alpha_c)
print("beta Y,Z at end of cell ",betayz)
print("alpha Y,Z at end of cell ",alphayz)
print("gamma Y,Z at end of cell ",gammayz)


#extract s coordinate and beta_y from twiss_profiles and create figure beta_y_profile.png
s = twiss_profiles['s']
beta_y = twiss_profiles['beta_y']
beta_z = twiss_profiles['beta_z']
alpha_y = twiss_profiles['alpha_y']
alpha_z = twiss_profiles['alpha_z']

plot_data_xy_multi(s,[beta_y,beta_z], 'beta_profiles', labels=["beta profile","s [m]","beta [m]"],style=['k+','b.'],legend=['beta_y','beta_z'])
plot_data_xy_multi(s,[alpha_y,alpha_z], 'alpha_profiles', labels=["alpha profile","s [m]","alpha [m]"],style=['k+','b.'],legend=['alpha_y','alpha_z'])

#plot horizontal phase advance mu_y
mu_y = twiss_profiles['mu_y']
plot_data_xy_multi(s, mu_y, 'mu_y_profile', labels=["mu_y profile","s [m]","mu_y [rad]"],style='b+')

#plot dispersion in the horizontal plane
disp_y = twiss_profiles['disp_y']
plot_data_xy_multi(s, disp_y, 'disp_profile', labels=["Dispersion profile","s [m]","D [m]"],style='b+')
