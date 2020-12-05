emma_cell = Line('emma')
xpas = (10,10,10)
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
#field
fb = -6.695 * fr * T
db = 4.704 * dr * T


emma_cell.add(DRIFT('ld', XL=ld/cm/2))
emma_cell.add(CHANGREF(ALE=angle))
emma_cell.add(CHANGREF(YCE=d_offset/cm))
emma_cell.add(QUADRUPO('defoc', XL=dq/cm, R_0=dr/cm, B_0=db/kgauss, XPAS=xpas, KPOS=1))
emma_cell.add(CHANGREF(YCE=-d_offset/cm))
emma_cell.add(DRIFT('sd', XL=sd/cm))
emma_cell.add(CHANGREF(YCE=f_offset/cm))
emma_cell.add(QUADRUPO('foc', XL=fq/cm, R_0=fr/cm, B_0=fb/kgauss, XPAS=xpas, KPOS=1))
emma_cell.add(CHANGREF(YCE=-f_offset/cm))
emma_cell.add(DRIFT('ld', XL=ld/cm/2))


data = gcp.get_cell_properties(cell=emma_cell, min_ke=10e6, max_ke=20e6, ke_steps=11, particle='e')
print(gcp.cell_properties_table(data, ["KE", "stable", "Y", "T", "NU_Y", "NU_Z"]))

assert(numpy.all(data['KE']/1e6 == [10,11,12,13,14,15,16,17,18,19,20])) # check energies
assert(numpy.all(data['stable'])) # check stable

expected = {}
expected["Y"] = [ 0.45632163,  0.45643818,  0.42289249,  0.35681577,  0.25933353,
        0.13155897, -0.02541233, -0.21050668, -0.4226756 , -0.66089852,
       -0.92418469]

expected["T"] = [-38.17345815, -31.58627656, -25.13740116, -18.82605897,
       -12.65156118,  -6.61314898,  -0.70991473,   5.05923372,
        10.6955845 ,  16.20061784,  21.57599406]

expected["NU_Y"] = [ 0.35508189,  0.30735083,  0.27438171,  0.24933841,  0.22936626,
        0.21294045,  0.19913457,  0.18733871,  0.1771285 ,  0.16819642,
        0.1603126 ]
        
expected["NU_Z"] = [ 0.27119052,  0.238363  ,  0.21307242,  0.19273991,  0.17592354,
        0.16172517,  0.14954476,  0.13896095,  0.12966663,  0.12143133,
        0.11407792]

#print data['T'].__repr__()

for prop in ["Y", "T", "NU_Y", "NU_Z"]:
	print("Checking", prop)
	max_error = (data[prop] - expected[prop]).max()
	print("  max error", max_error)
	assert ( max_error < 1e-8)
