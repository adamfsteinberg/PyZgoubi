make_line('line')

# create and OBJET2 and add 1 electron to it
ob = OBJET2()
ob.set(BORO=ke_to_rigidity(10e6, ELECTRON_MASS))
ob.add(Y=0, T=0.1, D=1)
add(ob)

add(FAISCNL(FNAME='zgoubi.fai')) # record point
d1 = DRIFT(XL=50)
q1 = QUADRUPO(XL=20, B_0=5, R_0=10, XPAS=(10,20,10))
d2 = DRIFT(XL=50)
q2 = QUADRUPO(XL=20, B_0=-5, R_0=10, XPAS=(10,20,10))
add(d1, q1, d2, q2)
assert d1.XL == 50, "can't set and read parameters"
add(FAISCNL(FNAME='zgoubi.fai')) # record point
add(END())

out_text =  output() # file fed to zgoubi (zgoubi.dat)
run() # run the line in zgoubi
assert run_success(), "zgoubi did not complete sucessfully"
res_text = res() # show the zgoubi.res file

fai_data =  get_all('fai') # so the data recorded with FAISCNL





print("basic test suscessful")
