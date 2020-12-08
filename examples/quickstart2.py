myline = Line('line')

# create and OBJET2 and add 1 electron to it
ob = OBJET2()
ob.set(BORO=ke_to_rigidity(10e6, ELECTRON_MASS))
ob.add(Y=0, T=0.1, D=1)
myline.add(ob)

myline.add(FAISCNL(FNAME='zgoubi.fai')) # record point
d1 = DRIFT(XL=50)
q1 = QUADRUPO(XL=20, B_0=5, R_0=10, XPAS=(10,20,10), KPOS=1)
d2 = DRIFT(XL=50)
q2 = QUADRUPO(XL=20, B_0=-5, R_0=10, XPAS=(10,20,10), KPOS=1)
myline.add(d1, q1, d2, q2)
myline.add(FAISCNL(FNAME='zgoubi.fai')) # record point
myline.add(END())

print(myline.output()) # file fed to zgoubi (zgoubi.dat)
myresults = myline.run() # run the line in zgoubi
print(myresults.res()) # show the zgoubi.res file
print(myresults.get_all('fai')) # so the data recorded with FAISCNL


