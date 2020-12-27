make_line('line')

ob = OBJET2()
ob.set(BORO=ke_to_rigidity(10e6, ELECTRON_MASS))
ob.add(Y=0, T=0.1, D=1)
add(ob)

#test values
tv = [x+0.1234567890123456789 for x in range(1,10)]

q1 = QUADRUPO(XL=tv[0], R_0=tv[1],B_0=tv[2],  XPAS=(10,20,10))
add(q1)
add(END())

output_lines =  output().split('\n')
quad_line = [a.strip() for a in output_lines].index("'QUADRUPO'")
print("Looking for values in pyzgoubi output")
print("  ", output_lines[quad_line])
print("  ", output_lines[quad_line + 1])
print("  >>>", output_lines[quad_line + 2], " <<<")
print("  ", output_lines[quad_line + 3])
print("  ", output_lines[quad_line + 4])
bits = output_lines[quad_line + 2].split()

print()
for x in range(3):
	error = tv[x]- float(bits[x])
	print(repr(tv[x]), bits[x], error)
	if error > 1e-12:
		print("value does not match to 12 decimals")
		raise ValueError("error to big: %s" % error)


run()
res_lines = res().split('\n')


# find QUADRUPO output in res file, skipping over the echo section
for n, l in enumerate(res_lines[quad_line+1:]):
	if "QUADRUPO" in l:
		quad_line = n+quad_line+1
		break

print("Looking for values in Zgoubi output")
print("  ", res_lines[quad_line])
print("  ", res_lines[quad_line + 1])
print("  ", res_lines[quad_line + 2])
print("  ", res_lines[quad_line + 3])
print("  ", res_lines[quad_line + 4])
print("  ", res_lines[quad_line + 5])
print("  ", res_lines[quad_line + 6])
print()

rv=[-1,-1,-1]
for line in res_lines[quad_line:quad_line + 7]:
	if "Length  of  element" in line:
		rv[0] = float(line.partition("=")[2].split()[0])
	if "Bore  radius" in line:
		rv[1] = float(line.partition("=")[2].split()[0])
	if "B-QUADRUPOLE" in line:
		rv[2] = float(line.partition("=")[2].split()[0])

print(rv)

for x in range(3):
	error = tv[x]-rv[x]
	print(repr(tv[x]), rv[x], error)
	if error > 1e-3:
		print("value does not match to 3 decimals")
		raise ValueError("error to big: %s" % error)



