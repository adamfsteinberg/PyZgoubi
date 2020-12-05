
def check(line, types, names):
	el_types = []
	el_names = []
	for el in line.elements():
		el_types.append(el._zgoubi_name)
		el_names.append(el.label1)

	print(line)
	print(el_types)
	print(types)
	assert(el_types == types)
	print(el_names)
	print(names)
	assert(el_names == names)
	print()



cell1 = Line("cell1")
cell2 = Line("cell2")

cell1.add(DRIFT("d1"))
cell1.add(QUADRUPO("q1"))

cell2.add(DRIFT("d2"))
cell2.add(QUADRUPO("q2"))


check(cell1, ["DRIFT","QUADRUPO"], ["d1", "q1"])
check(cell2, ["DRIFT","QUADRUPO"], ["d2", "q2"])

check(-cell1, ["QUADRUPO", "DRIFT"], ["q1", "d1"])


check(2 * cell2, ["DRIFT","QUADRUPO","DRIFT","QUADRUPO"], ["d2","q2","d2","q2"])
check(cell1 + cell2, ["DRIFT","QUADRUPO","DRIFT","QUADRUPO"], ["d1","q1","d2","q2"])

cell3 = Line("cell3")
cell3.add(cell1)
cell3.add(cell2)


check(cell3, ["DRIFT","QUADRUPO","DRIFT","QUADRUPO"], ["d1","q1","d2","q2"])
