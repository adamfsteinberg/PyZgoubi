from zgoubi.test_macros import *


b = BEND()
assert(b.XL == 0)


b = BEND(XL=1)
assert(b.XL == 1)
assert(b.get("XL") == 1)

b.set(XL=2)
assert(b.XL == 2)
assert(b.get("XL") == 2)

# FIXME should also make this work properly
#b.XL = 3
#assert(b.XL == 3)
#print b.get("XL")
#assert(b.get("XL") == 3)


b = BEND(XL=4, IL=1)
assert(b.XL == 4)
assert(b.IL == 1)



#check some things that should not work
assert_eval_raises("b = BEND(XX=1)", locals())
assert_eval_raises("b = BEND(XX=1)", locals(), NameError)

assert_eval_raises("b = BEND(XL=1, XX=1)", locals(), NameError)


# and the more complex 
d = DIPOLES()
d = DIPOLES(AT=5)
d = DIPOLES(IL=1, AT=5)
assert(d.AT == 5)
d.set(AT=6)
assert(d.IL == 1)
assert(d.AT == 6)

d.add(ACN=1)
d.add(ACN=1, B_0=2)
assert_eval_raises("d.add(XX=2)", locals(), ValueError)
assert_eval_raises("d.add(ACN=1, XX=2)", locals(), ValueError)



