"""Contains functions useful in testing.

This duplicates some functions that you would find built in to the python test libraries. However PyZgoubi's testing is at quite a high level (check that this script works) rather than proper unit tests (check that this function works), so existing libraries are not a great fit.

"""
from __future__ import division, print_function

def assert_eval_raises(code, loc=None, exc="any"):
	"""Check that the code rasies an error. Usually you want to pass the local variables:

		assert_eval_raises("b = BEND(XX=1)", locals())

	If you care that it is a specific error set exc eg:

		assert_eval_raises("b = BEND(XX=1)", locals(), ValueError)

	"""
	if loc is None:
		loc = locals()
	try:
		exec(code, loc)
	except Exception as inst:
		raised_exc = inst
	else:
		raised_exc = None

	
	# no exception, wanted any
	# no exception, wanted specific
	# exception, wanted any  # good
	# exception, wrong
	# exception, correct     # good

	# no exception
	if raised_exc is None:
		raise AssertionError('Code: "%s" did not raise any exception, expected "%s"'%(code, exc))
	# exception
	else:
		# if we are happy to acept any error type
		if exc=="any":
			print("any exception rasied")
			return
		# otherwsie check it is the correct error
		elif isinstance(raised_exc, exc):
			print("correct exception raised")
			return
		else:
			raise AssertionError('Code: "%s" rasied exception %s, expected "%s"'%(code,type(raised_exc) ,exc))

