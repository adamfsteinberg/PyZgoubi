#!/usr/bin/env python3
from __future__ import print_function
import tempfile
import shutil
import os
import sys
import subprocess
import time
import datetime
import getopt

def usage():
	print("Usage:")
	print(sys.argv[0])
	print(sys.argv[0], "[testname [testname2 ...]]")
	print(sys.argv[0], "--zgoubi=/path/to/zgoubi")
	print(sys.argv[0], "--zgoubi=/path/to/zgoubi1,/path/to/zgoubi2")
	print(sys.argv[0], "--zgoubi=~/bin")
	print(sys.argv[0], "--logfile=testlog.txt")



try:
	opts, args = getopt.getopt(sys.argv[1:], "h", ["help", "logfile=", "zgoubi="])
except getopt.GetoptError as err:
# print help information and exit:
	print(str(err)) # will print something like "option -a not recognized"
	usage()
	sys.exit(2)


orig_dir = os.getcwd()


log_file_path = "test-" + datetime.datetime.today().strftime("%Y%m%d-%H%M") + ".log"
zgoubi_bins = [None]

for o, a in opts:
	if o in ("-h", "--help"):
		usage()
		sys.exit()
	elif o in ("--logfile"):
		log_file_path = a
	elif o in ("--zgoubi"):
		zgoubi_bins = []
		zgoubi_paths = [os.path.expanduser(x) for x in  a.split(',')]
		for path in zgoubi_paths:
			if os.path.isdir(path):
				files = [os.path.join(path, x) for x in os.listdir(path) if x.startswith("zgoubi")]
				zgoubi_bins += files
			else:
				zgoubi_bins.append(path)

		print(zgoubi_bins)

	

log = open(log_file_path, "w")
print("writing test log to:", log_file_path)

#install pyzoubi to temp folder
install_dir = tempfile.mkdtemp(prefix='pyzgoubi_test_inst_')
print("installing to", install_dir)

print("installing to", install_dir, file=log)
#install_res = os.system("python setup.py install --prefix=%s"%install_dir)
#subprocess.Popen(["./setup.py", "clean", "--all"])
clean_proc = subprocess.Popen(["python3", "./setup.py", "clean", "--all"], stdout=log, stderr=subprocess.STDOUT)
clean_proc.wait()
install_res = subprocess.Popen(["python3", "./setup.py", "install", "--single-version-externally-managed", "--prefix=%s"%install_dir], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding="utf-8")

for line in install_res.communicate()[0].split('\n'):
	print(line, file=log)
	if line.startswith('export PYTHONPATH='):
		env_pythonpath = line.rpartition('=')[2].rpartition(':')[0]
		os.environ["PYTHONPATH"] = os.pathsep + env_pythonpath
	if line.startswith('export PATH='):
		env_path = line.rpartition('=')[2].rpartition(':')[0]
		os.environ["PATH"] = os.pathsep + env_path
	if line:
		last_line = line.strip()

if sys.platform == "win32":
	pyzgoubi_cmd = os.path.join(last_line, "pyzgoubi.bat")
else:
	pyzgoubi_cmd = "pyzgoubi"
print(last_line)
	
if install_res.returncode != 0:
	print("ERROR: install failed")
	print("ERROR: install failed", file=log)
	print("If there were permission errors, try running 'sudo ./setup.py clean --all' and 'sudo rm install.log zgoubi/version.py'")
	print("If there were permission errors, try running 'sudo ./setup.py clean --all' and 'sudo rm install.log zgoubi/version.py'", file=log)
	sys.exit(1)

pyzgoubi_cmd += " --debug"

print()
print(file=log)
log.flush()
pyzgoubi_cmd = pyzgoubi_cmd.replace(r"\$", "$")
proc = subprocess.Popen(pyzgoubi_cmd+" --version", shell=True, stderr=subprocess.STDOUT, stdout=log)
proc.wait()

# move to another temp dir for running tests
run_dir = tempfile.mkdtemp(prefix='pyzgoubi_test_run_')
print("running tests from", run_dir)
print("running tests from", run_dir, file=log)
os.chdir(run_dir)

test_dir = os.path.join(install_dir, 'share', 'pyzgoubi','test')

tests = os.listdir(test_dir)
tests.sort()

number_of_tests = len(tests)
tests_run = 0
tests_sucess = []
tests_fail = []
tot_time = 0

for test_file in tests:
	print()
	print("\n", "="*40, file=log)
	log.flush()
	if len(args) > 0:
		if test_file not in args:
			print("skipping", test_file)
			print("skipping", test_file, file=log)
			continue
	full_test_file = os.path.join(test_dir, test_file)
	tests_run += 1
	print("running test %s, %d of %d"%(test_file, tests_run, number_of_tests))
	print("running test %s, %d of %d"%(test_file, tests_run, number_of_tests), file=log)
	log.flush()
	for zgoubi_bin in zgoubi_bins:
		if zgoubi_bin is None:
			command = pyzgoubi_cmd + " " + full_test_file
			test_name = test_file
		else:
			command = pyzgoubi_cmd + " --zgoubi="+ zgoubi_bin + " " + full_test_file
			print("Using zgoubi:", zgoubi_bin)
			print("Using zgoubi:", zgoubi_bin, file=log)
			test_name = test_file + " " + zgoubi_bin

		t0 = time.time()
		#	result = os.system(command)
		print(command)
		print(command, file=log)
		log.flush()
		proc = subprocess.Popen(command, shell=True, stderr=subprocess.STDOUT, stdout=log)
		proc.wait()
		result = proc.returncode
		t1 = time.time()
		t = t1 - t0
		log.flush()
		print()
		print("\n", "="*40, file=log)
		log.flush()
		if result == 0:
			print("PASS:", test_name)
			print("PASS:", test_name, file=log)
			tests_sucess.append(test_name)
		else:
			print("FAIL:", test_name)
			print("FAIL:", test_name, file=log)
			tests_fail.append(test_name)
		print("Took %s sec"%t)
		print("Took %s sec"%t, file=log)
		tot_time += t


print("\nSummary:")
print("\nSummary:", file=log)
if zgoubi_bins[0] is not None:
	print("Using zgoubi binaries:")
	print("Using zgoubi binaries:", file=log)
	for zgoubi_bin in zgoubi_bins:
		print("\t", zgoubi_bin)
		print("\t", zgoubi_bin, file=log)


print("Ran %d tests"%tests_run)
print("Ran %d tests"%tests_run, file=log)
print("Pass %d"%len(tests_sucess))
print("Pass %d"%len(tests_sucess), file=log)
print("Fail %d"%len(tests_fail))
print("Fail %d"%len(tests_fail), file=log)
if len(tests_fail) != 0:
	print("Failed tests:")
	print("Failed tests:", file=log)
	for t in tests_fail:
		print(t)

print("Took %s sec"%tot_time)
print("Took %s sec"%tot_time, file=log)

os.chdir(orig_dir)

shutil.rmtree(install_dir)
shutil.rmtree(run_dir)

print("Test written log to:", log_file_path)

if len(tests_fail) == 0:
	sys.exit(0)
else:
	sys.exit(1)

