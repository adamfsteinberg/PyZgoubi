#!/usr/bin/env python

"""Functions to download and build zgoubi. To use run::

	pyzgoubi --install-zgoubi
	or
	pyzgoubi --install-zgoubi version

To so a list of avaliable versions run::

	pyzgoubi --install-zgoubi list


"""

from __future__ import division, print_function
import os
import sys
import shutil
import subprocess
import urllib2
from zgoubi import common

class ZgoubiBuildError(Exception):
	"Problem building Zgoubi"
	pass


zgoubi_build_dir = os.path.join(os.path.expanduser("~"), ".pyzgoubi", "build")
zgoubi_build_dir2 = os.path.join(zgoubi_build_dir, "zgoubi-trunk")
zgoubi_install_dir = os.path.join(os.path.expanduser("~"), ".pyzgoubi", "bin")
#zgoubi_svn_address = "https://zgoubi.svn.sourceforge.net/svnroot/zgoubi/trunk"
#zgoubi_svn_address = "http://svn.code.sf.net/p/zgoubi/code/trunk"
zgoubi_svn_address = "svn://svn.code.sf.net/p/zgoubi/code/trunk"

if sys.platform == "win32": exe_ext = ".exe"
else: exe_ext = ""

def check_for_programs():
	"Check that required programs are installed"
	devnull = open(os.devnull, "w")
	try:
		subprocess.call(['svn', '--version'], stdout=devnull)
	except OSError:
		raise ZgoubiBuildError("svn not found: install subversion")
	try:
		subprocess.call(['patch', '--version'], stdout=devnull)
	except OSError:
		raise ZgoubiBuildError("patch not found: install patch")
	try:
		subprocess.call(['make', '--version'], stdout=devnull)
	except OSError:
		raise ZgoubiBuildError("make not found: install make")
	try:
		subprocess.call(['gfortran', '--version'], stdout=devnull)
	except OSError:
		raise ZgoubiBuildError("gfortran not found: install gfortran")


def get_zgoubi_svn():
	"Download zgoubi from SVN"
	if os.path.isdir(zgoubi_build_dir2):
		print("Zgoubi build folder already exists:", zgoubi_build_dir)
		ret = subprocess.call(['svn', 'info'], cwd=zgoubi_build_dir2)
		if ret == 0:
			print("Zgoubi SVN already present at:", zgoubi_build_dir)
			print("Reverting local changes")
			ret2 = subprocess.call(['svn', 'revert', '-R', '.'], cwd=zgoubi_build_dir2)
			return
		
		if ret != 0 or ret2 != 0:
			raise ZgoubiBuildError("%s already exists, but does not contain a working checkout. Please remove it."%zgoubi_build_dir)

	print("Downloading Zgoubi SVN:", zgoubi_svn_address)
	common.mkdir_p(zgoubi_build_dir)
	ret = subprocess.call(['svn', 'co', zgoubi_svn_address, "zgoubi-trunk"], cwd=zgoubi_build_dir)
	if ret != 0:
		raise ZgoubiBuildError("SVN download failed")


def set_zgoubi_version(version=None):
	"Set downloaded SVN to a given version. or, if no version given, to latest version"
	ret = subprocess.call(['svn', 'revert', '-R', '.'], cwd=zgoubi_build_dir2)
	if version is None:
		ret = subprocess.call(['svn', 'switch', '^/trunk'], cwd=zgoubi_build_dir2)
		ret = subprocess.call(['svn', 'update', '--non-interactive'], cwd=zgoubi_build_dir2)
	elif "svnr" in version:
		rev =  version["svnr"]
		ret = subprocess.call(['svn', 'switch', '^/trunk'], cwd=zgoubi_build_dir2)
		ret = subprocess.call(['svn', 'update', '-r', '%s'%rev, '--non-interactive'], cwd=zgoubi_build_dir2)
	elif "svntag" in version:
		tag = version["svntag"]
		ret = subprocess.call(['svn', 'switch', '^/tags/'+tag], cwd=zgoubi_build_dir2)
		ret = subprocess.call(['svn', 'update', '--non-interactive'], cwd=zgoubi_build_dir2)
	if ret != 0:
		raise ZgoubiBuildError("SVN update failed")
	ret = subprocess.call(['svn', 'revert', '-R', '.'], cwd=zgoubi_build_dir2)
	if ret != 0:
		raise ZgoubiBuildError("SVN update failed at final revert")


def apply_zgoubi_patches(patches):
	"Download and apply a set of patches"
	for patch in patches:
		print("applying", patch)
		patchname = patch.rpartition("/")[2]
		pf = open(os.path.join(zgoubi_build_dir2, patchname), "w")
		try:
			pf.write(urllib2.urlopen(patch).read())
		except urllib2.HTTPError as e:
			raise ZgoubiBuildError("Patch download failed (Error %s): %s" % (e.code, patch))
			
		pf.close()
		ret = subprocess.call('patch -p0 < %s' % patchname, cwd=zgoubi_build_dir2, shell=True)
		if ret != 0:
			raise ZgoubiBuildError("Patch application failed: %s" % patch)

include_files = "zgoubi/PARIZ.H include/FILPLT.H include/FILSPN.H include/MAXTRA.H include/MXFS.H include/MXLD.H include/MXSCL.H include/MXSTEP.H".split()

def edit_includes(new_params):
	# read each inlcude file
	for ifile_name in include_files:
		if not os.path.exists(os.path.join(zgoubi_build_dir2, ifile_name)): continue
		org_file_content = open(os.path.join(zgoubi_build_dir2, ifile_name)).readlines()
		newfile = open(os.path.join(zgoubi_build_dir2, ifile_name), "w")
		for line in org_file_content:
			# find the parameter lines
			if not line.lower().startswith("c") and line.strip().startswith("PARAMETER"):
				params = line.strip()[9:]
				params = params.lstrip("( ").rstrip(") ")
				params = params.split(",")
				params_d = {}
				# get the old parameters
				for p in params:
					k,dummy,v = p.partition("=")
					k = k.strip()
					v = v.strip()
					params_d[k] = v
				comments = []
				# update with new_params
				for k,v in params_d.items():
					if k in new_params.keys():
						params_d[k] = new_params[k]
						comments.append("changed %s from %s to %s"%(k,v,new_params[k]))
						print("in", ifile_name ,comments[-1])
				# if a change was made write it to the file
				if comments:
					params_s = ",".join(["%s=%s"%(k,v) for k,v in params_d.items()])
					newline = "      PARAMETER (%s)\n" % params_s
					for comment in comments:
						newfile.write("C   pyzgoubi build script:" + comment+"\n")
					newfile.write("C"+line)
					newfile.write(newline)
				else:
					# otherwise write the original line
					newfile.write(line)
			else:
				# non active parameter lines are passed through
				newfile.write(line)

def make_zgoubi(makecommands, makecleancommands=None, threads=2):
	"Build zgoubi source code"
	print("building zgoubi")
	if makecleancommands:
		for makecleancommand in makecleancommands:
			command = makecleancommand.split()
			ret = subprocess.call(command, cwd=zgoubi_build_dir2)
			if ret != 0:
				raise ZgoubiBuildError("Make clean failed")
	for makecommand in makecommands:
		command = makecommand.split()+ ['-j%d'%threads]
		ret = subprocess.call(command, cwd=zgoubi_build_dir2)
		if ret != 0:
			raise ZgoubiBuildError("Building zgoubi failed:" + " ".join(command))

def install_zgoubi(suffix="", install_zpop=True):
	"Install zgoubi into ~/.pyzgoubi folder"
	common.mkdir_p(zgoubi_install_dir)

	shutil.copy(os.path.join(zgoubi_build_dir2, "zgoubi", "zgoubi"+exe_ext),
	                os.path.join(zgoubi_install_dir, "zgoubi%s"%suffix+exe_ext))
	if install_zpop:
		shutil.copy(os.path.join(zgoubi_build_dir2, "zpop", "zpop"),
		            os.path.join(zgoubi_install_dir, "zpop%s"%suffix))


zgoubi_versions = {}
zgoubi_versions["261+patches"] = dict(svnr=261,
patches=[
"http://www.hep.man.ac.uk/u/sam/pyzgoubi/zgoubipatches/mcmodel-fix.diff",
"http://www.hep.man.ac.uk/u/sam/pyzgoubi/zgoubipatches/zgoubi_parallel_build.diff",
"http://www.hep.man.ac.uk/u/sam/pyzgoubi/zgoubipatches/objet3.diff",
"http://www.hep.man.ac.uk/u/sam/pyzgoubi/zgoubipatches/kobj301.diff",
],
)

zgoubi_versions["329+patches"] = dict(svnr=329,
patches=[
"http://www.hep.man.ac.uk/u/sam/pyzgoubi/zgoubipatches/build_tweaks2.diff",
"http://www.hep.man.ac.uk/u/sam/pyzgoubi/zgoubipatches/kobj301_3.diff",
],
makecommands=["make -f Makefile_zgoubi_gfortran"],
makecommands_zpop=["make -f Makefile_zpop_gfortran"],
includes={"MXSTEP":10000},
)

zgoubi_versions["360+patches"] = dict(svnr=360,
patches=[
"http://www.hep.man.ac.uk/u/sam/pyzgoubi/zgoubipatches/build_tweaks2.diff",
"http://www.hep.man.ac.uk/u/sam/pyzgoubi/zgoubipatches/kobj301_4.diff",
],
makecommands=["make -f Makefile_zgoubi_gfortran"],
makecommands_zpop=["make -f Makefile_zpop_gfortran"],
includes={"MXSTEP":10000},
)

zgoubi_versions["365"] = dict(svnr=365,
patches=[
"http://www.hep.man.ac.uk/u/sam/pyzgoubi/zgoubipatches/build_tweaks2.diff",
#"http://www.hep.man.ac.uk/u/sam/pyzgoubi/zgoubipatches/kobj301_4.diff",
],
makecommands=["make -f Makefile_zgoubi_gfortran"],
makecommands_zpop=["make -f Makefile_zpop_gfortran"],
includes={"MXSTEP":10000},
)

#on 32bit some arrays need shrinking
zgoubi_versions["365_32bit"] = dict(svnr=365,
patches=[
"http://www.hep.man.ac.uk/u/sam/pyzgoubi/zgoubipatches/build_tweaks2.diff",
#"http://www.hep.man.ac.uk/u/sam/pyzgoubi/zgoubipatches/kobj301_4.diff",
],
makecommands=["make -f Makefile_zgoubi_gfortran"],
makecommands_zpop=["make -f Makefile_zpop_gfortran"],
includes={"MXSTEP":10000, "MXX":801, "MXY":29},
)
if sys.platform == "win32":
	zgoubi_versions["365_32bit"]["patches"][0] = "http://www.hep.man.ac.uk/u/sam/pyzgoubi/zgoubipatches/build_tweaks2_windows.diff"
	zgoubi_versions["365"]["patches"][0] = "http://www.hep.man.ac.uk/u/sam/pyzgoubi/zgoubipatches/build_tweaks2_windows.diff"

zgoubi_versions["427"] = dict(svnr=427,
patches=[
"http://www.hep.man.ac.uk/u/samt/pyzgoubi/zgoubipatches/build_tweaks3.diff",
"http://www.hep.man.ac.uk/u/samt/pyzgoubi/zgoubipatches/r430_impmod.diff",
"http://www.hep.man.ac.uk/u/samt/pyzgoubi/zgoubipatches/r430_rndm.diff",
],
makecommands=["make -f Makefile_zgoubi_gfortran"],
makecommands_zpop=["make -f Makefile_zpop_gfortran"],
includes={"MXSTEP":10000},
)

zgoubi_versions["431"] = dict(svnr=431,
patches=[
"http://www.hep.man.ac.uk/u/samt/pyzgoubi/zgoubipatches/build_tweaks_r431.diff",
"http://www.hep.man.ac.uk/u/samt/pyzgoubi/zgoubipatches/makefile_funcd.diff",
"http://www.hep.man.ac.uk/u/samt/pyzgoubi/zgoubipatches/r430_impmod.diff",
"http://www.hep.man.ac.uk/u/samt/pyzgoubi/zgoubipatches/r430_rndm.diff",
],
makecommands=["make -f Makefile_zgoubi_gfortran"],
makecommands_zpop=["make -f Makefile_zpop_gfortran"],
includes={"MXSTEP":10000},
)

zgoubi_versions["437"] = dict(svnr=437,
patches=[
"http://www.hep.man.ac.uk/u/samt/pyzgoubi/zgoubipatches/build_tweaks_r437.diff",
#"http://www.hep.man.ac.uk/u/samt/pyzgoubi/zgoubipatches/makefile_funcd.diff",
#"http://www.hep.man.ac.uk/u/samt/pyzgoubi/zgoubipatches/r430_impmod.diff",
#"http://www.hep.man.ac.uk/u/samt/pyzgoubi/zgoubipatches/r430_rndm.diff",
],
makecommands=["make -f Makefile_zgoubi_gfortran"],
makecommands_zpop=["make -f Makefile_zpop_gfortran"],
includes={"MXSTEP":10000},
)

zgoubi_versions["437_nonative"] = dict(svnr=437,
patches=[
"http://www.hep.man.ac.uk/u/samt/pyzgoubi/zgoubipatches/build_tweaks_r437_nonative.diff",
#"http://www.hep.man.ac.uk/u/samt/pyzgoubi/zgoubipatches/makefile_funcd.diff",
#"http://www.hep.man.ac.uk/u/samt/pyzgoubi/zgoubipatches/r430_impmod.diff",
#"http://www.hep.man.ac.uk/u/samt/pyzgoubi/zgoubipatches/r430_rndm.diff",
],
makecommands=["make -f Makefile_zgoubi_gfortran"],
makecommands_zpop=["make -f Makefile_zpop_gfortran"],
includes={"MXSTEP":10000},
)

zgoubi_versions["535"] = dict(svnr=535,
patches=[
"http://www.hep.man.ac.uk/u/samt/pyzgoubi/zgoubipatches/build_tweaks_r532.diff",
],
makecommands=["make -f Makefile_zgoubi_gfortran"],
makecleancommands=["make -f Makefile_zgoubi_gfortran clean"],
makecommands_zpop=["make -f Makefile_zpop_gfortran"],
includes={"MXSTEP":10000},
)

zgoubi_versions["570"] = dict(svnr=570,
patches=[
"http://www.hep.man.ac.uk/u/samt/pyzgoubi/zgoubipatches/build_tweaks_r558.diff",
"http://www.hep.man.ac.uk/u/samt/pyzgoubi/zgoubipatches/zpop_pltopt_write.diff",
"http://www.hep.man.ac.uk/u/samt/pyzgoubi/zgoubipatches/zgoubi_impfai_inquire.diff",
],
makecommands=["make -f Makefile_zgoubi_gfortran"],
makecleancommands=["make -f Makefile_zgoubi_gfortran clean"],
makecommands_zpop=["make -f Makefile_zpop_gfortran"],
includes={"MXSTEP":10000},
)

zgoubi_versions["967"] = dict(svnr=967,
patches=[
"http://www.hep.man.ac.uk/u/samt/pyzgoubi/zgoubipatches/build_tweaks_r962.diff",
],
makecommands=["make"],
makecleancommands=["make clean"],
makecommands_zpop=["make MakeFiles/Makefile_zpop_gfortran"],
includes={"MXSTEP":10000},
)

zgoubi_versions["1451"] = dict(svnr=1451,
patches=[
"http://www.hep.man.ac.uk/u/samt/pyzgoubi/zgoubipatches/build_tweaks_r1451.diff",
],
makecommands=["make"],
makecleancommands=["make clean"],
makecommands_zpop=["make MakeFiles/Makefile_zpop_gfortran"],
includes={"MXSTEP":10000},
)

zgoubi_versions["6.0.1"] = dict(svntag="zgoubi-6.0.1",
patches=[
"http://www.hep.man.ac.uk/u/samt/pyzgoubi/zgoubipatches/build_tweaks_r558.diff",
#"http://www.hep.man.ac.uk/u/samt/pyzgoubi/zgoubipatches/zpop_pltopt_write.diff",
#"http://www.hep.man.ac.uk/u/samt/pyzgoubi/zgoubipatches/zgoubi_impfai_inquire.diff",
],
makecommands=["make -f Makefile_zgoubi_gfortran"],
makecleancommands=["make -f Makefile_zgoubi_gfortran clean"],
makecommands_zpop=["make -f Makefile_zpop_gfortran"],
includes={"MXSTEP":10000},
)

zgoubi_versions["6.0.2"] = dict(svntag="zgoubi-6.0.2",
patches=[
"http://www.hep.man.ac.uk/u/samt/pyzgoubi/zgoubipatches/build_tweaks_r558.diff",
#"http://www.hep.man.ac.uk/u/samt/pyzgoubi/zgoubipatches/zpop_pltopt_write.diff",
#"http://www.hep.man.ac.uk/u/samt/pyzgoubi/zgoubipatches/zgoubi_impfai_inquire.diff",
],
makecommands=["make -f Makefile_zgoubi_gfortran"],
makecleancommands=["make -f Makefile_zgoubi_gfortran clean"],
makecommands_zpop=["make -f Makefile_zpop_gfortran"],
includes={"MXSTEP":10000},
)

def install_zgoubi_all(version="6.0.2", include_opts=None):
	"This currently install a version of zgoubi known to work with pyzgoubi"
	if include_opts is None: include_opts = {}
	check_for_programs()
	if sys.platform.startswith('linux'):
		build_zpop = True
	else:
		build_zpop = False

	if version in ['list', 'help']:
		print("Available versions:", " ".join(zgoubi_versions.keys()))
		exit(0)
	if version not in zgoubi_versions:
		raise ZgoubiBuildError("Unknown version: "+ version+ "\nTry "+ " ".join(zgoubi_versions.keys()))
	print("Preparing to install zgoubi:", version)
	get_zgoubi_svn()
	set_zgoubi_version(zgoubi_versions[version])
	apply_zgoubi_patches(zgoubi_versions[version]['patches'])

	# edit any of the include files
	include_opts_all = {}
	if "includes" in zgoubi_versions[version]:
		include_opts_all.update(zgoubi_versions[version]["includes"])
	include_opts_all.update(include_opts)
	if include_opts_all:
		edit_includes(include_opts_all)
	
	make_zgoubi(zgoubi_versions[version].get("makecommands", ['make']), zgoubi_versions[version].get("makecleancommands", ["make clean"]) )
	if build_zpop:
		makecommands_zpop = zgoubi_versions[version].get("makecommands_zpop", None)
		if makecommands_zpop:
			make_zgoubi(makecommands_zpop)


	install_zgoubi("_"+version, build_zpop)
	print("\nInstalled zgoubi into ", zgoubi_install_dir)
	print("Add the following line to ~/.pyzgoubi/settings.ini")
	print("zgoubi_path =", os.path.join(zgoubi_install_dir, "zgoubi_"+version+exe_ext))

if __name__ == "__main__":
	install_zgoubi_all()
