======================
PyZgoubi Installation
======================

PyZgoubi is developed on Linux, but can also be installed on MacOS X and Windows.

Requirements
------------

PyZgoubi follows standard Python install methods, but requires a few packages to be installed first.


General Requirements
^^^^^^^^^^^^^^^^^^^^

PyZgoubi requires Python 3.7 or newer. It also requires the python math and science packages, Numpy and SciPy. For some graphical functions it uses Matplotlib.

Linux
"""""
On linux Python 3 is usually already installed. Other packages can be installed using the package manager, for example on Fedora/Red Hat/Centos::

  sudo dnf install python3-numpy python3-scipy python3-matplotlib

Or on Ubuntu/Debian::

  sudo apt-get install python3-numpy python3-scipy python3-matplotlib

MacOS X
"""""""

On MacOS X the simplest method to install the python and the scientific packages is with Miniconda from https://docs.conda.io/en/latest/miniconda.html


Windows
"""""""

On Windows the simplest method to install the python and the scientific packages is with Miniconda from https://docs.conda.io/en/latest/miniconda.html



Additional Requirements
^^^^^^^^^^^^^^^^^^^^^^^

Some additional packages are useful if you want to run development versions of PyZgoubi and Zgoubi.

* GCC/GFortran, for building Zgoubi
* Subversion, source code management for Zgoubi
* Git, source code management for PyZgoubi
* MinGW, on windows required for building Zgoubi

Linux
"""""

GCC/GFortran, Subversion and Git can be installed with the package manager,for example on Fedora/Red Hat/Scientific::

  sudo yum install gcc gcc-gfortran subversion git patch

Or on Ubuntu/Debian::

  sudo apt-get install gcc gfortran subversion git

Windows
"""""""

You need Mingw to compile Zgoubi on windows. Get the file mingw-get-setup.exe from http://sourceforge.net/projects/mingw/ or http://www.mingw.org/

From the MinGW Installation manager make sure you have::

  mingw-developer-toolkit
  mingw32-base
  mingw32-gcc-fortran

You need subversion to get development versions of Zgoubi https://subversion.apache.org/
Windows packages can be downloaded from http://www.sliksvn.com/en/download

You need Bazaar to get development versions of PyZgoubi
http://bazaar.canonical.com/


Zgoubi
^^^^^^

Zgoubi is the tracking engine at the heart of PyZgoubi. Compiled binaries and source code of Zgoubi can be downloaded from the `Zgoubi website <https://sourceforge.net/projects/zgoubi/>`_

Some features in PyZgoubi require features in the development branch of Zgoubi. The development source can be downloaded using SVN::

  svn checkout svn://svn.code.sf.net/p/zgoubi/code/trunk zgoubi-code

.. _automaticzgoubiinstall:

Automatic Zgoubi install
""""""""""""""""""""""""

PyZgoubi includes an automatic script for building a development version of Zgoubi that has been confirmed to be compatible. Once PyZgoubi and the additional requirements above have been installed Zgoubi can be built with the command::

  pyzgoubi --install-zgoubi

You will then need to edit your PyZgoubi configuration file to use this build. Edit the settings.ini file in the .pyzgoubi folder in your home directory, and add the path line given by the previous command.

Specific versions can be built::

  pyzgoubi --install-zgoubi 570
  pyzgoubi --install-zgoubi list   # to see a list of avaliable versions

Additional build options can be given::

  pyzgoubi --install-zgoubi MXSTEP=1000

This allows adjusting some of the compile time parameters found in the include/ files, for example

+----------------+-------------------------------------------+
| Parameters     |                                           |
+================+===========================================+
| MXSTEP         | Maximum steps per magnet                  |
+----------------+-------------------------------------------+
| MXL            | Maximum elements in line                  |
+----------------+-------------------------------------------+
| MXT            | Maximum number of particles               |
+----------------+-------------------------------------------+
| MXX, MXY       | Maximum steps in X and Y for field maps   |
+----------------+-------------------------------------------+
| IZ             | Maximum steps in Z for field maps         |
+----------------+-------------------------------------------+


Automatic Zgoubi install Windows
""""""""""""""""""""""""""""""""

You need to have the following already installed PyZgoubi, Subversion and minGW.

To activate the minGW tools for this terminal you need to temporarily add to the PATH variable.::

  set PATH=%PATH%;C:\MinGW\bin;C:\MinGW\msys\1.0\bin

then build with::

  pyzgoubi --install-zgoubi

if you are have 32bit MinGW you might need to do::

  pyzgoubi --install-zgoubi 365_32bit

You will then need set this zgoubi in the pyzgoubi settings. edit the settings.ini file in .pyzgoubi in your home folder, and add a line like::

  zgoubi_path = C:\Users\sam\.pyzgoubi\bin\zgoubi_365_32bit.exe


Manual Zgoubi install windows
"""""""""""""""""""""""""""""

In C:\MinGW\msys\1.0 you will find msys.bat. Double clicking this will give you a terminal with the GNU tools required for building Zgoubi.

Make sure the minGW path is set properly in msys (see http://www.mingw.org/wiki/Getting_Started ). For me I had to rename the file fstab.sample to fstab in C:\MinGW\msys\1.0\etc


Manual install of Zgoubi: (see below for auto install)
To build Zgoubi download the source zgoubi-5.1.0.zip from https://sourceforge.net/projects/zgoubi/
Extract the source folder, and cd into it from the msys terminal::

    cd /c/Users/sam/Documents/zgoubi-5.1.0/

edit the Makefile in this folder. You need to add a '#' to the lines that contain 'cd zpop', as zpop cannot be built on Windows.

Then type::

  make

If this completes without errors, try running::

  zgoubi/zgoubi.exe

You should see a message::

  PGM ZGOUBIL error open file zgoubi.dat

This means that Zgoubi has compiled.



Installation
------------

Once the requirements are installed, PyZgoubi itself can be installed.

Installing from the Python Package Index
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

PyZgoubi can be installed with the tool pip::

  pip install --user pyzgoubi


Linux and Mac OS X
^^^^^^^^^^^^^^^^^^

I recommend that the PyZgoubi source and install are kept together in a folder. Create a new folder, for example::

  mkdir ~/zgoubi
  cd ~/zgoubi
  mkdir install
  mkdir src

To install are release version download the tar.gz file from the `PyZgoubi website <http://www.hep.manchester.ac.uk/u/samt/pyzgoubi/>`_. Download it into the ~/zgoubi/src directory, and unzip it::

  cd ~/zgoubi/src
  tar -xf pyzgoubi-0.4.91.tar.gz
  cd pyzgoubi-0.4.91

Or to get a development version::

  cd ~/zgoubi/src
  git clone https://git.code.sf.net/p/pyzgoubi/code-git pyzgoubi-trunk
  cd pyzgoubi-trunk

Then to install::

  ./setup.py install --prefix=~/zgoubi/install

You will then nee to edit you bash set up use PyZgoubi. Edit you .bashrc, add the lines as instructed by the output of the previous command, e.g. ::

  export PYTHONPATH=/home/sam/zgoubi/install/lib/python2.7/site-packages:$PYTHONPATH
  export PATH=/home/sam/zgoubi/install/bin:$PATH

To check your install open a new terminal and run::

  pyzgoubi --version

If you want to use PyZgoubi's auto install script for Zgoubi check the 'Automatic Zgoubi Install' section now.

You may need to edit you settings.ini file in the .pyzgoubi folder, to adjust settings and to set the path to your Zgoubi install.

Windows
^^^^^^^

To keep everything neat its is work making a folder to keep PyZgoubi source code and install in.::

  C:\Users\sam>mkdir pyzgoubi
  C:\Users\sam>cd pyzgoubi
  C:\Users\sam\pyzgoubi>mkdir source
  C:\Users\sam\pyzgoubi>mkdir install
  C:\Users\sam\pyzgoubi>cd source

To get the current developement (you need Bazaar installed) version run::

  C:\Users\sam\pyzgoubi\source>git clone https://git.code.sf.net/p/pyzgoubi/code-git pyzgoubi-trunk
  C:\Users\sam\pyzgoubi>cd pyzgoubi-trunk
  C:\Users\sam\pyzgoubi\source\pyzgoubi-trunk>python setup.py install --prefix=C:\Users\sam\pyzgoubi\install

The installer will prompt you to add another path to your PATH, eg::

  C:\Users\sam\pyzgoubi\install\Scripts

To do this go to the environment variables control panel, and add to the end of the Path variable (using a semicolon ';' to separate it from the existing entries).

If you then open a new command prompt, and run::

  pyzgoubi --version

you should see some output showing which version of PyZgoubi and its dependencies you are running.



