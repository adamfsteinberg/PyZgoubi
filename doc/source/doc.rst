.. include:: include.txt

======================
PyZgoubi Documentation
======================

The latest version of PyZgoubi and this document can be found at http://pyzgoubi.org/

Introduction
------------

`Zgoubi <http://sourceforge.net/projects/zgoubi>`_ is a particle tracking code maintained by François Méot. PyZgoubi is an interface to Zgoubi written in python. It aims to make input files that are easy to read and can contain calculations, loops, and any other python program feature.

A basic knowledge of `Python <http://python.org/>`_ is useful to use PyZgoubi. Python is a general purpose, high level, interpreted programming language. I recommend the `Official Python Tutorial <http://docs.python.org/tutorial/>`_.

Quick Start
-----------

Put the following in a text file named quickstart.py::

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
    add(FAISCNL(FNAME='zgoubi.fai')) # record point
    add(END())

    print(output()) # file fed to zgoubi (zgoubi.dat)
    run() # run the line in zgoubi
    print(res()) # show the zgoubi.res file
    print(get_all('fai')) # so the data recorded with FAISCNL

The file can also be found in examples/quickstart.py

Run it with::

    pyzgoubi quickstart.py

You should see the zgoubi data file, the output from zgoubi, the res file from zgoubi, and the coordinates at the 2 FAISCNL points displayed to the terminal.

For more advance use it may be better store the |Line| with a name, this allows you to work with multiple |Lines| simultaneously. The following shows the same file, but with a named |Line|::

    myline = Line('line')

    # create and OBJET2 and add 1 electron to it
    ob = OBJET2()
    ob.set(BORO=ke_to_rigidity(10e6, ELECTRON_MASS))
    ob.add(Y=0, T=0.1, D=1)
    myline.add(ob)

    myline.add(FAISCNL(FNAME='zgoubi.fai')) # record point
    d1 = DRIFT(XL=50)
    q1 = QUADRUPO(XL=20, B_0=5, R_0=10, XPAS=(10,20,10))
    d2 = DRIFT(XL=50)
    q2 = QUADRUPO(XL=20, B_0=-5, R_0=10, XPAS=(10,20,10))
    myline.add(d1, q1, d2, q2)
    myline.add(FAISCNL(FNAME='zgoubi.fai')) # record point
    myline.add(END())

    print(myline.output()) # file fed to zgoubi (zgoubi.dat)
    myresults = myline.run() # run the line in zgoubi
    print(myresults.res()) # show the zgoubi.res file
    print(myresults.get_all('fai')) # so the data recorded with FAISCNL

Generating zgoubi.dat files
---------------------------

The simplest level of Pyzgoubi is to use it to generate a zgoubi.dat file. You can then run Zgoubi with the file to get your results.

PyZgoubi works with |Lines|, which are made up of Elements. The Elements are mostly the same as those used in Zgoubi. Most Elements correspond to magnets and other beam line elements of a physical accelerator, e.g. DIPOLE, MUTIPOL and CAVITE. Some are used to define the beam e.g. OBJET2. Some are to control the execution of Zgoubi, e.g. FAISCNL, MATRIX and END. A number of these Elements can be added to a |Line|. This |Line| then has all the information needed to make a zgoubi.dat file. 

A |Line| can be created as follows::

    emma = Line('emma simulation')

emma is the name used to refer to the |Line| in the program. 'emma simulation' is the text that will be put at the start of the zgoubi.dat file.

To create an Element use::

    q1 = QUADRUPO('defoc', XL=20, R_0=2, B_0=2, XPAS=0.1)

This can now be added to the |Line| ::

    emma.add(q1)

These last 2 steps can be contracted::

    emma.add(QUADRUPO('defoc', XL=20, R_0=2, B_0=2, XPAS=0.1))

however this means that there is no reference to the element that could be used for modifying its parameters later, e.g.::

    q1.set(B_0=3)

The |Line| can now be used to output a zgoubi.dat using the output() function::

    print emma.output()

All these instructions can be put in a text file and run using the command pyzgoubi. (pyzgoubi is an alias set up by the installer). Below is a section of emma.py from the examples::

        emma = Line('emma')
        xpas = (20,20,20)

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

        fb = -6.695 * fr * T
        db = 4.704 * dr * T

        ob = OBJET2()
        emma.add(ob)

        emma.add(ELECTRON())

        emma.add(DRIFT('ld', XL=ld*cm_/2))
        emma.add(CHANGREF(ALE=angle))

        emma.add(CHANGREF(YCE=d_offset*cm_))
        emma.add(QUADRUPO('defoc', XL=dq*cm_, R_0=dr*cm_, B_0=db*kgauss_, XPAS=xpas))

        emma.add(CHANGREF(YCE=-d_offset*cm_))

        emma.add(DRIFT('sd', XL=sd*cm_))

        emma.add(CHANGREF(YCE=f_offset*cm_))
        emma.add(QUADRUPO('foc', XL=fq*cm_, R_0=fr*cm_, B_0=fb*kgauss_, XPAS=xpas))
        emma.add(CHANGREF(YCE=-f_offset*cm_))

        emma.add(DRIFT('ld', XL=ld*cm_/2))

        emma.add(FAISCNL(FNAME='zgoubi.fai'))

        emma.add(REBELOTE(K=99, NPASS=10))

        emma.add(END())

        rigidity = ke_to_rigidity(10e6, 0.51099892e6)
        ob.set(BORO=-rigidity)
        ob.add(Y=0, T=0, D=1)

        print(emma.output())

This can be run with the command::

    pyzgoubi emma.py


It will build a |Line| and print the zgoubi.dat input to the screen. The '.py' extension is not necessary, but will cause your text editor to use python syntax highlighting.

Particles / Particuls
"""""""""""""""""""""

Zgoubi allows a particles parameters such as charge, mass and half-life to be set with the PARTICUL element. These not needed for basic tracking, as the rigidity is given in the OBJET element. PyZgoubi offers some pre-made particles::

	ELECTRON
	PROTON
	MUON
	IMMORTAL_MUON

The IMMORTAL_MUON is a muon with an infinite lifetime, for use when decay is not needed. Anti particles with opposite can me made by negating a particle. eg::

	m = -MUON()
	my_line.add(m)
	#or
	my_line.add(-MUON())

The masses and charges are defined in zgoubi/constants.py


Advanced Elements
-----------------

Some Zgoubi element can represent a set of magnets, this can be useful defining doublets and triplets.

These elements are used slightly differently to simpler elements. The non looping section is defined normally.::


    triplet = FFAG('triplet', IL=0, AT=10 ... )

Then the looped part can be added::

    triplet.add(ACN = 6, BZ_0 = 0.5 ...)
    triplet.add(ACN = 4, BZ_0 = -0.5 ...)
    triplet.add(ACN = 6, BZ_0 = 0.5 ...)

N gets automatically set. All the looped parts can be removed using::

    triplet.clear()


There is also a FAKE_ELEM element. This allows you to put arbitrary text into the zgoubi.dat file. It is useful for using an Zgoubi element that pyzgoubi does not have a definition for. For example::

    change_txt = """'CHANGREF'
    5.0 0 10.0
    """
    change = FAKE_ELEM(change_txt)
    line.add(change)


Available Elements
""""""""""""""""""

- BEND
- CAVITE
- CHANGREF
- DRIFT
- ELECTRON
- END
- FAISCEAU
- FAISCNL
- FAISTORE
- FAKE_ELEM
- FFAG
- IMMORTAL_MUON
- MARKER
- MATRIX
- MULTIPOL
- MUON
- OBJET1
- OBJET2
- OBJET5
- PARTICUL
- PROTON
- QUADRUPO
- REBELOTE
- TOSCA

To find the full list of elements available in the current version run::

    pyzgoubi --help elements

To find the names of the parameters available for an element use::

    pyzgoubi --help element_name

e.g.::

    pyzgoubi --help MULTIPOL

Use this in combination with the Zgoubi manual. Most parameters have the same name. Greek letters in the manual are usually anglicised  (e.g. τ to tau). Subscripts are represented with underscores. When there is a entrance and exit version of a parameter E (entrée) and S (sortie) are used to distinguish.

Running Zgoubi
--------------

Once a |Line| has been created and had the needed elements added it can be run. PyZgoubi will take care of creating a temporary directory, creating the zgoubi.dat file and running Zgoubi. This is done to prevent zgoubi from overwriting any existing files. If you wish to keep any of the output files you must use the commands to copy these to where you want them.


The following example shows how to run a |Line|::


    #create line
    emma = Line('emma')

    #add elements
    emma.add( ...  )
    ...

    #run line
    emma_res = emma.run()

    #save output
    emma_res.save_res("emma.res")
    emma_res.save_plt("emma.plt")

Note that you will need to make sure your |Line| will actually create plt or fai files, otherwise you will receive a file not found error. See the Zgoubi manual for more information.

The |run| function can take several options. If you want to inspect the directory where zgoubi is run, or to use zpop, then set xterm=True. If you want to change the directory that zgoubi is run in you can use the tmp_prefix option. It is best to make sure this is a local disk (i.e. not a network/remote disk). The default directory can be set in the zgoubi_settings.py file.::

    emma.run(xterm=True)
    emma.run(tmp_prefix = '/var/tmp/sam/')
    emma.run(xterm=True, tmp_prefix = '/home/sam/tmp/')

If you want to do analysis of the simulation you can use the |Results| object that is returned by the |run| function.::

    res = emma.run()

See the |Results| Object chapter for more info.

Each time a |Line| is run a temporary director is created. These are normally automatically cleared up when PyZgoubi finishes (also the /tmp directory is usually emptied when a computer shuts down). However if you are making repeated calls to |run|, then you may want to manually clear away these files. This can be done with the |clean| function. Don't clean the |Results| until you have finished working with its output files.::

    res = emma.run()
    res.clean()


Results Object
--------------

When you run a |Line| it creates a |Results| object, that can be used to get information about the paths of the particles.::

    res = pamela.run()

|get_all| and |get_track|, let you get lists of the particle coordinates. They each need to be told if they should read the plt (points within the magnetic elements) or fai (beam at FAISCNL element). |get_all| returns a list of dictionaries, containing all the coordinates and information. |get_track| returns a list of lists of just the requested coordinates.::

    print res.get_all('plt')
    print res.get_all('fai')
    print res.get_track('fai', ['Y','T'])

Assuming that you are using a recent Zgoubi version get_all will return a numpy structured array (recarray) which has named columns. To see the names of the columns use::

	fai_data = res.get_all('fai')
	print fai_data.dtype.names

To find the units look inside a zgoubi.fai file.

Bunch Objects
-------------

New in Pyzgoubi 0.4

The |Bunch| object represents a bunch of particles. Each particle has coordinates D, Y, T, Z, P, S, tof, and X, and the whole bunch has the shared properties rigidity/kinetic energy, particle mass, particle charge. These are stored in m, rad, eV, s, and automatically converted to Zgoubi units by the associated  functions.

It can be used in 2 ways. Firstly with the OBJET_bunch element, which behaves similarly to the real OBJET elements. Secondly it can be used to drive Zgoubi in a more abstracted method.

To use with OBJET_bunch, first create a bunch, then create an OBJET_bunch, then add it to a |Line|::

	my_bunch = Bunch(nparticles=5, ke=1e6, mass=PROTON_MASS, charge=1)
	my_bunch.particles()[0]['Y'] = 3
	...

 	ob = OBJET_bunch(bunch=my_bunch)
	my_line = Line("a line")
	my_line.add(ob)
	my_line.add( ... )

This |Line| can then be run as before.

The second method is to create a |Line| that has no OBJET, PARTICUL or END elements, only beamline elements. Then the |Line.track_bunch| method can be called, which will return the tracked bunch::

	my_line = Line("a line")
	my_line.add( QUADRUPO( ... ) )
	my_line.add( DRIFT( ... ) )
	...
	my_bunch2 = my_line.track_bunch(my_bunch)

This allows tracking a bunch around multiple lattices. Suppose that you create 3 |Lines|: injecttion_line, ring and extraction_line. You can then take a bunch through each in turn::

	my_bunch = Bunch( ... )

	my_bunch = injecttion_line.track_bunch(my_bunch)
	for n in xrange(n_laps):
		my_bunch = ring.track_bunch(my_bunch)
	my_bunch = extraction_line.track_bunch(my_bunch)

Note that this method does not allow you to access the Result object.

If you have a multi-CPU or multi-core CPU, then you can swap |Line.track_bunch| for the multithreaded version |Line.track_bunch_mt|. The multithreaded version also has the advantage that it can track an arbitrarily large bunch (more than Zgoubi's max particles limit).

There are a number of generators for standard bunches, e.g Kapchinskij-Vladimirskij (KV), waterbag and Guassian::

    gen_gauss_x_xp_y_yp(npart, emit_y, emit_z, beta_y, beta_z, alpha_y, alpha_z, seed=None, ke=None, rigidity=0, mass=0, charge=1)
    gen_gauss_x_xp_y_yp_s_dp(npart, emit_y, emit_z, beta_y, beta_z, alpha_y, alpha_z, mom_spread=0, bunch_length=0, disp=0, disp_prime=0, seed=None, ke=None, rigidity=0, mass=0, charge=1)
    gen_halo_x_xp_y_yp(npart, emit_y, emit_z, beta_y, beta_z, alpha_y, alpha_z, seed=None, ke=None, rigidity=0, mass=0, charge=1)
    gen_kv_x_xp_y_yp(npart, emit_y, emit_z, beta_y, beta_z, alpha_y, alpha_z, seed=None, ke=None, rigidity=0, mass=0, charge=1)
    gen_waterbag_x_xp_y_yp(npart, emit_y, emit_z, beta_y, beta_z, alpha_y, alpha_z, seed=None, ke=None, rigidity=0, mass=0, charge=1)

For example to create a KV bunch with 1000 protons::

    my_bunch = Bunch.gen_kv_x_xp_y_yp(1000, 1e-3, 1e-3, 4, 5, 1e-3, 2e-2, ke=50e6, mass=PROTON_MASS, charge=1)



Complex Lines
-------------

It is often useful to break down a complex lattice into sections and cells. This can be done in several ways in PyZgoubi. 

Firstly you can use python features such as loops, e.g. to put 5 identical FODO cells you could use ::

    line = Line("example")

    line.add( ... )
    line.add( ... )
    
    for x in range(5):
        line.add(QUADRUPO( ... ))
        line.add(DRIFT( ... ))
        line.add(QUADRUPO( ... ))
        line.add(DRIFT( ... ))

    line.add( ... )

Note that the for loop block lasts as long as the code is indented.

If you want to make one iteration different then you can do a test based on the x ::

    for x in range(5):
        line.add(QUADRUPO( ... ))
        line.add(DRIFT( ... ))
        line.add(QUADRUPO( ... ))
        line.add(DRIFT( ... ))
        if (x == 2): #note x counts from zero
            line.add(DIPOLE( ... )

Another method is to make a |Line| for each cell, and then build the full lattice from the cells::
	
	cell1 = Line("cell1")
	cell1.add(...)
	cell1.add(...)

	cell2 = Line("cell2")
	cell2.add(...)
	cell2.add(...)

	full_line = Line("fullline")
	full_line.add(cell1)
	full_line.add(cell2)
	#or
	full_line.add(cell1, cell2)
	
The +, * and - operators are overloaded for the line to add, repeat and reverse them.For example, assuming that you have defined the Lines ``injection``, ``extraction`` and ``ring_cell``::

	full_line = injection + n * ring_cell + extraction

To use a cell with its elements in reverse order::

	full_line.add(-cell1)


Command line arguments
----------------------
There is a convenience function for using command line arguments. For example near the start of the code put::

    number_of_laps = int(get_cmd_param('laps', 10))

Then when creating the |Line| use::

    l.add(REBELOTE(NPASS=number_of_laps-1, KWRIT=1, K=99))

to access the variable. When you run the your simulation use the command line argument::

    pyzgoubi sim.py laps=50

I the second parameter to get_cmd_param() is the default. If you don't give a value as an argument then the default is used. If you don't give a default eg::

    particle_energy = float(get_cmd_param('energy'))

then you will receive and error if you don't provide a value as an argument.


Units
-----

PyZgoubi does not do any automatic unit conversion. When you give a parameter you must give the units that Zgoubi expects. However PyZgoubi does define some values to save some effort, in the conversion. Any of the following will set x to 2::

	x = 2 * m
	x = 200 * cm
	x = 2000 * mm

Then to output x in a different unit::

	print x * m_ , "m"
	print x * cm_ , "cm"
	print x * mm_ , "mm"

So for example DRIFT expects cm, but your lattice might use meters, so do::

	d_length = 0.6 * m
	DRIFT('d1', XL=d_length*cm_)

Currently the following units are defined::

	m
	cm
	mm
	T
	kgauss

but more can be added on request.

For conversion between degrees and radians use the python math functions::

	degrees(2*pi) # gives 360
	radians(180) # gives 3.1416...

Settings and Options
--------------------

PyZgoubi settings are stored in a file .pyzgoubi/settings.ini in your home folder. It is automatically created the first time pyzgoubi is used. It can be used to customise some PyZgoubi options.

The following keys can be set::

	extra_defs_files

Give the path of any additional definition files you want to be considered::

	tmp_dir

Where temporary files should be written, this is most likely /tmp/, but in some case you may wish to use /var/tmp/ or a ramdisk /dev/shm/::

	zgoubi_path 

The path to the zgoubi binary file. Note that this can also be set with the commandline option --zgoubi=/path/to/zgoubi.

Debugging and Profiling
"""""""""""""""""""""""
PyZgoubi can be run with pythons interactive mode (same as "python -i") so that in the event of an error the user is given a python prompt to inspect variables at the point of the exception.::

    pyzgoubi -i script.py

PyZgoubi can generate profiles using the cProfile module. use::

    pyzgoubi --profile script.py

The profiling information can be read with the python pstats module (part of cProfile). For example to see in which functions most time was spent run::

	python -c "import pstats; pstats.Stats('prof.log').sort_stats('cumulative').print_stats()"

or start a python shell and run::

	>>> import pstats
	>>> p = pstats.Stats('prof.log')
	>>> p.sort_stats('cumulative').print_stats()


.. _Logging:

Logging levels
""""""""""""""

The verbosity of PyZgoubi can be adjusted. By default the log_level is set to 'warn', so only warnings and error messages are printed. One can raise the level to 'debug' which will also show debug messages. To do this for a single run use::

	pyzgoubi --debug script
	or
	pyzgoubi --log_level=debug script

or you can adjust the value in the settings.ini file::

	log_level = debug

Other levels can also be set, see http://docs.python.org/library/logging.html for more information.


Upgrade Notes
-------------

Although it would be nice to have perfect backwards compatibility that sometimes interferes with progress, and things have to break. There should be no breaks between a version X.Y.Z and X.Y.Z+1, but there can be between X.Y.Z and X.Y+1.0.

0.8.x
"""""

PyZgoubi 0.8 uses Python 3 and no longer support Python 2.

If you are porting old pyzgoubi scripts you may find the following document useful https://docs.python.org/3/howto/pyporting.html

In many cases you may be able to do the conversion automatically with the 2to3 tool.

0.8 also uses yaml definitions of elements, so it is not currently possible for the user to add local element definitions. Instead new elements should be added to the https://github.com/PyZgoubi/zgoubi-metadata project.

0.3.x to 0.4.x
""""""""""""""

PyZgoubi 0.4 supports the new fai/plt output formats introduced into Zgoubi in early 2010. These have a header that labels the columns. Reading the new format was taken as an opportunity to use numpy more extensively. If an older version of Zgoubi (5.1 or 5.0) is being used then the old fai/plt reading code will be used, and data will be returned as python dictionaries and arrays. If a newer version of Zgoubi is being used (SVN r251 or newer), then Results.get_all() will return a `numpy structured array <http://docs.scipy.org/doc/numpy/user/basics.rec.html>`_ with the column names.

Some column names will change when using the new fai/plt files. This is because older versions of PyZgoubi muddled the S and X coordinates. Also 'D' and 'D0' are now the more accurate 'D-1' and 'D0-1'. The coordinate name is now taken from the fai/plt file directly.

When using the new fai/plt files labels are stored as a fixed length string, and so include any whitespace, e.g. 'foc' vs. 'foc'. To  get back to the short version use strip(), e.g.::

	label1 = 'foc     '
	label_short = label1.strip()
	#or
	labels1 = ['foc     ', 'defoc   ']
	labels2 = [x.strip() for x in labels1]

The 'tol' parameter in find_closed_orbit() now is a measure of convergence rather than area. This should give better results over a wide range of scales. However, a large 'tol' is needed to give the same degree of accuracy. If you previously had tol=1e-10, then it may no longer converge, but if you change it to tol=1e-6 you will get a similar result to before.

Some obsolete functions have been removed: Results.plt_to_csv(), Results.get_all_old(), Line.split_line()

Some obsolete functions have been marked as deprecated, and will give warnings when used. These will likely be removed in a future version of PyZgoubi, unless you let me know that they are worth keeping. Some of these functions are old and undocumented, and I suspect no one uses them. Some have just move, like the functions for directly accessing and saving res and fai files (moved from Line to Results). If you do not care about modifying you're code to prepare for a future version then pass the argument ``-W ignore::DeprecationWarning`` to python.

Tips
----

Python hints
""""""""""""
If there is a '#' character on a line, everything after it is treated as a comment.

Python uses whitespace to delimit blocks (instead of braces '{' and '}' in C/C++). The PyZgoubi code uses tabs, so it is best to use tabs in your input files. If you get an 'IndentationError' check that you have not mixed spaces and tabs, or accidentally started a line with a space/tab.

Python identifiers (variable, function, object names etc) are case sensitive, they must start with a letter and only contain letters, numbers and underscores.


Zgoubi hints
""""""""""""

Errors
~~~~~~

sfe: formatted io not allowed::

	Zgoubi, version 5.0.0.
	Job  started  on  05-Feb-09,  at  15:00:29 
	Copying  zgoubi.dat  into  zgoubi.res,
	numbering  and  labeling  elements...

	acc                              Zgoubi Version 5.0.0.  

	366/  367 REBELOTE                    sfe: formatted io not allowed
	apparent state: unit 21 named zgoubi.res
	lately writing sequential formatted external IO


This may mean that you have tried to write output to both ascii and binary files, eg zgoubi.fai and b_zgoubi.fai


Troubleshooting
---------------

There are several levels at which problems can occur. Errors in the input file, bugs in PyZgoubi, bugs in Zgoubi.

Debug mode
""""""""""

Some useful information is shown when debugging is enabled, for example warnings about common mistakes. See :ref:`Logging`.


Check the line
""""""""""""""

Printing the line instance will show what elements it contains::
	
	my_line = Line("a line")
	my_line.add( QUADRUPO( ... ) )
	my_line.add( DRIFT( ... ) )
	...
	print(my_line)


Element output
""""""""""""""

Check what an element is outputting to the zgoubi.dat file::

	q1 = QUADRUPO( ... )
	print q1.output()

or the whole line::
	
	print(my_line.output)

If they are not what you expect have a look in PyZgoubi's definition files.

Res file
""""""""

Read the zgoubi.res file. It shows how Zgoubi interpreted the zgoubi.dat file::

	results = my_line.run()
	print(results.res())

Check all the units, Zgoubi uses a range of different units.

xterm
"""""

Open an xterm in the temp working folder, and have a look at the files zgoubi has output::

	results = my_line.run(xterm=True)

From here you can modify the zgoubi.dat and run zgoubi again.


Memory Issues
"""""""""""""

Zgoubi uses statically allocated arrays with sizes set at compile time for storage of elements, particle, field maps etc. When it starts it will allocate several GB of memory that it may not actually use. If you are on a low memory machine (less that 4GB) you may see the message::

    Check that you have sufficient RAM or enable memory overcommit

If you are on Linux you maybe able to work around this by adjusting your memory overcommit settings. This allows programs to make allocations larger than the available memory, they will be stopped if they actually try to use all of this memory. Unless you are using large field maps you are unlikely to use all the memory that Zgoubi allocates. To allow over committing adjust the sysctl value vm.overcommit_memory with the command::

    sysctl vm.overcommit_memory=1

Or by editing the file /etc/sysctl.conf. See https://www.kernel.org/doc/Documentation/vm/overcommit-accounting for more details.

Matplotlib Issues
"""""""""""""""""

If you see the error::

    ImportError: No module named backend_tkagg

You need to install the matplotlib tk backend, for example the package python-matplotlib-tk. If this is not possible you can switch to a different default plotting backend by editing (or creating) a matplotlib config file (.config/matplotlib/matplotlibrc or .matplotlib/matplotlibrc on linux), for example add the line::

    backend : Agg

There are a range of backends to choose from, see http://matplotlib.org/faq/usage_faq.html#what-is-a-backend .


Limit issues
""""""""""""

Zgoubi stores many things in fixed sized arrays. If you reach the limit of these you might see an error such as::

    Y-dim of map is too large, max is  200   Occured in element # 4,  at pass # 1
    
    PROCEDURE STOPPED: too many elements

If so you can rebuild Zgoubi with larger settings, see :ref:`automaticzgoubiinstall`. 








