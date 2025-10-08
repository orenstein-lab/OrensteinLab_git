# OrensteinLab

A hackable python interface for running physics lab experiments. In addition to pre-built tools for orchestrating complex measurements, it provides a bare-bones syntax and protocol to quickly develop software to coordinate instruments and motors in a modern lab environment, with a high degree of code recyclability. The basic structure is very simple:


	OrentsteinLab_git/
	  - /devices - folder for storing lowest level files and folders which handle reading and writing to various instruments and motors. Every instrument has a file in this folder following the templates shown in instrument_template.py and motor_template.py
	  - measurement.py - file containing routines used frequently during experiments. For code to be recyclable, functions here should, as much as possible, be written so as to treat motors and instruments abstractly as opposed to directly referencing files in /control. This takes some effort but is worthwhile in the long run, and underlies the power of this package.
	  - *motors_and_instruments.py - file defining the important MOTOR_DICT, INSTRUMENT_DICT, ACTIVE_MOTORS, and ACTIVE_INSTRUMENTS objects, which link abstract motors ('x', 'y', temp', etc...) to specific hardware whose interface is defined in /devices. In general this is different for each cryostat and is easy to modify, so it is not tracket by git.
	  - *configuration.py - file containing CONFIG_DICT, which holds configuration information that is specific to each system.
	  - /users - folder containing folders for each user to write their own scripts and functions, a sandbox of sorts
	  - plotters.py - file containing plotting routines for plotting in real time during experiments
	  - helper.py - file with various helper functions
	  - concurrency.py - file with multithreading routines
	  
	  * not tracked by git so that each system can have its own copy

<p align="center">
<img width="622" height="372" alt="heirarchy" src="https://github.com/user-attachments/assets/f65b67bf-e4e3-4f4d-8a1a-fbecb3d0d753" />
</p>

The philosophy behind this structure is that the same set of routines in `measurement.py` can be used throughout the lab regardless of the specific hardware used at a given cryostat. To accomplish this, high-level functions in `measurement.py` are written to handle abstracted `motors` and `instruments`, which are linked to concrete low-level hardware functions in `/devices` through the `MOTOR_DICT` and `INSTRUMENT_DICT`. Although `MOTOR_DICT` and `INSTRUMENT_DICT` define the abstracted motors and instruments available for sure, the lists `ACTIVE_INSTRUMENTS` and `ACTIVE_MOTORS` define which instruments are actually active at runtime. The scripts in `measurement.py` only can interact with activated motors and instruments. For example, most data acquisition scripts will automatically keep metadata associated with active motors and instruments. All datafiles are generate within the folder of the script that runs this pacakge.

A distinction is made between abstracted `motors` and `instruments`. `motors` are equipment that have a `move` and `read` function, for example a temperature controller. `instruments` have only a `read` function, which can return a large stream of labeled measurements. The whole point of abstracting low-level control and high-level functionality is that the user can easily and quickly define abstracted `motors` and `instruments` at will. A single piece of physical hardware equipment can be both a `motor` and and `instrument`; for example, a Zurich lock-in amplifier might be a treated as an `instrument` when it comes to reading the output of a measurement, but a `motor` can be defined to control AUX-OUT, or change the oscillator frequency, etc.

In addition, users are encouraged to come up with their own abstracted `motors` and write more complex low-level functions to achieve them. For example, in `/devices/newport.py` there are motors defined as `corotateaxes12`, which coordinate motion of two Newport axes as a single motor, and enbles treating corotation birefringence measurements as an abstract `motor`.

Each different peice of equipment thus must have its own low-level file in `/devices`, which defines all low-level actions and also wrapper functions or other kinds of primary abstractions that are then used in `MOTOR_DICT` or `INSTRUMENT_DICT`. 

Below is a quick example of what a script might look like that uses this package to take a multidimensional map and generate data. See below for explanations of some of the most useful routines in `measurement.py`

	import sys, os
	filename_head = os.getcwd()
	
	# useful import
	import numpy as np
	import pyvisa
	import matplotlib.pyplot as plt
	import xarray as xr
	
	# import high-level functionality
	from OrensteinLab_git.motors import MOTOR_DICT, INSTRUMENT_DICT
	import OrensteinLab_git.measurement as meas
	
	# import low-level functions directly into namespace for easy command line use
	from OrensteinLab_git.control.montana_temp import *
	from OrensteinLab_git.control.attocube import *
	from OrensteinLab_git.control.lakeshore import *
	from OrensteinLab_git.control.opticool import *
	from OrensteinLab_git.control.zurich import *
	from OrensteinLab_git.control.newport import *
	from OrensteinLab_git.control.razorbill import *
	
	# read current position
	read_x(); read_y()
	
	# take an xy scan, specifiying which demodulator channel to display to view DC reflectivity
	meas.motor_scan({'x':(0,100,20,{}), 'y':(0,100,20,{})}, filename='testscan', vars=['Demod 1 r'])

	# move to a position and take a time-resolved measurement as a function of temperature,showcasing how the same function be used for very different types of scans. Now it displays two demodulator values
	meas.motor_scan({'time_delay':(0,3000,10,{}), 'temp':(300,200,10,{'wait_time':30})}, filename='test_pumpprobe', vars=['Demod 1 r', 'Demod 2 r'])

# Reference
	
Althought the goal is that this architecture give users maximum flexibility when writing code, the functionality hinges on a few 'protocol' like practices when designing code. These are outlined below:

### Guidelines for writing low-level hardware control files

Files in `/devices`  handling low-level control of hardware that are to represent `motors` have at minimum the following functions (can be named whatever but must have this structure):
1. `initialize()`: function that returns a handle `obj`, which is used to communicate with the instrument. All other functions rely on `obj` to work. Sometimes differnent motors must use the same handle (such as an Attocube controller). In this case care must be needed to keep track of handles externally so as to not overwrite handles when one handle needs to be linked to several motors (for example, `obj` objects can often be pickled so that they can be accessed be different functions or threads that don't share a namespace. I am sure that more complex systems could be used to make this more robust, but in the spirit of simplicity I leave it up to the user to decide how best to handle these situations. 
2. `read(obj=None, **kwargs)`: file that whose first argument is `obj`, which is a  `kwarg` that defaults to `None`. If `obj=None`, then the function must automatically call `initialize()` to genereate a local `obj`. Returns a `float` representing the measured value. 
3. `move(position, obj=None, **kwargs)`: file that whose first argument is `position`, the position to which the motor must move to, and whose second argument is `obj`, which is a  `kwarg` that defaults to `None`. If `obj=None`, then the function must automatically call `initialize()` to genereate a local `obj`. Doesn't return anything.
4. `close(obj)`: a function which properly shutsdown a connection via handle `obj`, if necessary. Many motors don't need to do much but the motor must have this function anyway for consistency.
 
Files in `/devices`  handling low-level control of hardware that are to represent `instruments` have at minimum the following functions (can be named whatever but must have this structure):
1. `initialize()`: function that returns a handle `obj`, which is used to communicate with the instrument. All other functions rely on `obj` to work. 
2. `read(obj=None, **kwargs)`: file that whose first argument is `obj`, which is a  `kwarg` that defaults to `None`. If `obj=None`, then the function must automatically call `initialize()` to genereate a local `obj`. Returns a `dict` representing the data, whose `key:value` pairs are the name of measured quantity (for example, 'Demod 1 x') and the associate measured value. Note that this is a very different output from a `motor`, which returns a `float`; herein lies another distinction 
4. `close(obj)`: a function which properly shutsdown a connection via handle `obj`, if necessary. Many motors don't need to do much but the motor must have this function anyway for consistency.

Note that `read` and `move` functions can also have additional `**kwargs` arguments, for example wait times or ramp rates or anything like that.

System specific information (hardware addresses, for example) are stored in `CONFIG_DICT` in `configuration.py`, and so files and functions in `\devices` should always pull such information from there. If such "hard-coded" information is needed for new hardware, new entries should be entered into `CONFIG_DICT`.

### Guidelines for defining `MOTOR_DICT`, `INSTRUMENT_DICT`, `ACTIVE_MOTORS`, and `ACTIVE_INSTRUMENTS`

`MOTOR_DICT` is a dictionary where each `motor` is defined abstractly by creating an entry of the following form:

`'motor_name':{'read':read_func, 'move':move_func, 'init':init_func, 'close':close_func, 'move_back':move_back, 'name':'Motor Name'}`

that is, each motor is `key:value` pair where `key` is the abstracted name for the motor which can be called in high-level functions, and `value` is another dictionary which links the abstracted `read`,`move`,`init`,`close` functions to low-level functions in `/control`. Here `'move_back'` is an entry which is used to specify how much to overshoot when returning to initial position during a multidimensional motor scan, and `'name'` is the `string` that labels the motor when written into data files. In principle, additional entries to this dictionary can be made if new functionalities are required in the future.

`INSTRUMENT_DICT` is similar, expect an `instrument` does not have a `'move'` function. There is also a `'name'` entry, which can be used (for example, to differentiate between two different lock-in amplifiers).

`ACTIVE_MOTORS`: This is a `list` of abstract motors from the `MOTOR_DICT`, and is used to specify what motors should be actively considered within the "system" during operation. Any `motor` in `ACTIVE_MOTORS` will be read during data acquisition and added as either metadata to file headers or just directly into the datastream. Care must be taken when writing files in `measurement.py` to properly use `ACTIVE_MOTORS`.

`ACTIVE_INSTRUMENTS`: Similar to `ACTIVE_MOTOR`, but for `instruments`. All active instruments are measured in `measurement.py` scripts.
### Guidelines for writing functions in `measurement.py`

This is the trickiest part of contributing to this package, as the architecture above must be leveraged to gain any advantages. Arguable, this file is also the place that needs to most work, as it is a bit of a mess. But it is a contained mess! Anyway, here are a few things to keep in mind when writing new functions:

(1) As much as possible, code that can be recycled for different kinds of motors (for example, scanning) should be written so as to work with any kind of abstract motor. The arguments to these functions should therefore rely on names of abstract motors above all else. 

(2)Make sure to make it easy to pass `**kwargs` to motor `read` and `move` functions, for maximal functionality (this makes it easy to change all sorts of paramters in your scanning functions without actually having a million entries in the function itself). See note on `mkwargs_dict` and `ikwargs_dict` arguments below.

(3) Follow some protocols for writing to files. Every data file should have the following structure. There are functions in `helper.py` to quickly write data files in a standard format:

	[Metadata]
	...
	...
	...
	[Data]
	Header
	...
	...
	...
	
where the elipses are there data entries go and `Header` labels the data columns (ie, the `'name'` entry for each `motor` or `instrument` that will be read). 

Importantly, when doing a motor scan the following protocol is adheared to: motor values that are being scanned are written to under `'name'`, whereas values that are actually measured at each point in the scan are written as `'name' (measured)`, so that there is not ambiguity as to whether the recorded value is nominal or actually measured.

(4) Use and contribute to helper functions that make it easier to do common tasks such as retrieving information about motors, creating file headers, writing to files, etc.

(5) Make use of functions in `plotters.py` for all plotting - these routines for initializing and updating plots is very flexible.

### usage of some key functions in 	`measurement.py`

This needs to get fleshed out, but for now I just want to make a few comments:

(1) all scanning/mapping functions work based of an argument called `map_dict`, which is a a dictionary where each entry is of the form `'motor_name:(start, stop, step, kwarg_dict)'`, or alternatively `'motor_name:(positions, kwarg_dict)'`, and where entries to the 	`kwarg_dict` are `key:value` pairs where `key` is the name of a kwarg to the`move` function of the given motor, and `value` is the value to set that kwarg to. For example, you might specify a wait time with `kwarg_dict={'wait_time':30}` to have the motor wait 30 seconds after each movement.

(2) `motor_sequence` works similarly, except that the `map_dict` is replaced by a `sequence_list` defining a set of consecutive motor sweeps. Here, each entry is a tuple `(motor_name, target, kwargs_move_dict, tolerance)` and the function will acquire data as it sweeps the motors in order from left to right.

(3) `timed_measurement` is used to just record the lockin for some specified time, can be useful

(4) most functions in `measurement.py` have the following optional arguments:

		- `mkwargs_read_dict`: dictionary where keys are `motor` names and values are a dictionary of kwargs for motor read functions. if left empty, defaults to no kwargs
		- `ikwargs_dict`: dictionary of key value pairs where keys are name of `instruments` and values are dictionary of kwargs for the instrument read function. Defaults to instruments in `ACTIVE_INSTRUMENTS` with no kwargs.
		- `mobj_dict`: dictionary of key:value pairs where keys are `motor` names and values are `motor` handle object. This makes it very easy to pass around all active motor objects.
		- `iobj_dict`: dictionary of key:value pairs where keys are `instrument` names and values are `instrument` handle object. This makes it very easy to pass around all active motor objects.
		- `vars`: list of variables to plot during measurement. Defaults to `DEFAULT_VARS`.

(5) `monitor_motors` plots the position of a motor as a function of time. can also be useful

(6) `filename_head` automatically is set to directory that script is run in. If you want to save data to another directory, change `filename_head`.


### Some other tips

(1) To find pymeasure instruments run the following:

	from pymeasure.instruments import list_resources
	list_resources()

