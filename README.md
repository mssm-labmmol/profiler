# Description

`profilerTools` is a suite of Python programs to perform molecular
torsional scans and parameterize the related torsional and nonbonded
parameters. It is based on (classical) molecular mechanics.

`profilerTools` is composed of two scripts:

* `profilerGen`, for obtaining torsional-scan trajectories and their
  energy profiles;

* `profilerOpt`, for the parameterization of torsional terms and the
  related 1,4 Lennard-Jones terms;

For implementation details, usage instructions and a tutorial, please
consult the [documentation](./doc/src/doc.pdf).

# Installation

`profilerTools` depends on `numpy`, `DEAP` and `scikit-learn`, which
should be automatically installed by following the instructions
below.

To install, run

	python setup.py install

After that, the two executables `profilerOpt` and `profilerGen` should
be available as commands in your PATH.

Note that, depending on your default installation path, you may need
to run this command with administrative privileges.

You may face issues when installing if your `numpy` installation did
not include the C header files (e.g. when using `pip`). In this case,
you may face compilation errors like

	fatal errror: numpy/arrayobject.h: No such file or directory

For users of the APT package manager, this can usually be remedied with

	sudo apt install python-numpy

# Requirements

* `python` (preferably >3.4)
* `numpy`
* `DEAP`
* `scikit-learn`

# Limitations

* The program is designed for small systems, and, for this reason,
  *all* nonbonded interactions are calculated explicitly (*i.e.*, no
  cutoff truncation).
  
* The current implementation of molecular-mechanics calculations is
  very far from optimal in terms of speed (even without taking the
  previous item into account). For this reason, the duration of
  `profilerOpt` jobs increase very fastly with the number and size of
  the systems and also with the size of the evolutionary-algorithm
  population.
  
# Utilities

To facilitate the use of `profilerTools`, we provide additional
programs to convert topologies from other formats to the .stp
format. They are

* `gromacs_itp_to_stp`, to convert from GROMACS .itp files;

* `gromos_top_to_stp`, to convert from GROMOS .top files;

These programs, each in a separate directory, can be found in the directory
`utils`. For usage instructions, please consult their README files.
