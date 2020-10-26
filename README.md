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

`profilerTools` is composed of Python scripts, so no installation is
necessary. To execute a script, pass it to a Python interpreter with
the appropriate arguments, as described in the
[documentation](./doc/src/doc.pdf).

# Requirements

* `python` (preferably >3.0)
* `numpy`

# Limitations

* The program is designed for small systems, and, for this reason,
  *all* nonbonded interactions are calculated explicitly (*i.e.*, no
  cutoff truncation).
  
* The current implementation of molecular-mechanics calculations is
  very far from optimal in terms of speed (even without taking the
  previous item into account). For this reason, the duration of
  `profilerOpt` jobs increase very fastly with the number and size of
  the systems and also with the size of the genetic-algorithm
  population.
  
* With parallelization, `profilerOpt` jobs with a given RNG seed are
  *not* reproducible (without parallelization, they are).

# Utilities

To facilitate the use of `profilerTools`, we provide additional
programs to convert topologies from other formats to the .stp
format. They are

* `gromacs_itp_to_stp`, to convert from GROMACS .itp files;

* `gromos_top_to_stp`, to convert from GROMOS .top files;

These programs, each in a separate directory, can be found in the
directory `utils`. For usage instructions, consult their README files.
