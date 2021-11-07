# Description

`profilerTools` is a suite of Python programs to perform molecular torsional
scans and parameterize the related torsional and nonbonded parameters. It is
based on (classical) molecular mechanics.

`profilerTools` is composed of two scripts:

* `profilerGen`, for obtaining torsional-scan trajectories and their energy
  profiles;

* `profilerOpt`, for the parameterization of torsional terms and the related 1,4
  Lennard-Jones terms;

For implementation details, usage instructions and a tutorial, please consult
the [documentation](./doc/src/doc.pdf).

# Installation

`profilerTools` depends on [numpy](https://github.com/numpy/numpy),
[DEAP](https://github.com/DEAP/deap) and
[scikit-learn](https://github.com/scikit-learn/scikit-learn), which
should be automatically installed by `setuptools` when following the
instructions below.

## Recommended Installation

The recommended installation also requires the `numpy` C header files
and the GNU Scientific Library (GSL) header files. Therefore, you will
face issues if you do not have GSL installed or if your `numpy`
installation does not include the header files. For the APT package
manager, this can be remedied with

	sudo apt install python-numpy libgsl-dev

For other package managers (e.g., YUM and RPM), the GSL may be
available under the name `gsl-devel`. Alternatively, you can
install these libraries from their official repositories.

Once these requirements are taken care of, you can install
`profilerTools` with the command

	python setup.py build install

After that, the two executables `profilerOpt` and `profilerGen` should
be available as commands in your PATH.

Note that, depending on your default installation path, you may need
to run this command with administrative privileges.

## Suboptimal Installation

If you have any problems obtaining the `numpy` and/or GSL header
files, you can install `profilerTools` without relying on them with
the command

    python setup.py build install --without-geometrypy --without-gsl

After that, the two executables `profilerOpt` and `profilerGen` should
be available as commands in your PATH.

Note that, depending on your default installation path, you may need
to run this command with administrative privileges.

The consequences of this installation when compared with the
recommended one are:

- a worse performance of the molecular-mechanics calculations
- no support for regularized regression in the optimization

# Requirements

* `python` (>=3.6)
* [numpy](https://github.com/numpy/numpy)
* [DEAP](https://github.com/DEAP/deap)
* [scikit-learn](https://github.com/scikit-learn/scikit-learn)
* [GNU Scientific Library](https://www.gnu.org/software/gsl/)

# Limitations

* The program is designed for small systems, and, for this reason, *all*
  nonbonded interactions are calculated explicitly (*i.e.*, no cutoff
  truncation).

* The current implementation of molecular-mechanics calculations is not
  optimal in terms of speed (even without taking the previous item
  into account).

# Utilities

To facilitate the use of `profilerTools`, we provide additional programs to
convert topologies from other formats to the .stp format. They are

* `gromacs_itp_to_stp`, to convert from GROMACS .itp files;

* `gromos_top_to_stp`, to convert from GROMOS .top files;

These programs, each in a separate directory, can be found in the directory
`utils`. For usage instructions, please consult their README files.
