# Description

`gromacs_itp_to_stp` is a Python script to convert GROMACS files with the .itp extension to the .stp extension supported by `profilerTools`. It is intended to facilitate the usage of `profilerTools` to GROMACS users.

# Dependencies

Required:

* `python` (preferably >3.0)
* `numpy`

Optional:

* `cpp` for interpretation of C-like preprocessing directives (e.g., `#include`).

Also, bear in mind that that the script has only been tested in Linux systems.

# Disclaimer

While every effort has been made to cover the most common functional forms and force fields and also to provide an accurate file type conversion, the script may not apply to some systems and, like any other computer program, may have bugs. We strongly encourage the user to read the limitations below before using the script, and also to verify the correctness of the output files before using them for other purposes.

# Limitations

* The script can only process absolute paths in `#include` directives.

* The script can only process the GROMACS directives listed below. If your .itp file contains more than one molecule definition (e.g., it involves solute and solvent) or a feature that is not supported (e.g., position restraints), please modify it first.
    * `[ defaults ]`
    * `[ atomtypes ]`
    * `[ bondtypes ]`
    * `[ angletypes ]`
    * `[ dihedraltypes ]`
    * `[ moleculetype ]` (only one in the entire .itp file)
    * `[ atoms ]`
    * `[ bonds ]`
    * `[ angles ]`
    * `[ dihedrals ]`
    * `[ pairs ]`
    * `[ exclusions ]`

* The script supports only the functional forms that are also supported in `profilerTools`, i.e. those with the following type codes:
    * Bonds: 1, 2.
    * Angles: 1, 2, 5.
    * Dihedrals: 1, 2, 3, 4, 5, 9.
  
# Usage

The usage of the script can be consulted by running it with the help option `-h`, which outputs the message below.
Optional arguments are enclosed in square brackets.

```
usage: gromacs_itp_to_stp.py [-h] -f <input.itp> [-o <output.stp>] [--no-inter]

arguments:
  -f <input.itp>   input *.itp file
  -o <output.stp>  output *.stp file
  --no-inter       suppress interface to write profilerTools blocks [ refdihedral ] and [ opt* ]
```
