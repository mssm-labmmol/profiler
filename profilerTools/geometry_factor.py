#
# This file is part of the profilerTools suite (see
# https://github.com/mssm-labmmol/profiler).
#
# Copyright (c) 2020 mssm-labmmol
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

from functools import partial
import numpy as np


def geom_c6(r):
    return -(r**(-6))


def geom_c12(r):
    return (r**(-12))


def geom_tors_force_constant(m, phi):
    return np.cos(m*np.radians(phi))


def geom_tors_phase(m, phi):
    return np.sin(m*np.radians(phi))


def geom_tors_ryckaert(m, phi):
    return (np.cos(np.radians(phi))) ** (m-1)


def select_geometry_function(selstring, dihtype):
    """Selects an appropriate geometry function based on the string
    `selstring` and on the dihedral type `dihtype`."""
    if (selstring == 'cs6'):
        return geom_c6
    elif (selstring == 'cs12'):
        return geom_c12
    elif (selstring.startswith('k_')):
        m = int(selstring.split('_')[1])
        if (dihtype == 'standard'):
            return partial(geom_tors_force_constant, m)
        elif (dihtype == 'ryckaert'):
            return partial(geom_tors_ryckaert, m)
    elif (selstring.startswith('phi_')):
        m = int(selstring.split('_')[1])
        return partial(geom_tors_phase, m)


def calculate_dofs(conf, atom_idxs, selstring):
    """Calculates interatomic distances or dihedral-angle values based on the
    configuration `conf`, on a matrix of atom indexes `atom_idxs` (either with 2
    or with 4 columns) and on a string `selstring` identifying the type of
    parameter.
    """
    if atom_idxs.size != 0: # do something only if list is not empty
        if (selstring == 'cs6') or (selstring == 'cs12'):
            if atom_idxs.shape[1] != 2:
                raise ValueError
            return conf.getDistances(atom_idxs[:,0], atom_idxs[:,1])
        else:
            if atom_idxs.shape[1] != 4:
                raise ValueError
            return conf.getDihedrals(atom_idxs[:,0], atom_idxs[:,1], atom_idxs[:,2],
                                    atom_idxs[:,3])
    else:
        return []
    

def calculate_geometry_factor(conf, selstring, dof_atom_idxs, dihtype):
    """Calculates geometry factors for configuration `conf` and optimized
    parameter `p`. For instance, if `p` represents a CS6 parameter, this returns
    a list of -r_{ij}**(-1/6) for all pairs (i,j) that are modelled by this
    parameter (r is the interatomic distance).
    
    :param conf: A :class:`~profilerTools.configuration.configuration`.
    
    :param selstring: A string identifying the nature of the parameter (e.g.,
                      'cs6', 'k_3', 'phi_2',...).
    
    :param dof_atom_idxs: A matrix of indexes of the atoms involved in the DOFs
                          corresponding to this parameter; this is either a
                          matrix of pairs (ai, aj), with two columns, or a
                          matrix of quadruples (ai, aj, ak, al), with four
                          columns.
    
    :param dihtype: String indicating the dihedral type ('standard' or
                    'ryckaert').
               
    :returns: A :class:`~np.ndarray` with the geometric factors for all
              optimized DOFs corresponding to parameter `p` and configuration
              `conf`.

    """
    geom_func = select_geometry_function(selstring, dihtype)
    dofs = calculate_dofs(conf, dof_atom_idxs, selstring)
    return np.array( [geom_func(dof) for dof in dofs] )

