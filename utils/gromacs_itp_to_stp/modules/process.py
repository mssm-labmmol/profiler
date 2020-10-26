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

import re
import sys
import copy

from numpy import sqrt
from .gen_pairs import gen_pairs
from .interface import write_optblocks_via_interface

class BlockLine(Exception): pass
class EmptyLine(Exception): pass
class EOFLine(Exception):   pass

def advance(stream):
    line = stream.readline()
    if (line == ''):
        raise EOFLine
    
    line = line.strip()
    if (line == ''):
        raise EmptyLine
    line = line.split(';')[0]
    flds = line.split()
    if (flds == []):
        # comment line
        return advance(stream)
    if (flds[0] == '['):
        raise BlockLine
    return flds

def dof_has_parameters (dof):
    return 'pars' in dof

def doftype_get_parameters(doftype):
    """Returns a dictionary of parameters.
    Works for bonds, angles, dihedrals and pairs."""
    doftype_copy = copy.deepcopy(doftype)
    remove_idxs = ['ai', 'aj', 'ak', 'al', 'func']
    for idx in remove_idxs:
        if idx in doftype_copy.keys():
            doftype_copy.pop(idx)
    return doftype_copy

def compare_bond_atomtypes (types_bond1, types_bond2):
    check1 = types_bond1 == types_bond2
    check2 = types_bond1 == types_bond2[::-1]
    return check1 or check2

def compare_angle_atomtypes (types_angle1, types_angle2):
    if types_angle1[1] == types_angle2[1]:
        l1 = [types_angle1[0], types_angle1[2]]
        l2 = [types_angle2[0], types_angle2[2]]
        return compare_bond_atomtypes(l1, l2)
    else:
        return False

def compare_dihedral_atomtypes (types_dih1, types_dih2):
    check_list_1 = [d1 == d2 for d1,d2 in zip(types_dih1, types_dih2)]
    # reverse order
    check_list_2 = [d1 == d2 for d1,d2 in zip(types_dih1, types_dih2[::-1])]
    return ((False not in check_list_1)) or ((False not in check_list_2))

# note: dih2 is the one that will contain the 'X' wildcard
def compare_generic_dihedral_atomtypes (types_dih1, types_dih2_in_ff):
    check_list_1 = [(d1 == d2) or (d2 == 'X') for d1,d2 in zip(types_dih1, types_dih2_in_ff)]
    # reverse order
    check_list_2 = [(d1 == d2) or (d2 == 'X') for d1,d2 in zip(types_dih1, types_dih2_in_ff[::-1])]
    return ((False not in check_list_1)) or ((False not in check_list_2))

def dof_compare(type_1, type_2, generic=False):
    nbodies = len(type_1)
    if nbodies == 2:
        return compare_bond_atomtypes(type_1, type_2)
    if nbodies == 3:
        return compare_angle_atomtypes(type_1, type_2)
    if nbodies == 4:
        if generic:
            return compare_generic_dihedral_atomtypes(type_1, type_2)
        else:
            return compare_dihedral_atomtypes(type_1, type_2)

# goto next block and return block name or -1 if EOF is reached
def goto_next_block (stream):
    while True:
        line = stream.readline()
        if line == '':
            break
        m = re.match("\[ (\w+) \]", line)
        if m:
            return m.group(1)
    return -1

def check_supported_defaults (nbfunc):
    if (nbfunc != 1):
        raise Exception("error: only LJ interactions are supported so far")
    return True

def check_supported_bondtypes (bondtype):
    if (bondtype == 1) or (bondtype == 2):
        return True
    raise Exception("error: bondtype %d not supported" % bondtype)

def check_supported_angletypes (angletype):
    if angletype in [1,2,5]:
        return True
    raise Exception("error: angletype %d not supported" % angletype)

def check_supported_dihedraltypes (dihedraltype):
    if (dihedraltype) in [1, 2, 3, 4, 5, 9]:
        return True
    raise Exception("error: dihedraltype %d not supported" % dihedraltype)

def get_angletype_shape(angletype):
    """Identifies an angletype as 'harmonic' or 'urey-bradley'"""
    check_supported_angletypes(angletype)
    if angletype == 1 or angletype == 2:
        return 'harmonic'
    elif angletype == 5:
        return 'urey-bradley'

def get_number_of_multiple_dihedrals(dihedral):
    """Get number of parameter sets for a particular dihedral."""
    return len(list(dihedral['pars'].values())[0])

def get_dihedraltype_shape (dihedraltype):
    """Identifies a dihedraltype as 'periodic' (requires multiplicity), 'harmonic', 'ryckaert' or 'fourier'."""
    check_supported_dihedraltypes(dihedraltype)
    if (dihedraltype) in [1, 4, 9]:
        return 'periodic'
    elif (dihedraltype) == 2:
        return 'harmonic'
    elif (dihedraltype) == 3:
        return 'ryckaert'
    elif (dihedraltype) == 5:
        return 'fourier'

def get_dihedraltype_type(dihedraltype):
    """Identifies a dihedraltype as proper or improper."""
    check_supported_dihedraltypes(dihedraltype)
    if dihedraltype in [1, 3, 5, 9]:
        return 'proper'
    elif dihedraltype in [2,4]:
        return 'improper'

def parse_dihedraltype_definition(dihedraltype, number_of_fields):
    """Given a dihedraltype code and a number of fields, returns a list
    with keys for each field in order."""
    if dihedraltype in [1,9]:
        if number_of_fields == 8:
            return ['ai', 'aj', 'ak', 'al', 'func', 'phi0', 'k', 'm']
        elif number_of_fields == 6:
            return ['aj', 'ak', 'func', 'phi0', 'k', 'm']
        else:
            raise Exception("Can't parse dihedraltype, incompatible number of fields.")
    elif dihedraltype == 4:
        if number_of_fields == 8:
            return ['ai', 'aj', 'ak', 'al', 'func', 'phi0', 'k', 'm']
        elif number_of_fields == 6:
            return ['ai', 'al', 'func', 'phi0', 'k', 'm']
        else:
            raise Exception("Can't parse dihedraltype, incompatible number of fields.")
    elif dihedraltype == 2:
        if number_of_fields == 7:
            return ['ai', 'aj', 'ak', 'al', 'func', 'phi0', 'k']
        elif number_of_fields == 5:
            return ['ai', 'al', 'func', 'phi0', 'k']
        else:
            raise Exception("Can't parse dihedraltype, incompatible number of fields.")
    elif dihedraltype == 3:
        if number_of_fields == 11:
            return ['ai', 'aj', 'ak', 'al', 'func', 'c0', 'c1', 'c2', 'c3', 'c4', 'c5']
        elif number_of_fields == 9:
            return ['aj', 'ak', 'func', 'c0', 'c1', 'c2', 'c3', 'c4', 'c5']
    elif dihedraltype == 5:
        if number_of_fields == 9:
            return ['ai', 'aj', 'ak', 'al', 'func', 'f1', 'f2', 'f3', 'f4']
        elif number_of_fields == 7:
            return ['aj', 'ak', 'func', 'f1', 'f2', 'f3', 'f4']
    else:
        raise Exception("Can't parse dihedraltype, unexpected type code %d." % dihedraltype)
        
def read_defaults_from_stream (stream):
    defaults = {}
    while True:
        try:
            flds = advance(stream)
            defaults['nbfunc'] = int(flds[0])
            defaults['comb-rule'] = int(flds[1])
            defaults['gen-pairs'] = flds[2]
            defaults['fudgeLJ'] = float(flds[3])
            defaults['fudgeQQ'] = float(flds[4])
        except EOFLine:
            raise Exception("Premature EOF when reading .itp file.")
        except BlockLine:
            break
        except EmptyLine:
            break
    check_supported_defaults (defaults['nbfunc'])
    return defaults
        
def read_exclusions_from_stream (stream):
    exclusions = []
    while True:
        try:
            flds = advance(stream)
            for x in range(1, len(flds)):
                exclusions.append( [flds[0], flds[x]] )
        except EOFLine:
            raise Exception("Premature EOF when reading .itp file.")
        except BlockLine:
            break
        except EmptyLine:
            break
    return exclusions

def read_atomtypes_from_stream (stream):
    atomtypes = []
    while True:
        try:
            flds = advance(stream)
            this_atomtype = {}
            if (len(flds) == 7):
                # no bond_type column
                this_atomtype['name'] = flds[0]
                this_atomtype['btype'] = flds[0]
                this_atomtype['at-num'] = int(flds[1])
                this_atomtype['mass'] = float(flds[2])
                this_atomtype['charge'] = float(flds[3])
                this_atomtype['ptype'] = flds[4]
                this_atomtype['v'] = float(flds[5])
                this_atomtype['w'] = float(flds[6])
            elif (len(flds) == 8):
                # this is specific to OPLS 
                this_atomtype['name'] = flds[0] 
                this_atomtype['btype'] = flds[1] 
                this_atomtype['at-num'] = int(flds[2])
                this_atomtype['mass'] = float(flds[3])
                this_atomtype['charge'] = float(flds[4])
                this_atomtype['ptype'] = flds[5]
                this_atomtype['v'] = float(flds[6])
                this_atomtype['w'] = float(flds[7])
            else:
                raise Exception("strange format in [ atomtypes ]")
            atomtypes.append( this_atomtype )
        except EOFLine:
            raise Exception("Premature EOF when reading .itp file.")
        except BlockLine:
            break
        except EmptyLine:
            break
    return atomtypes

def read_bondtypes_from_stream (stream):
    bondtypes = []
    while True:
        try:
            flds = advance(stream)
            this_bondtype = {}
            this_bondtype['ai'] = flds[0]
            this_bondtype['aj'] = flds[1]
            this_bondtype['func'] = int(flds[2])
            this_bondtype['b0'] = float(flds[3])
            this_bondtype['kb'] = float(flds[4])
            check_supported_bondtypes ( this_bondtype['func'] )
            bondtypes.append( this_bondtype )
        except EOFLine:
            raise Exception("Premature EOF when reading .itp file.")
        except BlockLine:
            break
        except EmptyLine:
            break
    return bondtypes

def read_angletypes_from_stream (stream):
    angletypes = []
    while True:
        try:
            flds = advance(stream)
            this_angletype = {}
            this_angletype['ai'] = flds[0]
            this_angletype['aj'] = flds[1]
            this_angletype['ak'] = flds[2]
            this_angletype['func'] = int(flds[3])
            this_angletype['th0'] = float(flds[4])
            this_angletype['k'] = float(flds[5])
            try:
                this_angletype['r13'] = float(flds[6])
                this_angletype['kub'] = float(flds[7])
            except IndexError:
                this_angletype['r13'] = 0.00
                this_angletype['kub'] = 0.00
            check_supported_angletypes ( this_angletype['func'] )
            angletypes.append( this_angletype )
        except EOFLine:
            raise Exception("Premature EOF when reading .itp file.")
        except BlockLine:
            break
        except EmptyLine:
            break
    return angletypes

# reads both proper and improper dihedrals
def read_dihedraltypes_from_stream (stream):
    dihtypes = []
    while True:
        try:
            flds = advance(stream)
            this_dihtype = {}
            this_dihtype['ai'] = 'X'
            this_dihtype['aj'] = 'X'
            this_dihtype['ak'] = 'X' 
            this_dihtype['al'] = 'X'
            if str(flds[4]).isnumeric():
                flds_keys = parse_dihedraltype_definition(int(flds[4]), len(flds))
            else:
                flds_keys = parse_dihedraltype_definition(int(flds[2]), len(flds))
            for i, key in enumerate(flds_keys):
                this_dihtype[key] = flds[i]
            this_dihtype['func'] = int(this_dihtype['func'])
            check_supported_dihedraltypes(this_dihtype['func'])
            shape = get_dihedraltype_shape(this_dihtype['func'])
            if (shape == 'periodic'):
                this_dihtype['phi0'] = float(this_dihtype['phi0'])
                this_dihtype['k'] = float(this_dihtype['k'])
                this_dihtype['m'] = int(this_dihtype['m'])
            elif (shape == 'harmonic'):
                this_dihtype['phi0'] = float(this_dihtype['phi0'])
                this_dihtype['k'] = float(this_dihtype['k'])
            elif (shape == 'ryckaert'):
                this_dihtype['c0'] = float(this_dihtype['c0'])
                this_dihtype['c1'] = float(this_dihtype['c1'])
                this_dihtype['c2'] = float(this_dihtype['c2'])
                this_dihtype['c3'] = float(this_dihtype['c3'])
                this_dihtype['c4'] = float(this_dihtype['c4'])
                this_dihtype['c5'] = float(this_dihtype['c5'])
            elif (shape == 'fourier'):
                this_dihtype['f1'] = float(this_dihtype['f1'])
                this_dihtype['f2'] = float(this_dihtype['f2'])
                this_dihtype['f3'] = float(this_dihtype['f3'])
                this_dihtype['f4'] = float(this_dihtype['f4'])
            dihtypes.append( this_dihtype )
        except EOFLine:
            raise Exception("Premature EOF when reading .itp file.")
        except BlockLine:
            break
        except EmptyLine:
            break
    return dihtypes

def read_nonbond_params_from_stream (stream):
    nonbond_params = []
    while True:
        try:
            flds = advance(stream)
            this_nonbond_param = {}
            this_nonbond_param['ai'] = flds[0]
            this_nonbond_param['aj'] = flds[1]
            this_nonbond_param['func'] = int(flds[2])
            check_supported_defaults ( this_nonbond_param['func'] )
            this_nonbond_param['v'] = float(flds[3])
            this_nonbond_param['w'] = float(flds[4])
            nonbond_params.append( this_nonbond_param )
        except EOFLine:
            raise Exception("Premature EOF when reading .itp file.")
        except BlockLine:
            break
        except EmptyLine:
            break
    return nonbond_params

def read_pairtypes_from_stream (stream):
    return read_nonbond_params_from_stream (stream)

def read_partypes_from_stream (stream):
    atomtypes = []
    bondtypes = []
    angletypes = []
    dihedraltypes = []
    nonbond_params = []
    pairtypes = []
    partypes = {'atomtypes': atomtypes,
                'bondtypes': bondtypes,
                'angletypes': angletypes,
                'dihedraltypes': dihedraltypes,
                'nonbond_params': nonbond_params,
                'pairtypes': pairtypes}

    while True:
        block = goto_next_block (stream)
        if (block == "moleculetype"):
            # this means partypes are done
            break
        if (block == -1):
            break
        if (block == "atomtypes"):
            atomtypes += read_atomtypes_from_stream (stream)
        if (block == "bondtypes"):
            bondtypes += read_bondtypes_from_stream (stream)
        if (block == "angletypes"):
            angletypes += read_angletypes_from_stream (stream)
        if (block == "dihedraltypes"):
            dihedraltypes += read_dihedraltypes_from_stream (stream)
        if (block == "nonbond_params"):
            nonbond_params += read_nonbond_params_from_stream (stream)
        if (block == "pairtypes"):
            pairtypes += read_pairtypes_from_stream (stream)
    return partypes

def read_defaults_partypes_from_stream (stream):
    while True:
        block = goto_next_block (stream)
        if (block == -1):
            break
        if (block == "defaults"):
            defaults = read_defaults_from_stream (stream)
            break
    partypes = read_partypes_from_stream (stream)
    return [defaults, partypes]

def read_atoms_from_stream (stream, atomtypes):
    atoms = []
    while True:
        try:
            flds = advance(stream)
            this_atom = {}
            this_atom['type'] = flds[1]
            this_atom['resnr'] = int(flds[2])
            this_atom['resname'] = flds[3]
            this_atom['name'] = flds[4]
            this_atom['cgnr'] = int(flds[5])
            this_atom['q'] = float(flds[6])
            this_atom['m'] = float(flds[7])
            # also include 'btype' 
            for at in atomtypes:
                if at['name'] == this_atom['type']:
                    this_atom['btype'] = at['btype']
            atoms.append( this_atom )
        except EOFLine:
            break
        except BlockLine:
            break
        except EmptyLine:
            break    
    return atoms

def read_pairs_from_stream (stream):
    pairs = []
    while True:
        try:
            flds = advance(stream)
            this_pair = {}
            this_pair['ai'] = int(flds[0])
            this_pair['aj'] = int(flds[1])
            this_pair['func'] = int(flds[2])
            if (this_pair['func'] != 1):
                raise Exception("error: in [ pairs ] block, only func = 1 is accepted.")
            if (len(flds) > 3):
                this_pair['pars'] = [float(flds[3]), float(flds[4])]
            pairs.append( this_pair )
        except EOFLine:
            break
        except BlockLine:
            break
        except EmptyLine:
            break    

    return pairs
     
def read_bonds_from_stream (stream):
    bonds = []
    while True:
        try:
            flds = advance(stream)
            this_bond = {}
            this_bond['ai'] = int(flds[0])
            this_bond['aj'] = int(flds[1])
            this_bond['func'] = int(flds[2])
            if (len(flds) >= 5):
                this_bond['pars'] = { 'kb': float(flds[4]), 'b0': float(flds[3]) }
            bonds.append( this_bond )
        except EOFLine:
            break
        except BlockLine:
            break
        except EmptyLine:
            break
    return bonds

def read_angles_from_stream (stream):
    angles = []
    while True:
        try:
            flds = advance(stream)
            this_angle = {}
            this_angle['ai'] = int(flds[0])
            this_angle['aj'] = int(flds[1])
            this_angle['ak'] = int(flds[2])
            this_angle['func'] = int(flds[3])
            if len(flds) > 4:
                try:
                    this_angle['pars'] = {'k': float(flds[5]), 'th0': float(flds[4]), 'r13': float(flds[6]), 'kub': float(flds[7])}
                except IndexError:
                    this_angle['pars'] = {'k': float(flds[5]), 'th0': float(flds[4]), 'r13': 0.00, 'kub': 0.00}
            angles.append( this_angle )
        except EOFLine:
            break
        except BlockLine:
            break
        except EmptyLine:
            break
    return angles

def read_dihedrals_from_stream (stream):
    dihedrals = []
    while True:
        try:
            flds = advance(stream)
            this_dihedral = {}
            this_dihedral['ai'] = int(flds[0])
            this_dihedral['aj'] = int(flds[1])
            this_dihedral['ak'] = int(flds[2])
            this_dihedral['al'] = int(flds[3])
            this_dihedral['func'] = int(flds[4])
            check_supported_dihedraltypes(this_dihedral['func'])
            if len(flds) > 5:
                shape = get_dihedraltype_shape(this_dihedral['func'])
                if (shape == 'periodic'):
                    this_dihedral['pars'] = { 'phi0': float(flds[5]), 'k': float(flds[6]), 'm': float(flds[7]) }
                elif (shape == 'harmonic'):
                    this_dihedral['pars'] = { 'phi0': float(flds[5]), 'k': float(flds[6]) }
                elif (shape == 'ryckaert'):
                    this_dihedral['pars'] = {
                        'c0': float(flds[5]),
                        'c1': float(flds[6]),
                        'c2': float(flds[7]),
                        'c3': float(flds[8]),
                        'c4': float(flds[9]),
                        'c5': float(flds[10])}
                elif (shape == 'fourier'):
                    this_dihedral['pars'] = {
                        'f1': float(flds[5]),
                        'f2': float(flds[6]),
                        'f3': float(flds[7]),
                        'f4': float(flds[8])}
            dihedrals.append( this_dihedral )            
        except EOFLine:
            break
        except BlockLine:
            break
        except EmptyLine:
            break
    return dihedrals


def reshape_dihedrals_for_uniqueness (dihedrals):
    """This function will take any dihedral ai-aj-ak-al that appears more
    than once in the list of dihedrals and fuse each of them into a
    single dihedral member.  Example: if dihedral 1-2-3-4 appears more
    than once in the input file, such as:
    
    ...
    1 2 3 4   9    phi1    k1    m1
    ....
    1 2 3 4   9    phi2    k2    m2
    
    it will be reshaped to a single dihedral 1-2-3-4 with
             'pars' = { 'phi0': [phi1,phi2], 'k': [k1,k2], 'm': [m1,m2] }

    """
    dihedrals_record = [] # keep track of unique dihedrals
    new_dihedrals = [] # output combining all dihedrals and their parameters
    spec_comment = "specified explicitly or in macro"
    for d in dihedrals:
        idxs = [d[x] for x in ['ai', 'aj', 'ak', 'al']]
        if not (idxs in dihedrals_record):
            dihedrals_record.append(idxs)
            # copy dihedral
            new_dihedral = copy.deepcopy(d)
            shape = get_dihedraltype_shape(d['func'])
            new_dihedral['comment'] = []
            if dof_has_parameters(d):
                # if explicit parameters, copy them as list
                for partype in d['pars'].keys():
                    new_dihedral['pars'][partype] = [d['pars'][partype]]
                new_dihedral['comment'] = [spec_comment]
            else:
                # if not explicit parameters, create empty lists to be filled later
                new_dihedral['pars'] = {}
                if (shape == 'periodic'):
                    new_dihedral['pars']['phi0'] = []
                    new_dihedral['pars']['k'] = []
                    new_dihedral['pars']['m'] = []
                elif (shape == 'harmonic'):
                    new_dihedral['pars']['phi0'] = []
                    new_dihedral['pars']['k'] = []
                elif (shape == 'ryckaert'):
                    new_dihedral['pars']['c0'] = []
                    new_dihedral['pars']['c1'] = []
                    new_dihedral['pars']['c2'] = []
                    new_dihedral['pars']['c3'] = []
                    new_dihedral['pars']['c4'] = []
                    new_dihedral['pars']['c5'] = []
                elif (shape == 'fourier'):
                    new_dihedral['pars']['f1'] = []
                    new_dihedral['pars']['f2'] = []
                    new_dihedral['pars']['f3'] = []
                    new_dihedral['pars']['f4'] = []
            new_dihedrals.append(new_dihedral)
        else:
            # raise error in case func is not 9
            if (d['func'] != 9):
                raise Exception("error: multiple dihedrals must be of type 9")
            where_is_it = dihedrals_record.index(idxs)
            # stack parameters of this dihedral
            if dof_has_parameters(d):
                new_dihedral['comment'].append(spec_comment)
                for partype in d['pars'].keys():
                    new_dihedral['pars'][partype].append(d['pars'][partype])
            else:
                raise Exception("Strange occurrence of multiple dihedrals in .itp file.")
    return new_dihedrals

def read_itp_from_stream (stream, atomtypes):
    atoms = []
    bonds = []
    angles = []
    dihedrals = []
    pairs = []
    exclusions = []

    itp = { 'atoms': atoms,
            'bonds': bonds,
            'angles': angles,
            'dihedrals': dihedrals,
            'exclusions': exclusions,
            'pairs': pairs}

    # read [ moleculetype ]
    while True:
        try:
            flds = advance(stream)
            itp['nexcl'] = int(flds[1])
        except EOFLine:
            raise Exception("Premature EOF when reading .itp file.")
        except BlockLine:
            break
        except EmptyLine:
            break

    # read other blocks
    while True:
        block = goto_next_block (stream)
        if (block == -1):
            break
        if (block == "atoms"):
            atoms += read_atoms_from_stream (stream, atomtypes)
        if (block == "bonds"):
            bonds += read_bonds_from_stream (stream)
        if (block == "angles"):
            angles += read_angles_from_stream (stream)
        if (block == "dihedrals"):
            dihedrals += read_dihedrals_from_stream (stream)
        if (block == "pairs"):
            pairs += read_pairs_from_stream (stream)
        if (block == "exclusions"):
            exclusions += read_exclusions_from_stream (stream)

    itp['dihedrals'] = reshape_dihedrals_for_uniqueness (dihedrals)
    return itp

def assign_bondtype_for_bond (bond, atoms, bondtypes):
    if (dof_has_parameters(bond)):
        bond['comment'] = "specified explicitly or in macro"
        return
    idx1 = bond['ai'] - 1
    idx2 = bond['aj'] - 1
    type1 = atoms[idx1]['btype']
    type2 = atoms[idx2]['btype']
    for bondtype in bondtypes:
        if dof_compare([type1,type2], [bondtype['ai'], bondtype['aj']]) and bond['func'] == bond['func']:
            bond['comment'] = " %s - %s " % (bondtype['ai'], bondtype['aj'])
            bond['pars'] = doftype_get_parameters(bondtype)
            return
    # oops!
    raise Exception("error: no bondtype found for bond %s - %s" % (type1, type2))

def assign_angletype_for_angle (angle, atoms, angletypes):
    if dof_has_parameters(angle):
        angle['comment'] = "specified explicitly or in macro"
        return
    idx1 = angle['ai'] - 1
    idx2 = angle['aj'] - 1
    idx3 = angle['ak'] - 1
    type1 = atoms[idx1]['btype']
    type2 = atoms[idx2]['btype']
    type3 = atoms[idx3]['btype']
    l1 = [type1, type2, type3]
    for angletype in angletypes:
        l2 = [angletype['ai'], angletype['aj'], angletype['ak']]
        if dof_compare(l1,l2) and angle['func'] == angletype['func']:
            angle['comment'] = " %s - %s - %s " % (l2[0],l2[1],l2[2])
            angle['pars'] = doftype_get_parameters(angletype)
            return
    # oops!
    raise Exception("error: no angletype found for angle %s - %s - %s" % (type1, type2, type3))

# IMPORTANT:
# ----------
#
# Dihedrals work in a very specific manner.
# Let's say you read a dihedral 1-2-3-4 under the [ dihedrals ]
# directive of your *.itp file. Regarding the assignment
# of parameters to this dihedral, there are the following
# possibilities:
#
#   1 - If type is not 9, you proceed like any other bonded type. If
#       parameters are specified in the *.itp line, assign
#       them. Otherwise, read the UNIQUE parameters from the
#       force-field bonded file.
#
#   2 - If type is 9 (multiple dihedrals), you have to check all
#       possible assignments. This means that one line in the *.itp
#       file might be converted into more than one line in the *.stp
#       file. In this case:
#
#         a - Regardless of parameters being specified in the *.itp
#             line, find non-generic (i.e. not involving X)
#             dihedraltypes that match the atom types of the dihedral.
#             Assign all compatible non-generic dihedrals.
#
#         b - If no non-generic dihedrals are found, you are allowed
#             to search among the generic dihedrals involving 'X'
#             wildcards. Assign all compatible generic dihedrals.
#
#         c - Regardless of the dihedraltypes assigned (generic or
#             non-generic), also assign the specified parameters (in
#             the same line of the *.itp file) if there are any.
#
#         d - Watch out! If the dihedral has more than one set of
#             specified parameters, it will appear more than once in
#             the original *.itp file. When this happens, you have to
#             recognize to avoid stacking the same dihedrals over and
#             over again.
#
# Function returns True if assignment was possible.
# This is important to evaluate whether generic dihedraltypes will be attempted.
def assign_nongeneric_dihedraltype_for_dihedral (dihedral, atoms, dihedraltypes):
    idx1 = dihedral['ai'] - 1
    idx2 = dihedral['aj'] - 1
    idx3 = dihedral['ak'] - 1
    idx4 = dihedral['al'] - 1
    type1 = atoms[idx1]['btype']
    type2 = atoms[idx2]['btype']
    type3 = atoms[idx3]['btype']
    type4 = atoms[idx4]['btype']
    l1 = [type1, type2, type3, type4]
    return_value = False
    for dihedraltype in dihedraltypes:
        l2 = [dihedraltype['ai'], dihedraltype['aj'], dihedraltype['ak'], dihedraltype['al']]
        if dof_compare(l1,l2) and (dihedral['func'] == dihedraltype['func']):
            dihedral['comment'].append(" %s - %s - %s - %s" % (l2[0],l2[1],l2[2],l2[3]))
            # dictionary of pars
            this_pars = doftype_get_parameters(dihedraltype)
            # now append to dihedral parameters (initialized when reshaped, so I can just append)
            for k in this_pars.keys():
                dihedral['pars'][k].append(this_pars[k])
            return_value = True
    # is true if any assingment was made
    return return_value

def assign_generic_dihedraltype_for_dihedral (dihedral, atoms, dihedraltypes):
    idx1 = dihedral['ai'] - 1
    idx2 = dihedral['aj'] - 1
    idx3 = dihedral['ak'] - 1
    idx4 = dihedral['al'] - 1
    type1 = atoms[idx1]['btype']
    type2 = atoms[idx2]['btype']
    type3 = atoms[idx3]['btype']
    type4 = atoms[idx4]['btype']
    l1 = [type1, type2, type3, type4]
    for dihedraltype in dihedraltypes:
        l2 = [dihedraltype['ai'], dihedraltype['aj'], dihedraltype['ak'], dihedraltype['al']]
        if dof_compare (l1,l2, generic=True) and (dihedral['func'] == dihedraltype['func']):
            dihedral['comment'].append(" %s - %s - %s - %s" % (l2[0],l2[1],l2[2],l2[3]))
            # dictionary of pars
            this_pars = doftype_get_parameters(dihedraltype)
            # now append to dihedral parameters (initialized when reshaped, so I can just append)
            for k in this_pars.keys():
                dihedral['pars'][k].append(this_pars[k])
    # oops!
    # if no parameters were assigned...
    if (list(dihedral['pars'].values())[0] == []):
        raise Exception("error: no dihedraltype found for dihedral %s - %s - %s - %s" % (type1, type2, type3, type4))

def is_pair_in_pairs (pair, pairs):
    for p in pairs:
        check_1 = [pair['ai'],pair['aj']] == [p['ai'],p['aj']]
        check_2 = [pair['aj'],pair['ai']] == [p['ai'],p['aj']]
        if (check_1) or (check_2):
            return True
    return False

def assign_pairtype_for_pair (pair, atoms, pairs, partypes, comb_rule, gen_pairs, fudgeLJ, fudgeQQ):

    q_i = atoms[pair['ai']-1]['q']
    q_j = atoms[pair['aj']-1]['q']
    q_ij = q_i * q_j

    # prepare default LJ values
    type_list = [at['name'] for at in partypes['atomtypes']]
    v_list = [at['v'] for at in partypes['atomtypes']]
    w_list = [at['w'] for at in partypes['atomtypes']]
    # get parameters
    v_i = v_list[ type_list.index(atoms[pair['ai']-1]['type']) ]
    w_i = w_list[ type_list.index(atoms[pair['ai']-1]['type']) ]
    v_j = v_list[ type_list.index(atoms[pair['aj']-1]['type']) ]
    w_j = w_list[ type_list.index(atoms[pair['aj']-1]['type']) ]
    if (comb_rule == 1) or (comb_rule == 3):
        v_ij = sqrt(v_i * v_j)
        w_ij = sqrt(w_i * w_j)
    elif (comb_rule == 2):
        v_ij = 0.5 * (v_i + v_j)
        w_ij = sqrt(w_i * w_j)
    
    # if 1-4 
    if (is_pair_in_pairs(pair, pairs)):
        q_ij *= fudgeQQ
        # if parameters were defined explicitly in .itp, use them
        if dof_has_parameters(pair):
            v_ij = pair['pars'][0]
            w_ij = pair['pars'][1]
        else:
            # if 'gen-pairs', use fudge for LJ
            if (gen_pairs == 'yes'):
                if (comb_rule == 1):
                    v_ij *= fudgeLJ
                    w_ij *= fudgeLJ
                else:
                    # multiply only epsilon
                    v_ij *= 1.0
                    w_ij *= fudgeLJ 
            # else, read from pairtypes
            elif (gen_pairs == 'no'):
                idx = -1
                for i,pt in enumerate(partypes['pairtypes']):
                    check_1 = [atoms[pair['ai']-1]['type'],atoms[pair['aj']-1]['type']] == [pt['ai'],pt['aj']]
                    check_2 = [atoms[pair['aj']-1]['type'],atoms[pair['ai']-1]['type']] == [pt['ai'],pt['aj']]
                    if (check_1) or (check_2):
                        idx = i
                        v_ij = partypes['pairtypes'][idx]['v']
                        w_ij = partypes['pairtypes'][idx]['w']
                        break
                # if not found, oops!
                if (idx == -1):
                    raise Exception("Error: no pairtype found for pair %s - %s and gen-pairs is set to 'no'." % (atoms[pair['ai']-1]['type'],atoms[pair['aj']-1]['type']))
    # if normal, override default parameters if defined in nonbond_params
    else:
        # search nonbond_params
        idx = -1
        for i,pt in enumerate(partypes['nonbond_params']):
            check_1 = [atoms[pair['ai']-1]['type'],atoms[pair['aj']-1]['type']] == [pt['ai'],pt['aj']]
            check_2 = [atoms[pair['aj']-1]['type'],atoms[pair['ai']-1]['type']] == [pt['ai'],pt['aj']]
            if (check_1) or (check_2):
                v_ij = partypes['nonbond_params'][i]['v']
                w_ij = partypes['nonbond_params'][i]['w']
                idx = i
                break
        
        

    # finally, assign
    pair['v'] = v_ij
    pair['w'] = w_ij
    pair['q'] = q_ij

def assign_pars_to_itp (itp, partypes, defaults):

    for pair in itp['pairs']:
        assign_pairtype_for_pair (pair, itp['atoms'], itp['pairs'], partypes, defaults['comb-rule'], defaults['gen-pairs'], defaults['fudgeLJ'], defaults['fudgeQQ'])

    for bond in itp['bonds']:
        assign_bondtype_for_bond (bond, itp['atoms'], partypes['bondtypes'])
        
    for angle in itp['angles']:
        assign_angletype_for_angle (angle, itp['atoms'], partypes['angletypes'])

    for dihedral in itp['dihedrals']:
        if not (assign_nongeneric_dihedraltype_for_dihedral (dihedral, itp['atoms'], partypes['dihedraltypes'])):
            assign_generic_dihedraltype_for_dihedral (dihedral, itp['atoms'], partypes['dihedraltypes'])

def get_atom_parameters (atomtype, ff_data):
    v = 0
    w = 0
    vs = 0
    ws = 0
    for x in ff_data['atomtypes']:
        if (x['name'] == atomtype):
            v = x['v']
            w = x['w']
            break
    for x in ff_data['pairtypes']:
        if (x['ai'] == atomtype) and (x['aj'] == atomtype):
            vs = x['v']
            ws = x['w']
            break
    return v, w, vs, ws

def write_defaults_to_stream(defaults, stream):
    stream.write("[ defaults ]\n")
    gen_pairs_dict = {'yes': 'fudge', 'no': 'no'}
    stream.write("; nbfunc comb-rule gen-pairs fudgeLJ fudgeQQ\n")
    stream.write("%8d%10d%10s%8.4f%8.4f\n\n" % (defaults['nbfunc'], defaults['comb-rule'], gen_pairs_dict[defaults['gen-pairs']], defaults['fudgeLJ'], defaults['fudgeQQ']))

def itp_to_stp (itp_fn, stp_fn, use_interface):
    if not (stp_fn.endswith('.stp')):
        raise Exception("output file must have extension .stp")
    if not (itp_fn.endswith('.itp')):
        raise Exception("input file must have extension .itp")

    fp = open(itp_fn, 'r')
    # read force field
    [defaults, ff_data]  = read_defaults_partypes_from_stream (fp)
    # read topology
    itp_data = read_itp_from_stream (fp, ff_data['atomtypes'])
    fp.close()

    # assign parameters based on force field and topology
    assign_pars_to_itp (itp_data, ff_data, defaults)

    # generate list of pairs taking exclusions into account
    atom_list_for_genpairs = [i + 1 for i, x in enumerate(itp_data['atoms'])]
    bond_list_for_genpairs = [[x['ai'], x['aj']] for x in itp_data['bonds']]
    pair_list = gen_pairs (atom_list_for_genpairs, bond_list_for_genpairs, itp_data['nexcl'], itp_data['exclusions'])
    pair_list_of_dicts = []
    # make list of pair data: [{'ai': ?, 'aj': ?, 'v': ?, 'w': ?, 'q': ?}, ...]
    for pl in pair_list:
        this_pair = {'ai':pl[0], 'aj':pl[1]}
        assign_pairtype_for_pair (this_pair, itp_data['atoms'], itp_data['pairs'], ff_data, defaults['comb-rule'], defaults['gen-pairs'], defaults['fudgeLJ'], defaults['fudgeQQ'])
        pair_list_of_dicts.append( this_pair )

    # now write to stp file
    fp = open(stp_fn, 'w')

    # defaults
    write_defaults_to_stream(defaults, fp)

    # atoms
    fp.write("[ atoms ]\n")
    for at in itp_data['atoms']:
        v, w, vs, ws = get_atom_parameters(at['type'], ff_data)
        if (defaults['gen-pairs'] == 'no'):
            fp.write("%-5s%18.7e%18.7e%18.7e%18.7e%18.5e ; %s \n" % (at['type'], v, w, vs, ws, at['q'], at['name']))
        else:
            fp.write("%-5s%18.7e%18.7e%18.5e ; %s \n" % (at['type'], v, w, at['q'], at['name']))
    fp.write("\n")

    fp.write("[ nbpairs ]\n")
    # first the standard pairs    
    idx = 1
    for pair in pair_list_of_dicts:
        tp_i = itp_data['atoms'][pair['ai']-1]['type']
        tp_j = itp_data['atoms'][pair['aj']-1]['type']
        fp.write("%5d%5d     1%18.8e%18.8e ; %s\n" % (pair['ai'], pair['aj'], pair['v'], pair['w'], "%d -> standard pair %s - %s" % (idx,tp_i,tp_j)))
        idx += 1
    # then the 1-4 pairs
    for pair in itp_data['pairs']:
        tp_i = itp_data['atoms'][pair['ai']-1]['type']
        tp_j = itp_data['atoms'][pair['aj']-1]['type']
        fp.write("%5d%5d     2%18.8e%18.8e ; %s\n" % (pair['ai'], pair['aj'], pair['v'], pair['w'], "%d -> 1-4 pair %s - %s" % (idx, tp_i, tp_j)))
        idx += 1
    fp.write("\n")

    # bonds
    fp.write("[ bonds ]\n")
    for bond in itp_data['bonds']:
        fp.write ("%5d%5d%5d%18.6e%18.6e" % (bond['ai'],bond['aj'],bond['func'],bond['pars']['b0'],bond['pars']['kb']))
        if 'comment' in bond:
            fp.write(" ; " + bond['comment'] + "\n")
        else:
            fp.write("\n")
    fp.write("\n")

    # angles
    fp.write("[ angles ]\n")
    for angle in itp_data['angles']:
        shape = get_angletype_shape(angle['func'])
        if (shape == 'harmonic'):
            fp.write ("%5d%5d%5d%5d%18.6e%18.6e" % (angle['ai'],angle['aj'],angle['ak'],angle['func'],angle['pars']['th0'],angle['pars']['k']))
        elif (shape == 'urey-bradley'):
            fp.write ("%5d%5d%5d%5d%18.6e%18.6e%18.6e%18.6e" % (angle['ai'],angle['aj'],angle['ak'],angle['func'],angle['pars']['th0'],angle['pars']['k'],
                                                                angle['pars']['r13'], angle['pars']['kub']))
        if 'comment' in angle:
            fp.write(" ; " + angle['comment'] + "\n")
        else:
            fp.write("\n")
    fp.write("\n")

    # proper dihedrals
    idx = 0
    fp.write("[ dihedrals ]\n")
    for dihedral in itp_data['dihedrals']:
        dihtype = get_dihedraltype_type(dihedral['func'])
        shape = get_dihedraltype_shape(dihedral['func'])
        npars = get_number_of_multiple_dihedrals(dihedral)
        if (dihtype == 'proper'):
            for i in range(npars):
                if shape == 'periodic':
                    fp.write ("%5d%5d%5d%5d%5d%18.6e%18.6e%5d" % (dihedral['ai'],dihedral['aj'],dihedral['ak'],dihedral['al'],dihedral['func'],
                            dihedral['pars']['phi0'][i],dihedral['pars']['k'][i],dihedral['pars']['m'][i]))
                    idx += 1
                elif shape == 'ryckaert':
                    fp.write ("%5d%5d%5d%5d%5d%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f" % (dihedral['ai'],dihedral['aj'],dihedral['ak'],dihedral['al'],dihedral['func'],
                            dihedral['pars']['c0'][i],
                            dihedral['pars']['c1'][i],
                            dihedral['pars']['c2'][i],
                            dihedral['pars']['c3'][i],
                            dihedral['pars']['c4'][i],
                            dihedral['pars']['c5'][i]))
                    idx += 1
                elif shape == 'fourier':
                    fp.write ("%5d%5d%5d%5d%5d%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f" % (dihedral['ai'],dihedral['aj'],dihedral['ak'],dihedral['al'],dihedral['func'],
                            dihedral['pars']['f1'][i],
                            dihedral['pars']['f2'][i],
                            dihedral['pars']['f3'][i],
                            dihedral['pars']['f4'][i]))
                    idx += 1
                try:
                    fp.write((" ; %d -> " % idx) + dihedral['comment'][i] + "\n")
                except IndexError:
                    fp.write("\n")
    fp.write("\n")

    # improper dihedrals
    fp.write("[ dihedrals ]\n")
    fp.write("; improper dihedrals\n")
    for dihedral in itp_data['dihedrals']:
        dihtype = get_dihedraltype_type(dihedral['func'])
        shape = get_dihedraltype_shape(dihedral['func'])
        npars = get_number_of_multiple_dihedrals(dihedral)
        if (dihtype == 'improper'):
            for i in range(npars):
                if (shape == 'harmonic'):
                    fp.write ("%5d%5d%5d%5d%5d%18.6e%18.6e" % (dihedral['ai'],dihedral['aj'],dihedral['ak'],dihedral['al'],dihedral['func'],
                        dihedral['pars']['phi0'][i], dihedral['pars']['k'][i]))
                elif (shape == 'periodic'):
                    fp.write ("%5d%5d%5d%5d%5d%18.6e%18.6e%5d" % (dihedral['ai'],dihedral['aj'],dihedral['ak'],dihedral['al'],dihedral['func'],
                            dihedral['pars']['phi0'][i],dihedral['pars']['k'][i],dihedral['pars']['m'][i]))
                try:
                    fp.write(" ; " + dihedral['comment'][i] + "\n")
                except IndexError:
                    fp.write("\n")
    fp.write("\n")
    if (use_interface):
        write_optblocks_via_interface(fp, itp_data)
    fp.close()

