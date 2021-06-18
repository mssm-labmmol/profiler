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

import copy


def is_distance_one(atom_1, atom_2, list_of_bonds):
    for bond in list_of_bonds:
        ld = [atom_1, atom_2]
        lr = ld[::-1]
        if (ld in list_of_bonds) or (lr in list_of_bonds):
            return True
    return False


def is_distance_two(atom_1, atom_2, list_of_bonds):
    bond_copy = copy.deepcopy(list_of_bonds)
    atoms_connected_to_1_or_2 = []
    for bond in bond_copy:
        if (atom_1 in bond) and (atom_2 in bond):
            continue
        if (atom_1 in bond):
            del (bond[bond.index(atom_1)])
            if (bond[0]) in atoms_connected_to_1_or_2:
                return True
            else:
                atoms_connected_to_1_or_2.append(bond[0])
        elif (atom_2 in bond):
            del (bond[bond.index(atom_2)])
            if (bond[0]) in atoms_connected_to_1_or_2:
                return True
            else:
                atoms_connected_to_1_or_2.append(bond[0])
    return False


def is_distance_three(atom_1, atom_2, list_of_bonds):
    if (is_distance_one(atom_1, atom_2, list_of_bonds)):
        return False
    if (is_distance_two(atom_1, atom_2, list_of_bonds)):
        return False
    for bond in list_of_bonds:
        if atom_1 == bond[0]:
            if is_distance_two(bond[1], atom_2, list_of_bonds):
                return True
        elif atom_1 == bond[1]:
            if is_distance_two(bond[0], atom_2, list_of_bonds):
                return True
        elif atom_2 == bond[0]:
            if is_distance_two(bond[1], atom_1, list_of_bonds):
                return True
        elif atom_2 == bond[1]:
            if is_distance_two(bond[0], atom_1, list_of_bonds):
                return True
    return False


def gen_pairs(list_of_atoms, list_of_bonds, nexcl, exclusions):
    """
    list_of_atoms:
    list_of_bonds:
    nexcl: int
    exclusions: [[ai_1, aj_1], [ai_2, aj_2], ...]"""
    if (nexcl > 3):
        raise Exception("error: only nexcl <= 3 supported")
    natoms = len(list_of_atoms)
    # generate list of all possible pairs and then remove the excluded
    # ones
    pairs = []
    for i in range(natoms):
        for j in range(i + 1, natoms):
            test_pair = [list_of_atoms[i], list_of_atoms[j]]
            if (test_pair not in exclusions):
                pairs.append([list_of_atoms[i], list_of_atoms[j]])
    funcs = [is_distance_one, is_distance_two, is_distance_three]
    for dist in range(nexcl):
        if (nexcl > dist - 1):
            idxs_to_be_removed = []
            for idx, [a1, a2] in enumerate(pairs):
                if (funcs[dist](a1, a2, list_of_bonds)):
                    idxs_to_be_removed.append(idx)
            pairs = [
                pairs[x] for x in range(len(pairs))
                if x not in idxs_to_be_removed
            ]
    return pairs
