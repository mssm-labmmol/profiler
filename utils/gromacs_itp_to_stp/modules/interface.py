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

class InterfaceError(Exception): pass

def input_loop(max=0, end='0', finish='-1', convert=lambda x: int(x)):
    output = []
    while True:
        answer = input("")
        if answer == end:
            return output, 'end' # this means the selection of a type has ended
        elif answer == finish:
            return output, 'finish' # this means the entire selection has ended
        else:
            output.append(convert(answer))
            if (max != 0) and (len(output) >= max):
                return output

def print_dihedral_list(itp_data):
    from .process import get_dihedraltype_type
    for i, dih in enumerate(itp_data['dihedrals']):
        if (get_dihedraltype_type(dih['func']) == 'proper'):
            names = {
                'ai': itp_data['atoms'][dih['ai'] - 1]['name'],
                'aj': itp_data['atoms'][dih['aj'] - 1]['name'],
                'ak': itp_data['atoms'][dih['ak'] - 1]['name'],
                'al': itp_data['atoms'][dih['al'] - 1]['name'],
            }
            types = {
                'ai': itp_data['atoms'][dih['ai'] - 1]['type'],
                'aj': itp_data['atoms'][dih['aj'] - 1]['type'],
                'ak': itp_data['atoms'][dih['ak'] - 1]['type'],
                'al': itp_data['atoms'][dih['al'] - 1]['type'],
            }
            print("\t{} - {}-{}-{}-{}, types: {}-{}-{}-{}, names: {}-{}-{}-{}".format(
                i+1, dih['ai'], dih['aj'], dih['ak'], dih['al'],
                types['ai'], types['aj'], types['ak'], types['al'],
                names['ai'], names['aj'], names['ak'], names['al']))

def write_refdihedral_block(stream, itp_data):
    print_dihedral_list(itp_data)
    print("Select the dihedrals to use as reference for the scan from the list above.")
    print("For each dihedral desired, type the corresponding number and then Enter.")
    print("  0 starts a new [ refdihedrals ] block.")
    print(" -1 ends the selection of reference dihedrals.")
    dihs, flag = input_loop()
    stream.write('[ refdihedrals ]\n')
    for dih in dihs:
        stream.write("{:<4d}{:<4d}{:<4d}{:<4d}\n".format(
            itp_data['dihedrals'][dih - 1]['ai'],
            itp_data['dihedrals'][dih - 1]['aj'],
            itp_data['dihedrals'][dih - 1]['ak'],
            itp_data['dihedrals'][dih - 1]['al']))
    stream.write('\n')
    if (flag == 'end'):
        # start a new type
        write_refdihedral_block(stream, itp_data)
    elif (flag == 'finish'):
        # finish this sub-interface
        return
    else:
        # unexpected! I don't know what to do!
        raise InterfaceError("Unexpected flag ", flag)

def write_optdihedrals_block(stream, itp_data):
    print_dihedral_list(itp_data)
    print("Select the dihedrals you want to optimize from the list below.")
    print("For each dihedral desired, type the corresponding number and then Enter.")
    print("  0 starts a new [ optdihedrals ] block.")
    print(" -1 ends the selection of dihedrals.")
    dihs, flag = input_loop()
    stream.write('[ optdihedrals ]\n')
    for dih in dihs:
        stream.write("{:<4d}{:<4d}{:<4d}{:<4d}\n".format(
            itp_data['dihedrals'][dih - 1]['ai'],
            itp_data['dihedrals'][dih - 1]['aj'],
            itp_data['dihedrals'][dih - 1]['ak'],
            itp_data['dihedrals'][dih - 1]['al']))
    stream.write('\n')
    if (flag == 'end'):
        # start a new type
        write_optdihedrals_block(stream, itp_data)
    elif (flag == 'finish'):
        # finish this sub-interface
        return
    else:
        # unexpected! I don't know what to do!
        raise InterfaceError("Unexpected flag ", flag)

def write_optpairs_block(stream, itp_data):
    for i, pair in enumerate(itp_data['pairs']):
        print("\t{} - {}-{}, types: {}-{}, names: {}-{}".format(
            i+1, pair['ai'], pair['aj'],
            itp_data['atoms'][pair['ai'] - 1]['type'], itp_data['atoms'][pair['aj'] - 1]['type'],
            itp_data['atoms'][pair['ai'] - 1]['name'], itp_data['atoms'][pair['aj'] - 1]['name']))
    print("Select the 1,4 pairs you want to optimize from the list above.")
    print("For each pair desired, type the corresponding number and then Enter.")
    print("  0 starts a new [ optpairs ] block.")
    print(" -1 ends the selection of pairs.")
    pairs, flag = input_loop()
    stream.write('[ optpairs ]\n')
    for pair in pairs:
        stream.write('{:<4d}{:<4d}\n'.format(
            itp_data['pairs'][pair - 1]['ai'], itp_data['pairs'][pair - 1]['aj']))
    stream.write('\n')
    if (flag == 'end'):
        # start a new type
        write_optpairs_block(stream, itp_data)
    elif (flag == 'finish'):
        # finish this sub-interface
        return
    else:
        # unexpected! I don't know what to do!
        raise InterfaceError("Unexpected flag ", flag)

def write_optatoms_block(stream, itp_data):
    for i, atom in enumerate(itp_data['atoms']):
        print("\t{} - type: {} name: {}".format(i+1, atom['type'], atom['name']))
    print("Select the atoms you want to optimize from the list above.")
    print("For each each atom desired, type its number and then Enter.")
    print("  0 starts a new [ optatoms ] block.")
    print(" -1 ends the selection of atoms.")
    atoms, flag = input_loop()
    stream.write('[ optatoms ]\n')
    for at in atoms:
        stream.write('{:<4d}'.format(int(at)))
    stream.write('\n\n')
    if (flag == 'end'):
        # start a new type
        write_optatoms_block(stream, itp_data)
    elif (flag == 'finish'):
        # finish this sub-interface
        return
    else:
        # unexpected! I don't know what to do!
        raise InterfaceError("Unexpected flag ", flag)

def choose_atoms_or_pairs(stream):
    while True:
        answer = input("Do you want to optimize an atom type (1), specific 1,4 pairs (2), or none of them (0)?\nType 0, 1 or 2, then Enter.\n")
        if (int(answer) == 1):
            return "atoms"
        if (int(answer) == 2):
            return "pairs"
        if (int(answer) == 0):
            return ""

def choose_refdihedrals(stream):
    while True:
        answer = input("Do you want to select reference dihedrals?\nType 0 (No), 1 (Yes), then Enter.\n")
        if (int(answer) == 1):
            return True
        if (int(answer) == 0):
            return False

def choose_optdihedrals(stream):
    while True:
        answer = input("Do you want to select dihedrals to optimize?\nType 0 (No), 1 (Yes), then Enter.\n")
        if (int(answer) == 1):
            return True
        if (int(answer) == 0):
            return False

def write_optblocks_via_interface(stream, itp_data):
    if (choose_refdihedrals(stream)):
        write_refdihedral_block(stream, itp_data)

    if (choose_optdihedrals(stream)):
        write_optdihedrals_block(stream, itp_data)

    choice = choose_atoms_or_pairs(itp_data)
    if (choice == 'atoms'):
        write_optatoms_block(stream, itp_data)
    elif (choice == 'pairs'):
        write_optpairs_block(stream, itp_data)
