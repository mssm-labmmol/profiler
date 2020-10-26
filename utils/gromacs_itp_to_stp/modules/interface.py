def input_loop(max=0, end='0', convert=lambda x: int(x)):
    output = []
    while True:
        answer = input("")
        if answer == end:
            return output
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
    print("Select the dihedral to use as reference for the scan from the list below.")
    print_dihedral_list(itp_data)
    dihs = input_loop(max=1)
    stream.write('[ refdihedral ]\n')
    for dih in dihs:
        stream.write("{:<4d}{:<4d}{:<4d}{:<4d}\n".format(
            itp_data['dihedrals'][dih - 1]['ai'],
            itp_data['dihedrals'][dih - 1]['aj'],
            itp_data['dihedrals'][dih - 1]['ak'],
            itp_data['dihedrals'][dih - 1]['al']))
    stream.write('\n')

def write_optdihedrals_block(stream, itp_data):
    print("Select the dihedrals you want to optimize from the list below.")
    print("For each dihedral desired, type the corresponding number and then Enter.")
    print("Typing 0 ends the selection.")
    print_dihedral_list(itp_data)
    dihs = input_loop()
    stream.write('[ optdihedrals ]\n')
    for dih in dihs:
        stream.write("{:<4d}{:<4d}{:<4d}{:<4d}\n".format(
            itp_data['dihedrals'][dih - 1]['ai'],
            itp_data['dihedrals'][dih - 1]['aj'],
            itp_data['dihedrals'][dih - 1]['ak'],
            itp_data['dihedrals'][dih - 1]['al']))
    stream.write('\n')

def write_optpairs_block(stream, itp_data):
    print("Select the 1,4 pairs you want to optimize from the list below.")
    print("For each pair desired, type the corresponding number and then Enter.")
    print("Typing 0 ends the selection.")
    for i, pair in enumerate(itp_data['nbpairs']):
        if (pair['func'] == 2): # 1,4 pair
            print("\t{} - {}-{}, types: {}-{}, names: {}-{}".format(
                i+1, pair['ai'], pair['aj'],
                itp_data['atoms'][pair['ai'] - 1]['type'], itp_data['atoms'][pair['aj'] - 1]['type'],
                itp_data['atoms'][pair['ai'] - 1]['name'], itp_data['atoms'][pair['aj'] - 1]['name']))
    pairs = input_loop()
    stream.write('[ optpairs ]\n')
    for pair in pairs:
        stream.write('{:<4d}{:<4d}\n'.format(
            itp_data['nbpairs'][pair - 1]['ai'], itp_data['nbpairs'][pair - 1]['aj']))
    stream.write('\n')

def write_optatoms_block(stream, itp_data):
    print("Select the atoms you want to optimize from the list below.")
    print("For each each atom desired, type its number and then Enter.")
    print("Typing 0 ends the selection.")
    for i, atom in enumerate(itp_data['atoms']):
        print("\t{} - type: {} name: {}".format(i+1, atom['type'], atom['name']))
    atoms = input_loop()
    stream.write('[ optatoms ]\n')
    for at in atoms:
        stream.write('{:<4d}'.format(at))
    stream.write('\n')

def choose_atoms_or_pairs(stream):
    while True:
        answer = input("Do you want to optimize an atom type (1) or specific 1,4 pairs (2)?\nType 1 or 2, then Enter.\n")
        if (int(answer) == 1):
            return "atoms"
        if (int(answer) == 2):
            return "pairs"

def write_optblocks_via_interface(stream, itp_data):
    write_refdihedral_block(stream, itp_data)
    write_optdihedrals_block(stream, itp_data)
    choice = choose_atoms_or_pairs(itp_data)
    if (choice == 'atoms'):
        write_optatoms_block(stream, itp_data)
    elif (choice == 'pairs'):
        write_optpairs_block(stream, itp_data)
    else:
        raise Exception("Unexpected choice: {}".format(choice))
