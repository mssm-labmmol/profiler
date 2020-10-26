import numpy as np 

def parseXYZtraj (fn):
    fp = open(fn, 'r')
    traj = []
    # file loop
    while (True):
        line = fp.readline()
        if (line == ''):
            break
        fields = line.split()
        natoms = int(fields[0])
        line = fp.readline() # title - ignored
        xs  = []
        ys  = []
        zs  = []
        els = []
        for i in range(natoms):
            line = fp.readline()
            fields = line.split()
            try:
                ele = int(fields[0])
                ele = z2name (ele)
            except ValueError:
                ele = fields[0]
                pass
            els.append(ele)
            x = float(fields[1]) / 10.0
            y = float(fields[2]) / 10.0
            z = float(fields[3]) / 10.0 
            xs.append(x)
            ys.append(y)
            zs.append(z)
        traj.append(np.array([xs,ys,zs]))
        
    fp.close()
    return (traj, els)

def z2name (ele):
    if ele == 1:
        return 'H'
    if ele == 2:
        return 'He'
    if ele == 3:
        return 'Li'
    if ele == 4:
        return 'Be'
    if ele == 5:
        return 'B'
    if ele == 6:
        return 'C'
    if ele == 7:
        return 'N'
    if ele == 8:
        return 'O'
    if ele == 9:
        return 'F'
    if ele == 10:
        return 'Ne'
    if ele == 11:
        return 'Na'
    if ele == 12:
        return 'Mg'
    if ele == 13:
        return 'Al'
    if ele == 14:
        return 'Si'
    if ele == 15:
        return 'P'
    if ele == 16:
        return 'S'
    if ele == 17:
        return 'Cl'
    if ele == 18:
        return 'Ar'
    if ele == 19:
        return 'K'
    if ele == 20:
        return 'Ca'
    if ele == 21:
        return 'Sc'
    if ele == 22:
        return 'Ti'
    if ele == 23:
        return 'V'
    if ele == 24:
        return 'Cr'
    if ele == 25:
        return 'Mn'
    if ele == 26:
        return 'Fe'
    if ele == 27:
        return 'Co'
    if ele == 28:
        return 'Ni'
    if ele == 29:
        return 'Cu'
    if ele == 30:
        return 'Zn'
    if ele == 31:
        return 'Ga'
    if ele == 32:
        return 'Ge'
    if ele == 33:
        return 'As'
    if ele == 34:
        return 'Se'
    if ele == 35:
        return 'Br'
    if ele == 36:
        return 'Kr'
    if ele == 37:
        return 'Rb'
    if ele == 38:
        return 'Sr'
    if ele == 39:
        return 'Y'
    if ele == 40:
        return 'Zr'
    if ele == 41:
        return 'Nb'
    if ele == 42:
        return 'Mo'
    if ele == 43:
        return 'Tc'
    if ele == 44:
        return 'Ru'
    if ele == 45:
        return 'Rh'
    if ele == 46:
        return 'Pd'
    if ele == 47:
        return 'Ag'
    if ele == 48:
        return 'Cd'
    if ele == 49:
        return 'In'
    if ele == 50:
        return 'Sn'
    if ele == 51:
        return 'Sb'
    if ele == 52:
        return 'Te'
    if ele == 53:
        return 'I'
    if ele == 54:
        return 'Xe'
    if ele == 55:
        return 'Cs'
    if ele == 56:
        return 'Ba'
    if ele == 57:
        return 'La'
    if ele == 58:
        return 'Ce'
    if ele == 59:
        return 'Pr'
    if ele == 60:
        return 'Nd'
    if ele == 61:
        return 'Pm'
    if ele == 62:
        return 'Sm'
    if ele == 63:
        return 'Eu'
    if ele == 64:
        return 'Gd'
    if ele == 65:
        return 'Tb'
    if ele == 66:
        return 'Dy'
    if ele == 67:
        return 'Ho'
    if ele == 68:
        return 'Er'
    if ele == 69:
        return 'Tm'
    if ele == 70:
        return 'Yb'
    if ele == 71:
        return 'Lu'
    if ele == 72:
        return 'Hf'
    if ele == 73:
        return 'Ta'
    if ele == 74:
        return 'W'
    if ele == 75:
        return 'Re'
    if ele == 76:
        return 'Os'
    if ele == 77:
        return 'Ir'
    if ele == 78:
        return 'Pt'
    if ele == 79:
        return 'Au'
    if ele == 80:
        return 'Hg'
    if ele == 81:
        return 'Tl'
    if ele == 82:
        return 'Pb'
    if ele == 83:
        return 'Bi'
    if ele == 84:
        return 'Po'
    if ele == 85:
        return 'At'
    if ele == 86:
        return 'Rn'
    if ele == 87:
        return 'Fr'
    if ele == 88:
        return 'Ra'
    if ele == 89:
        return 'Ac'
    if ele == 90:
        return 'Th'
    if ele == 91:
        return 'Pa'
    if ele == 92:
        return 'U'
    if ele == 93:
        return 'Np'
    if ele == 94:
        return 'Pu'
    if ele == 95:
        return 'Am'
    if ele == 96:
        return 'Cm'
    if ele == 97:
        return 'Bk'
    if ele == 98:
        return 'Cf'
    if ele == 99:
        return 'Es'
    if ele == 100:
        return 'Fm'
    if ele == 101:
        return 'Md'
    if ele == 102:
        return 'No'
    if ele == 103:
        return 'Lr'
    if ele == 104:
        return 'Rf'
    if ele == 105:
        return 'Db'
    if ele == 106:
        return 'Sg'
    if ele == 107:
        return 'Bh'
    if ele == 108:
        return 'Hs'
    if ele == 109:
        return 'Mt'
    if ele == 110:
        return 'Ds'
    if ele == 111:
        return 'Rg'
    if ele == 112:
        return 'Cn'
    if ele == 113:
        return 'Uut'
    if ele == 114:
        return 'Fl'
    if ele == 115:
        return 'Uup'
    if ele == 116:
        return 'Lv'
    if ele == 117:
        return 'Uus'
    if ele == 118:
        return 'Uuo'
