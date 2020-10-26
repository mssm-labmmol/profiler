import numpy as np 

def parseGROtraj (fn):
    fp = open(fn, 'r')
    traj = []
    # file loop
    while (True):
        line = fp.readline() # title - ignored
        if (line == ''):
            break
        line = fp.readline() # number of atoms
        fields = line.split()
        natoms = int(fields[0])
        xs  = []
        ys  = []
        zs  = []
        els = []
        for i in range(natoms):
            line = fp.readline()
            # gro is fixed-format 
            ele = line[10:15].strip()[0]
            x = float(line[20:28])
            y = float(line[28:36])
            z = float(line[36:44])
            els.append(ele)
            xs.append(x)
            ys.append(y)
            zs.append(z)
        line = fp.readline() # box 
        traj.append(np.array([xs,ys,zs]))
    fp.close()
    return (traj, els)
