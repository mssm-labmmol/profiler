#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 16:02:50 2020

@author: Yan M. H. Gonçalves, IQ-UFRJ, Brazil
"""
import  argparse
from    sys                              import  argv,            stdout,          stderr
from    modules.readopts                 import  preprocess_args
from    modules.stpParser                import  parseStpFile
from    modules.configuration            import  *
from    modules.coordParser              import  calcDihedral,    wrapAngleDegree
from    modules.energy_force_vectorized  import  MMCalculator
from    modules.minim                    import  steepestDescentsMinimizer, conjugateGradientMinimizer
from    datetime                         import  datetime
from    platform                         import  platform,        python_version
from    os.path                          import  isfile
import  numpy                            as      np

def minimize_conf (conf: configuration, elements: list, stpdict: dict, dih_idxs: tuple, phi_0: float, k: float,
        min_algo: int, dx0: float, dxm: float, dele: float, nsteps: int,
        out_ener: str, out_traj: str) -> configuration:
    conf_copy = conf.copy()
    mmcalc = MMCalculator ()
    mmcalc.createFromStpDictionary(stpdict)
    mmcalc.pushDihedralRestraint(*dih_idxs, phi_0, k)
    if (min_algo == 1):
        minimizer = steepestDescentsMinimizer (dx0=dx0, dxm=dxm, dele=dele, maxSteps=nsteps, mmCalc=mmcalc)
    if (min_algo == 2):
        minimizer = conjugateGradientMinimizer (calc=mmcalc, dx0=dx0, dxm=dxm, nsteps=nsteps, prec=dele)
    minimizer.runAndSave(conf_copy, elements, out_ener, out_traj)
    return conf_copy

def printheader (fp):
    version = '1.0'
    prog    = 'profilerGen'
    bars    = "-"*(len(prog)+len(version)+6)
    string = "                         *** {0} (v. {1}) ***\n".format(prog,version) 
    string+= "                             {0}\n\n".format(bars)
    string+= "Running in {0} with Python {1}.\n".format(platform(), python_version())
    string+= "Started execution at {0}.\n".format(datetime.strftime(datetime.now(),"%Y/%m/%d %H:%M:%S")) 
    print(string, file=fp)

def printfooter (fp):
    string = "Finished execution at {0}.\n".format(datetime.strftime(datetime.now(),"%Y/%m/%d %H:%M:%S")) 
    print(string, file=fp)

def checkfiles (fnlist):
    for fn in fnlist:
        if not (isfile(fn)):
            raise IOError ("File {0} does not exist!".format(fn))

def debug_msg (msg, fp, verbose):
    if (verbose):
        print(msg, file=fp)

def err_msg (msg):
    raise RuntimeError(msg)

if __name__ == '__main__':
    progdescr = """
    profilerGen computes torsional-scan trajectories and their corresponding energy profiles.
    """
    parser = argparse.ArgumentParser(description=progdescr, formatter_class=argparse.RawTextHelpFormatter)
    # hack to avoid printing "optional arguments:" in help message
    parser._optionals.title = "options"
    parser.add_argument('-c', dest='coords', required=True, type=str, help=
            "Torsional-scan trajectory or single molecular conformation (.xyz/.gro/.g96).")
    parser.add_argument('-t', metavar='STP', dest='pars', required=True, type=str, help=
            "Special-topology file (.stp).")
    parser.add_argument('-dr', metavar=('RFRST','RSTEP','RLST'), dest='dih_restr', type=float, required=True, nargs=3, help=
            "Torsional-scan angles (deg): from RFRST to RLST with a RSTEP step.")
    parser.add_argument('-dk', metavar='RFCT', dest='dih_k', required=False, default=5000, help=
            "Force constant for dihedral restraint (default = 5000 kJ/(mol.rad^2)).")
    parser.add_argument('-op', metavar='PREFIX', dest='out', required=True, type=str, help=
            "Prefix for output files.")
    parser.add_argument('-min', metavar='MALG', dest='min_alg', required=False, type=int, default=1, help=
            "Minimization algorithm: steepest descents (1 - default) or conjugate-gradient (2).")
    parser.add_argument('-dx0', metavar='DX0', dest='min_dx0', required=False, type=float, default=0.05, help=
            "Initial step size DX0 in minimization algorithm (default = 0.05 nm).")
    parser.add_argument('-dxm', metavar='DXM', dest='min_dxm', required=False, type=float, default=0.2, help=
            "Maximal step size DXM in minimization algorithm (default = 0.20 nm).")
    parser.add_argument('-dele', metavar='DELE', dest='min_dele', required=False, type=float, default=1.0e-09, help=
            "If |E_n+1 - E_n| <= DELE, minimization stops (default = 1.0e-09).")
    parser.add_argument('-nsteps', metavar='NMAX', dest='min_nsteps', required=False, type=int, default=50000, help=
            "Maximal number of steps in minimization (default = 50000).")
    # this is in case arguments are provided in a -f <file>
    args = parser.parse_args(preprocess_args(argv[1:]))

    # starting program
    printheader(stdout)
    # check if files exist or raise IOError
    checkfiles ([args.coords, args.pars])
    # build angles list
    dihs = [args.dih_restr[0] + i*args.dih_restr[1] for i in range(int((args.dih_restr[2] - args.dih_restr[0])/args.dih_restr[1])+1)]
    # get reference dihedral value from stp file
    stpdict = parseStpFile(args.pars)
    refdih_idx = stpdict['refdihedral']
    dih_idxs = stpdict['propers'][0][refdih_idx][:4]
    ens = ensemble([])
    ens.readFromTrajectory(args.coords)
    if (ens.size() > 1):
        optMode = 'trajectory'
        if ens.size() != len(dihs):
            err_msg("Number of configurations in input trajectory does not match required dihedral scan.")
        debug_msg("Optimizing from input trajectory.", stdout, True)
    else:
        optMode = 'configuration'
        debug_msg("Optimizing from single configuration.", stdout, True)

    if optMode == 'configuration':
        orig_dih = ens[0].getDihedral(dih_idxs[0]-1, dih_idxs[1]-1, dih_idxs[2]-1, dih_idxs[3]-1)
        diff_dih = [np.abs(wrapAngleDegree(orig_dih - dih)) for dih in dihs]
        closest_dih = diff_dih.index(min(diff_dih))
        
        # output list of configurations
        out_conf = [None] * len(dihs)
        
        # forward - this includes closest_dih
        last_conf = ens[0]
        for i in range(closest_dih, len(dihs)):
            debug_msg("Minimizing {} -> angle = {}".format(i,dihs[i]), stdout, True)
            out_ener = "{}_{}.dat".format(args.out, i)
            out_traj = "{}_{}.xyz".format(args.out, i)
            out_conf[i] = minimize_conf (last_conf, ens.elements, stpdict, dih_idxs, dihs[i], args.dih_k,
                    args.min_alg, args.min_dx0, args.min_dxm, args.min_dele, args.min_nsteps,
                    out_ener, out_traj)
            last_conf = out_conf[i]
        
        last_conf = ens[0]
        # backward - does NOT include closest_dih
        for k in range(0, closest_dih):
            i = closest_dih - k - 1
            out_ener = "{}_{}.dat".format(args.out, i)
            out_traj = "{}_{}.xyz".format(args.out, i)
            debug_msg("Minimizing {} -> angle = {}".format(i,dihs[i]), stdout, True)
            out_conf[i] = minimize_conf (last_conf, ens.elements, stpdict, dih_idxs, dihs[i], args.dih_k,
                    args.min_alg, args.min_dx0, args.min_dxm, args.min_dele, args.min_nsteps,
                    out_ener, out_traj)
            last_conf = out_conf[i]
    elif optMode == 'trajectory':
        # output list of configurations
        out_conf = [None] * len(dihs)
        for i, dih in enumerate(dihs):
            out_ener = "{}_{}.dat".format(args.out, i)
            out_traj = "{}_{}.xyz".format(args.out, i)
            debug_msg("Minimizing {} -> angle = {}".format(i,dihs[i]), stdout, True)
            out_conf[i] = minimize_conf (ens[i], ens.elements, stpdict, dih_idxs, dihs[i], args.dih_k,
                    args.min_alg, args.min_dx0, args.min_dxm, args.min_dele, args.min_nsteps,
                    out_ener, out_traj)

    ens.confs = out_conf
    ens.writeToFile(args.out + '.xyz', fmt='xyz')

    mmcalc = MMCalculator ()
    mmcalc.createFromStpDictionary(stpdict)
    mmcalc.calcForEnsembleAndSaveToFile(ens, "{}.dat".format(args.out))
    mmcalc.calcForEnsembleAndSaveToFile(ens, "{}_full.dat".format(args.out), shiftToZero=False, saveOnlyTotal=False)

    # end program
    printfooter(stdout)
