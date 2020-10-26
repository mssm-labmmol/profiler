//
// This file is part of the profilerTools suite (see
// https://github.com/mssm-labmmol/profiler).
//
// Copyright (c) 2020 mssm-labmmol
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#ifndef _TOPO_
#define _TOPO_
#include <stdio.h>

#define  RESNAME_MAX_SIZE   11
#define  N_ATOMS_MAX_SIZE   100
#define  N_BONDS_MAX_SIZE   400
#define  N_ANGLES_MAX_SIZE  400
#define  N_DIH_MAX_SIZE     100
#define  N_IDIH_MAX_SIZE    100
#define  N_EXCL_MAX_SIZE   1000
#define  N_PAIRS_MAX_SIZE   1000

/* Follows GROMOS syntax, where parameters are specified by a code.
 *
 * The parameters corresponding to each code can be retrieved from a 
 * force-field structure (defined/declared in forcefield.c). */

struct bond_t {
    int ai, aj, code;
};

struct angle_t {
    int ai, aj, ak, code;
};

struct dih_t {
    int ai, aj, ak, al, code;
};

struct idih_t {
    int ai, aj, ak, al, code;
};

/* exclusions */
struct excl_t {
    int ai, aj;
};

/* 1-4 neighbors considered in the GROMOS topology */
struct inei_t {
    int ai, aj;
};

struct topology_t {
    int    n_atoms;
    int    n_bonds;
    int    n_angles;
    int    n_dihedrals;
    int    n_impropers;
    int    n_pairs;
    int    n_excl;
    struct bond_t bonds[N_BONDS_MAX_SIZE]; 
    struct angle_t angles[N_ANGLES_MAX_SIZE];
    struct dih_t dihedrals[N_DIH_MAX_SIZE];
    struct idih_t impropers[N_IDIH_MAX_SIZE];
    struct excl_t exclusions[N_EXCL_MAX_SIZE];
    struct inei_t pairs[N_PAIRS_MAX_SIZE];
    double masses[N_ATOMS_MAX_SIZE];
    double charges[N_ATOMS_MAX_SIZE];
    int    iac[N_ATOMS_MAX_SIZE];
    int    cg_id[N_ATOMS_MAX_SIZE];
    char   resname[RESNAME_MAX_SIZE];
    char   atomname[N_ATOMS_MAX_SIZE][6];
};

/* convert cgid's from GROMOS standard to GROMACS standard  */
void correct_cgid (struct topology_t* topo_data)
{
    int current_group = 1;

    for (int i = 0; i < topo_data->n_atoms; i++)
    {
        int cg_id = topo_data->cg_id[i];

        if (cg_id == 0)
        {
            topo_data->cg_id[i] = current_group;
            continue;
        }
        else if (cg_id == 1)
        {
            topo_data->cg_id[i] = current_group;
            current_group++;
            continue;
        }
        else
        {
            fprintf(stderr, "error: found cgid different than 0 or 1, which is not GROMOS standard\n");
            exit(EXIT_FAILURE);
        }
    }

    return;
}


#endif
