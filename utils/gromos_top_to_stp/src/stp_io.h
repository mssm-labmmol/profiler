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

#ifndef _STP_IO_
#define _STP_IO_
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "topo.h"
#include "forcefield.h"
#include "regex_matcher.h"

/* from noMAD's answer in https://stackoverflow.com/questions/9562218/c-multiple-scanfs-when-i-enter-in-a-value-for-one-scanf-it-skips-the-second-s */
void seek_to_next_line( void )
{
    int c;
    while( (c = fgetc( stdin )) != EOF && c != '\n' );
}

void stp_write_atomic_data_to_stream (FILE* fp, struct topology_t* topo_data, struct forcefield_t* ff_data)
{
    const int resnr = 1;
    fprintf(fp, "[ atoms ]\n");
    fprintf(fp, "%-9s%-20s%-20s%-20s%-20s%9s\n", "; type", "c6", "c12", "cs6", "cs12", "q");
    for (int i = 0; i < topo_data->n_atoms; i++)
    {
        char type[10];
        struct pairtype_pars_t pairp;
        get_type_for_iac (ff_data, topo_data->iac[i], type);
        get_pairtype_parameters_for_code (ff_data, topo_data->iac[i], topo_data->iac[i], &pairp);
        fprintf(fp, "%-9s%-20.7e%-20.7e%-20.7e%-20.7e%9.4lf ; %d (%s)\n", type, pairp.c6, pairp.c12, pairp.cs6, pairp.cs12, topo_data->charges[i], i+1, topo_data->atomname[i]);
    }
    fprintf (fp, "\n");
    return;
}

int select_dihedrals_by_pattern (char* pattern, struct topology_t* topo_data, struct forcefield_t* ff_data, int* selected_idx)
{
    int count = 0;
    for (int i = 0; i < topo_data->n_dihedrals; i++)
    {
        char ni[6], nj[6], nk[6], nl[6];
        char ti[6], tj[6], tk[6], tl[6];
        strcpy(ni, topo_data->atomname[topo_data->dihedrals[i].ai - 1]);
        strcpy(nj, topo_data->atomname[topo_data->dihedrals[i].aj - 1]);
        strcpy(nk, topo_data->atomname[topo_data->dihedrals[i].ak - 1]);
        strcpy(nl, topo_data->atomname[topo_data->dihedrals[i].al - 1]);
        get_type_for_iac(ff_data, topo_data->iac[topo_data->dihedrals[i].ai - 1], ti);
        get_type_for_iac(ff_data, topo_data->iac[topo_data->dihedrals[i].aj - 1], tj);
        get_type_for_iac(ff_data, topo_data->iac[topo_data->dihedrals[i].ak - 1], tk);
        get_type_for_iac(ff_data, topo_data->iac[topo_data->dihedrals[i].al - 1], tl);
        if (compare_dih_with_pattern(ti, tj, tk, tl, pattern))
        {
            fprintf(stderr, "Automatically selected optdihedral %d -  %d-%d-%d-%d (%s-%s-%s-%s; %s-%s-%s-%s)\n",
                    i+1,
                    topo_data->dihedrals[i].ai, topo_data->dihedrals[i].aj, topo_data->dihedrals[i].ak, topo_data->dihedrals[i].al,
                    ni, nj, nk, nl,
                    ti, tj, tk, tl);
            selected_idx[count] = i+1;
            count++;
        }
    }
    return count;
}

int select_pairs_by_pattern (char* pattern, struct topology_t* topo_data, struct forcefield_t* ff_data, int* selected_idx)
{
    int count = 0;
    int parcount = 1;

    // first, loop over all possible pairs 
    for (int ai = 0; ai < topo_data->n_atoms; ai++)
    {
        for (int aj = ai+1; aj < topo_data->n_atoms; aj++)
        {
            // check if it is excluded
            int excl_flag = 0;
            for (int iexcl = 0; iexcl < topo_data->n_excl; iexcl++)
            {
                if (topo_data->exclusions[iexcl].ai == ai+1 &&
                        topo_data->exclusions[iexcl].aj == aj+1 )
                {
                    excl_flag = 1;
                    break;
                }
            }
            if (!excl_flag)
            {
	        // check if normal or 1-4
                int pair_flag = 0;
                for (int ipair = 0; ipair < topo_data->n_pairs; ipair++)
                {
                    if (topo_data->pairs[ipair].ai == ai+1 && 
                            topo_data->pairs[ipair].aj == aj+1)
                    {
                        pair_flag = 1;
                        break;
                    }
                }
                if (pair_flag)
                {
		  char ni[6], nj[6], ti[6], tj[6];
		  strcpy(ni, topo_data->atomname[ai]);
		  strcpy(nj, topo_data->atomname[aj]);
		  get_type_for_iac(ff_data, topo_data->iac[ai], ti);
		  get_type_for_iac(ff_data, topo_data->iac[aj], tj);
		  if (compare_pair_with_pattern(ti, tj, pattern)) {
		    fprintf(stderr, "Automatically selected 1-4 pair %d -  %d-%d (%s-%s; %s-%s)\n",
			    parcount,
			    ai+1, aj+1,
			    ni, nj,
			    ti, tj);
		    selected_idx[count] = parcount;
		    count++;
		  }
		}
		parcount++;
	    }
	}
    }
    return count;
}

int select_ref_by_inter (struct topology_t* topo_data, struct forcefield_t* ff_data)
{
    int answer;
    fprintf(stderr, "Which dihedral is your reference dihedral for the scan?\n");
    for (int i = 0; i < topo_data->n_dihedrals; i++)
    {
        char ni[6], nj[6], nk[6], nl[6];
        char ti[6], tj[6], tk[6], tl[6];
        strcpy(ni, topo_data->atomname[topo_data->dihedrals[i].ai - 1]);
        strcpy(nj, topo_data->atomname[topo_data->dihedrals[i].aj - 1]);
        strcpy(nk, topo_data->atomname[topo_data->dihedrals[i].ak - 1]);
        strcpy(nl, topo_data->atomname[topo_data->dihedrals[i].al - 1]);
        get_type_for_iac(ff_data, topo_data->iac[topo_data->dihedrals[i].ai - 1], ti);
        get_type_for_iac(ff_data, topo_data->iac[topo_data->dihedrals[i].aj - 1], tj);
        get_type_for_iac(ff_data, topo_data->iac[topo_data->dihedrals[i].ak - 1], tk);
        get_type_for_iac(ff_data, topo_data->iac[topo_data->dihedrals[i].al - 1], tl);
        fprintf(stderr, "%d -  %d-%d-%d-%d (%s-%s-%s-%s; %s-%s-%s-%s)\n",
                i+1,
                topo_data->dihedrals[i].ai, topo_data->dihedrals[i].aj, topo_data->dihedrals[i].ak, topo_data->dihedrals[i].al,
                ni, nj, nk, nl,
                ti, tj, tk, tl);
    }
    while (scanf(" %d", &answer) != EOF)
    {
        if ((answer > 0)&&(answer <= topo_data->n_dihedrals))
        {
            fprintf(stderr, "Selected %d.\n", answer);
            break;
        }
        else
            fprintf(stderr, "Invalid selection.\n");
    }
    seek_to_next_line();
    return answer;
}

int select_dihedrals_by_inter (struct topology_t* topo_data, struct forcefield_t* ff_data, int* selected_idx)
{
    int answer;
    int count = 0;
    fprintf(stderr, "Which dihedrals do you want to select for optimization?\n");
    fprintf(stderr, "(0 ends selection)\n\n");
    for (int i = 0; i < topo_data->n_dihedrals; i++)
    {
        char ni[6], nj[6], nk[6], nl[6];
        char ti[6], tj[6], tk[6], tl[6];
        strcpy(ni, topo_data->atomname[topo_data->dihedrals[i].ai - 1]);
        strcpy(nj, topo_data->atomname[topo_data->dihedrals[i].aj - 1]);
        strcpy(nk, topo_data->atomname[topo_data->dihedrals[i].ak - 1]);
        strcpy(nl, topo_data->atomname[topo_data->dihedrals[i].al - 1]);
        get_type_for_iac(ff_data, topo_data->iac[topo_data->dihedrals[i].ai - 1], ti);
        get_type_for_iac(ff_data, topo_data->iac[topo_data->dihedrals[i].aj - 1], tj);
        get_type_for_iac(ff_data, topo_data->iac[topo_data->dihedrals[i].ak - 1], tk);
        get_type_for_iac(ff_data, topo_data->iac[topo_data->dihedrals[i].al - 1], tl);
        fprintf(stderr, "%d -  %d-%d-%d-%d (%s-%s-%s-%s; %s-%s-%s-%s)\n",
                i+1,
                topo_data->dihedrals[i].ai, topo_data->dihedrals[i].aj, topo_data->dihedrals[i].ak, topo_data->dihedrals[i].al,
                ni, nj, nk, nl,
                ti, tj, tk, tl);
    }
    while (scanf(" %d", &answer) != EOF)
    {
        if (answer == 0)
            break;
        if ((answer > 0)&&(answer <= topo_data->n_dihedrals))
        {
            fprintf(stderr, "Selected %d.\n", answer);
            selected_idx[count] = answer;
            count++;
        }
        else
            fprintf(stderr, "Invalid selection.\n");
    }
    seek_to_next_line();
    return count;
}

int select_pairs_by_inter (struct topology_t* topo_data, struct forcefield_t* ff_data, int* selected_idx)
{
    int answer;
    int count = 0;
    int paircount = 1;
    fprintf(stderr, "Which 1-4 pairs do you want to select for optimization?\n");
    fprintf(stderr, "(0 ends selection)\n\n");
    // first, loop over all possible pairs 
    for (int ai = 0; ai < topo_data->n_atoms; ai++)
    {
        for (int aj = ai+1; aj < topo_data->n_atoms; aj++)
        {
            // check if it is excluded
            int excl_flag = 0;
            for (int iexcl = 0; iexcl < topo_data->n_excl; iexcl++)
            {
                if (topo_data->exclusions[iexcl].ai == ai+1 &&
                        topo_data->exclusions[iexcl].aj == aj+1 )
                {
                    excl_flag = 1;
                    break;
                }
            }
            if (!excl_flag)
            {
	        // check if normal or 1-4
                int pair_flag = 0;
                for (int ipair = 0; ipair < topo_data->n_pairs; ipair++)
                {
                    if (topo_data->pairs[ipair].ai == ai+1 && 
                            topo_data->pairs[ipair].aj == aj+1)
                    {
                        pair_flag = 1;
                        break;
                    }
                }
                if (pair_flag)
                {
		  char atomtype_i[FF_ATOMTYPE_BUFFER];
		  char atomtype_j[FF_ATOMTYPE_BUFFER];
		  get_type_for_iac(ff_data, topo_data->iac[ai], atomtype_i);
		  get_type_for_iac(ff_data, topo_data->iac[aj], atomtype_j);
		  fprintf(stderr, "%3d -   Name: %s-%s, Types: %s-%s, Idxs: %d-%d\n", paircount,
			  topo_data->atomname[ai],
			  topo_data->atomname[aj],
			  atomtype_i,
			  atomtype_j,
			  ai+1,
			  aj+1);
		}
		paircount++;
	    }
	}
    }
    while (scanf(" %d", &answer) != EOF)
    {
        if (answer == 0)
            break;
        if (answer > 0)
        {
            fprintf(stderr, "Selected %d.\n", answer);
            selected_idx[count] = answer;
            count++;
        }
        else
            fprintf(stderr, "Invalid selection.\n");
    }
    seek_to_next_line();
    return count;
}

int select_atoms_by_inter (struct topology_t* topo_data, struct forcefield_t* ff_data, int* selected_idx)
{
    int answer;
    int count = 0;
    fprintf(stderr, "Which atoms do you want to select for optimization?\n");
    fprintf(stderr, "(0 ends selection, -1 chooses pairs instead of atoms)\n\n");
    for (int i = 0; i < topo_data->n_atoms; i++)
    {
        char atomtype[FF_ATOMTYPE_BUFFER];
        get_type_for_iac(ff_data, topo_data->iac[i], atomtype);
        fprintf(stderr, "%d -  Name: %s, Type: %s\n", i+1, topo_data->atomname[i], atomtype);
    }
    while (scanf(" %d", &answer) != EOF)
    {
        if (answer == 0)
            break;
	if (answer == -1) {
	  count = -1;
	  break;
	}
	  
        if ((answer > 0)&&(answer <= topo_data->n_atoms))
        {
            fprintf(stderr, "Selected %d.\n", answer);
            selected_idx[count] = answer;
            count++;
        }
        else
            fprintf(stderr, "Invalid selection.\n");
    }
    seek_to_next_line();
    return count;
}

int select_atoms_by_type (char* type, struct topology_t* topo_data, struct forcefield_t* ff_data, int* selected_idx)
{
    int count = 0;
    for (int i = 0; i < topo_data->n_atoms; i++)
    {
        char atomtype[FF_ATOMTYPE_BUFFER];
        get_type_for_iac(ff_data, topo_data->iac[i], atomtype);
        if (!strcmp(atomtype,type))
        {
            selected_idx[count] = i+1;
            fprintf(stderr, "Automatically selected optatom %d, name = %s, type = %s.\n",
                    i+1, topo_data->atomname[i], atomtype);
            count++;
        }
    }
    return count;
}

void stp_write_pairs_data_to_stream (FILE* fp, struct topology_t* topo_data, struct forcefield_t* ff_data)
{
    const int funct = 1;
    int paircount = 1;
    fprintf (fp, "[ nbpairs ]\n");
    fprintf(fp, ";  ai   aj  type\n");
    // first, loop over all possible pairs 
    for (int ai = 0; ai < topo_data->n_atoms; ai++)
    {
        for (int aj = ai+1; aj < topo_data->n_atoms; aj++)
        {
            // check if it is excluded
            int excl_flag = 0;
            for (int iexcl = 0; iexcl < topo_data->n_excl; iexcl++)
            {
                if (topo_data->exclusions[iexcl].ai == ai+1 &&
                        topo_data->exclusions[iexcl].aj == aj+1 )
                {
                    excl_flag = 1;
                    break;
                }
            }
            if (!excl_flag)
            {
                double qij, c6ij, c12ij;
                char type_ai[FF_ATOMTYPE_BUFFER];
                char type_aj[FF_ATOMTYPE_BUFFER];

                struct pairtype_pars_t pairp;
                // get charges -- this does not depend on normal or 1-4
                // interactions
                qij = topo_data->charges[ai] * topo_data->charges[aj];
                get_pairtype_parameters_for_code (ff_data, topo_data->iac[ai], topo_data->iac[aj], &pairp);
                get_type_for_iac (ff_data, topo_data->iac[ai], type_ai);
                get_type_for_iac (ff_data, topo_data->iac[aj], type_aj);
		

                // check if normal or 1-4
                int pair_flag = 0;
                for (int ipair = 0; ipair < topo_data->n_pairs; ipair++)
                {
                    if (topo_data->pairs[ipair].ai == ai+1 && 
                            topo_data->pairs[ipair].aj == aj+1)
                    {
                        pair_flag = 1;
                        break;
                    }
                }
                if (pair_flag)
                {
                    c6ij = pairp.cs6;
                    c12ij = pairp.cs12;
                    fprintf(fp,"%5d%5d%5d%18.6e%18.6e     ; %d -> 1-4 pair %s-%s (%s-%s) \n", ai+1, aj+1, 2, paircount, c6ij, c12ij,
                            type_ai, type_aj, topo_data->atomname[ai], topo_data->atomname[aj]);
                }
                else
                {
                    c6ij = pairp.c6;
                    c12ij = pairp.c12;
                    fprintf(fp,"%5d%5d%5d%18.6e%18.6e     ; %d -> standard pair %s-%s (%s-%s) \n", ai+1, aj+1, 1, paircount, c6ij, c12ij,
			    type_ai, type_aj, topo_data->atomname[ai], topo_data->atomname[aj]);
                }

		paircount++;
            }
        }
    }
    fprintf (fp, "\n");
    return;
}


void stp_write_bond_data_to_stream (FILE* fp, struct topology_t* topo_data, struct forcefield_t* ff_data)
{
    const int funct = 2;
    fprintf (fp, "[ bonds ]\n");
    fprintf(fp, ";  ai   aj   type      l0          kb\n");
    for (int i = 0; i < topo_data->n_bonds; i++)
    {
        struct bond_pars_t bondp;
        /* Write bond parameters to a struct (defined in forcefield.c) */
        get_bond_parameters_for_code (ff_data, topo_data->bonds[i].code, &bondp);
        /* Write to stream */
        fprintf(fp,"%5d%5d%5d%12.4lf%12.4e ; type %d in GROMOS topology\n", topo_data->bonds[i].ai, topo_data->bonds[i].aj, funct, bondp.l, bondp.k, topo_data->bonds[i].code);
    }
    fprintf (fp, "\n");
    return;
}

void stp_write_angle_data_to_stream (FILE* fp, struct topology_t* topo_data, struct forcefield_t* ff_data)
{
    const int funct = 2;
    fprintf (fp, "[ angles ]\n");
    fprintf(fp, ";  ai   aj   ak   type   theta0      ka\n");
    for (int i = 0; i < topo_data->n_angles; i++)
    {
        struct angle_pars_t anglep;
        /* Write angle parameters to a struct (defined in forcefield.c) */
        get_angle_parameters_for_code (ff_data, topo_data->angles[i].code, &anglep);
        /* Write to stream */
        fprintf(fp,"%5d%5d%5d%5d%10.2lf%10.2lf ; type %d in GROMOS topology\n", topo_data->angles[i].ai, topo_data->angles[i].aj, topo_data->angles[i].ak, funct, anglep.theta, anglep.k,
                topo_data->angles[i].code);
    }
    fprintf (fp, "\n");
    return;
}

void stp_write_dihedral_data_to_stream (FILE* fp, struct topology_t* topo_data, struct forcefield_t* ff_data)
{
    /* For non-improper dihedrals, funct is 1. */
    const int funct = 1;
    fprintf (fp, "[ dihedrals ]\n");
    fprintf(fp, ";  ai   aj   ak   al   type      phi          k         mult\n");
    for (int i = 0; i < topo_data->n_dihedrals; i++)
    {
        struct dih_pars_t dihp;
        char type_ai[FF_ATOMTYPE_BUFFER];
        char type_aj[FF_ATOMTYPE_BUFFER];
        char type_ak[FF_ATOMTYPE_BUFFER];
        char type_al[FF_ATOMTYPE_BUFFER];

        get_type_for_iac (ff_data, topo_data->iac[topo_data->dihedrals[i].ai-1], type_ai);
        get_type_for_iac (ff_data, topo_data->iac[topo_data->dihedrals[i].aj-1], type_aj);
        get_type_for_iac (ff_data, topo_data->iac[topo_data->dihedrals[i].ak-1], type_ak);
        get_type_for_iac (ff_data, topo_data->iac[topo_data->dihedrals[i].al-1], type_al);

        /* Write dihedral parameters to a struct (defined in forcefield.c) */
        get_dihedral_parameters_for_code (ff_data, topo_data->dihedrals[i].code, &dihp);
        /* Write to stream */
        fprintf(fp,"%5d%5d%5d%5d%5d%12.4lf%12.4lf%9d ; type %d (%s-%s-%s-%s) in GROMOS topology\n", topo_data->dihedrals[i].ai, topo_data->dihedrals[i].aj, topo_data->dihedrals[i].ak,
               topo_data->dihedrals[i].al, funct, dihp.phi, dihp.k, dihp.mult, topo_data->dihedrals[i].code,
               type_ai, type_aj, type_ak, type_al);
    }
    fprintf (fp, "\n");
    return;
}

void stp_write_improper_data_to_stream (FILE* fp, struct topology_t* topo_data, struct forcefield_t* ff_data)
{
    /* For improper dihedrals, funct is 2. */
    const int funct = 2;
    fprintf (fp, "[ dihedrals ]\n; improper dihedrals\n");
    fprintf(fp, ";  ai   aj   ak   al   type      phi          k         \n");
    for (int i = 0; i < topo_data->n_impropers; i++)
    {
        struct idih_pars_t idihp;
        /* Write improper parameters to a struct (defined in forcefield.c) */
        get_improper_parameters_for_code (ff_data, topo_data->impropers[i].code, &idihp);
        /* Write to stream */
        fprintf(fp,"%5d%5d%5d%5d%5d%12.4lf%12.4lf ; type %d in GROMOS topology\n", topo_data->impropers[i].ai, topo_data->impropers[i].aj, topo_data->impropers[i].ak,
               topo_data->impropers[i].al, funct, idihp.phi, idihp.k, topo_data->impropers[i].code);
    }
    fprintf (fp, "\n");
    return;
}

void stp_write_selpairs_to_stream (FILE* fp, int* sel_pairs, int n_sel, struct topology_t* topo_data)
{
    fprintf(fp, "[ optpairs ]\n");
    fprintf(fp, "; idxs\n");
    for (int i = 0; i < n_sel; ++i)
    {
	int paircount = 0;
	// first, loop over all possible pairs 
	for (int ai = 0; ai < topo_data->n_atoms; ai++)
	  {
	      for (int aj = ai+1; aj < topo_data->n_atoms; aj++)
		{
		    // check if it is excluded
		    int excl_flag = 0;
		    for (int iexcl = 0; iexcl < topo_data->n_excl; iexcl++)
		      {
			  if (topo_data->exclusions[iexcl].ai == ai+1 &&
			      topo_data->exclusions[iexcl].aj == aj+1 )
			    {
				excl_flag = 1;
				break;
			    }
		      }
		    if (!excl_flag)
		      {
			  if (sel_pairs[i] == paircount + 1)
			    {
				fprintf(fp, "%-4d%-4d\n", ai+1, aj+1);
			    }
			  paircount++;
		      }
		}
	  }
	
    }
    fprintf(fp, "\n\n");
}

void stp_write_selatoms_to_stream (FILE* fp, int* sel_atoms, int n_sel)
{
    fprintf(fp, "[ optatoms ]\n");
    fprintf(fp, "; idxs\n");
    for (int i = 0; i < n_sel; ++i)
    {
        fprintf(fp, "%-4d", sel_atoms[i]);
    }
    fprintf(fp, "\n\n");
}

void stp_write_seldih_to_stream (FILE* fp, int* sel_dih, int n_sel, struct topology_t* topo_data)
{
    fprintf(fp, "[ optdihedrals ]\n");
    fprintf(fp, "; idxs\n");
    for (int i = 0; i < n_sel; ++i)
    {
	int ai = topo_data->dihedrals[sel_dih[i] - 1].ai;
	int aj = topo_data->dihedrals[sel_dih[i] - 1].aj;
	int ak = topo_data->dihedrals[sel_dih[i] - 1].ak;
	int al = topo_data->dihedrals[sel_dih[i] - 1].al;
        fprintf(fp, "%-4d%-4d%-4d%-4d\n", ai, aj, ak, al);
    }
    fprintf(fp, "\n\n");
}

static inline void stp_write_refdih_to_stream (FILE* fp, int refdih, struct topology_t* topo_data)
{
    int ai = topo_data->dihedrals[refdih - 1].ai;
    int aj = topo_data->dihedrals[refdih - 1].aj;
    int ak = topo_data->dihedrals[refdih - 1].ak;
    int al = topo_data->dihedrals[refdih - 1].al;
    fprintf(fp, "[ refdihedral ]\n; idx\n%-4d%-4d%-4d%-4d\n\n", ai, aj, ak, al);
}

static inline void stp_write_defaults_data_to_stream (FILE *fp)
{
    fprintf(fp, "[ defaults ]\n");
    fprintf(fp, "; nbfunc comb-rule gen-pairs fudgeLJ fudgeQQ\n");
    fprintf(fp, "       1         1    atomic     1.0     1.0\n\n");
}

void write_stp_to_stream (FILE* fp, struct topology_t* topo_data, struct forcefield_t* ff_data, int* sel_atoms, int n_sel, int* sel_dih, int n_sel_dih, int refdih, int optflag)
{
    /* Defaults info */
    stp_write_defaults_data_to_stream (fp);
    
    /* Atomic info -- */
    stp_write_atomic_data_to_stream (fp, topo_data, ff_data);

    stp_write_pairs_data_to_stream (fp, topo_data, ff_data);

    /* Bond info */
    stp_write_bond_data_to_stream (fp, topo_data, ff_data);

    /* Angle info */
    stp_write_angle_data_to_stream (fp, topo_data, ff_data);

    /* Dihedral info */
    stp_write_dihedral_data_to_stream (fp, topo_data, ff_data);

    /* Improper info */
    stp_write_improper_data_to_stream (fp, topo_data, ff_data);

    if (optflag) {
      /* Selected atoms */
      stp_write_selatoms_to_stream (fp, sel_atoms, n_sel);

      /* Selected dihedrals */
      stp_write_seldih_to_stream (fp, sel_dih, n_sel_dih, topo_data);

      /* Reference dihedral */
      stp_write_refdih_to_stream (fp, refdih, topo_data);
    }

    return;
}

void write_pairs_stp_to_stream (FILE* fp, struct topology_t* topo_data, struct forcefield_t* ff_data, int* sel_pairs, int n_sel, int* sel_dih, int n_sel_dih, int refdih, int optflag)
{
    /* Defaults info */
    stp_write_defaults_data_to_stream (fp);
    
    /* Atomic info -- */
    stp_write_atomic_data_to_stream (fp, topo_data, ff_data);

    stp_write_pairs_data_to_stream (fp, topo_data, ff_data);

    /* Bond info */
    stp_write_bond_data_to_stream (fp, topo_data, ff_data);

    /* Angle info */
    stp_write_angle_data_to_stream (fp, topo_data, ff_data);

    /* Dihedral info */
    stp_write_dihedral_data_to_stream (fp, topo_data, ff_data);

    /* Improper info */
    stp_write_improper_data_to_stream (fp, topo_data, ff_data);

    if (optflag) {
      /* Selected pairs */
      stp_write_selpairs_to_stream (fp, sel_pairs, n_sel, topo_data);

      /* Selected dihedrals */
      stp_write_seldih_to_stream (fp, sel_dih, n_sel_dih, topo_data);

      /* Reference dihedral */
      stp_write_refdih_to_stream (fp, refdih, topo_data);
    }

    return;
}

#endif
