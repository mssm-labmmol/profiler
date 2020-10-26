#ifndef _GROMOS_TOP_PARSER_
#define _GROMOS_TOP_PARSER_

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "topo.h"
#include "forcefield.h"

#define TOP_LINE_BUFFER 100
#define DEG_TO_RAD 0.0174533

/* puts next non-comment line in buffer variable */
/* returns -1 when it fails to read a line, 0 otherwise */
int top_read_next_line (FILE* fp, char* buffer)
{
    const char comment = '#';

    while (1)
    {
        if (!fgets(buffer, TOP_LINE_BUFFER, fp))
            return -1;
        if (buffer[0] != comment)
            return 0;
    }
}

/* checks if a line is the beginning of a block */
int is_block (char* buffer, const char* blockname)
{
    return !strncmp(buffer, blockname, strlen(blockname));
}

/* checks if a line is the end of a block */
int is_end_of_block (char* buffer)
{
    return !strncmp(buffer, "END", 3);
}

void parse_block_atomtypename (FILE* fp, topology_t* topo_data, forcefield_t* ff_data)
{
    int lines_read = 0;
    int n_atomtypes = -1;
    char line_buffer[TOP_LINE_BUFFER];

    while (1)
    {
        if (top_read_next_line (fp, line_buffer) != -1)
            lines_read++;

        else
        {
            fprintf(stderr, "unexpected end of file inside ATOMTYPENAME block\n");
            exit(EXIT_FAILURE);
        }

        if (lines_read == 1)
        {
            if (sscanf(line_buffer, "%d", &n_atomtypes) != 1)
            {
                fprintf(stderr, "error reading number of atomtypes in ATOMTYPENAME block\n");
                exit(EXIT_FAILURE);
            }
            ff_data->n_atomtypes = n_atomtypes;
        }

        if (lines_read > 1)
        {
            /* set iac */
            ff_data->atomtypes[lines_read-2].code = lines_read - 1;
            /* set type */
            sscanf(line_buffer, "%s", ff_data->atomtypes[lines_read-2].type);
            if (lines_read == n_atomtypes + 1)
            {
                /* finished parsing block */
                return;
            }
        }
    }
}

void parse_block_ljparameters  (FILE* fp, topology_t* topo_data, forcefield_t* ff_data)
{
    int lines_read = 0;
    int n_pairtypes = -1;
    char line_buffer[TOP_LINE_BUFFER];

    while (1)
    {
        if (top_read_next_line (fp, line_buffer) != -1)
            lines_read++;

        else
        {
            fprintf(stderr, "unexpected end of file inside LJPARAMETERS block\n");
            exit(EXIT_FAILURE);
        }

        if (lines_read == 1)
        {
            if (sscanf(line_buffer, "%d", &n_pairtypes) != 1)
            {
                fprintf(stderr, "error reading number of pair types in LJPARAMETERS block\n");
                exit(EXIT_FAILURE);
            }
            ff_data->n_pairtypes = n_pairtypes;
        }

        if (lines_read > 1)
        {
            int idx = lines_read - 2;
            /* set parameters */
            if (sscanf(line_buffer, "%d%d%lf%lf%lf%lf", &ff_data->pairtypes[idx].iac,&ff_data->pairtypes[idx].jac,
                        &ff_data->pairtypes[idx].c12, &ff_data->pairtypes[idx].c6,
                        &ff_data->pairtypes[idx].cs12,&ff_data->pairtypes[idx].cs6) != 6)
            {
                fprintf(stderr, "error: reading pair type in LJPARAMETERS block. unexpected format in line \n%s", line_buffer);
                exit(EXIT_FAILURE);
            }

            if (lines_read == n_pairtypes + 1)
            {
                /* finished parsing block */
                return;
            }
        }
    }
}

void parse_block_resname (FILE* fp, topology_t* topo_data, forcefield_t* ff_data)
{
    int lines_read = 0;
    int n_res = -1;
    char line_buffer[TOP_LINE_BUFFER];

    while (1)
    {
        if (top_read_next_line (fp, line_buffer) != -1)
            lines_read++;

        else
        {
            fprintf(stderr, "unexpected end of file inside RESNAME block\n");
            exit(EXIT_FAILURE);
        }

        if (lines_read == 1)
        {
            if (sscanf(line_buffer, "%d", &n_res) != 1)
            {
                fprintf(stderr, "error reading number of residues in RESNAME block\n");
                exit(EXIT_FAILURE);
            }

            /* only one residue supported */
            if (n_res != 1)
            {
                fprintf(stderr, "error: only one residue supported; your topology has %d!\n", n_res);
                exit(EXIT_FAILURE);
            }
        }

        if (lines_read == 2)
        {
            /* set resname */
            sscanf(line_buffer, "%s", topo_data->resname);

            /* finished parsing block */
            return;
        }
    }
}

void parse_block_soluteatom (FILE* fp, topology_t* topo_data, forcefield_t* ff_data)
{
    /* for this block, lines_read is actually "atoms read", since an atom can be
     * defined by more than one line */
    int lines_read = 0;
    int n_atoms = -1;
    char line_buffer[TOP_LINE_BUFFER];

    /* initialize n_pairs and n_excl just in case */
    topo_data->n_pairs = 0;
    topo_data->n_excl = 0;

    while (1)
    {
        if (top_read_next_line (fp, line_buffer) != -1)
            lines_read++;

        else
        {
            fprintf(stderr, "unexpected end of file inside SOLUTEATOM block\n");
            exit(EXIT_FAILURE);
        }

        if (lines_read == 1)
        {
            if (sscanf(line_buffer, "%d", &n_atoms) != 1)
            {
                fprintf(stderr, "error reading number of atoms in SOLUTEATOM block\n");
                exit(EXIT_FAILURE);
            }
            topo_data->n_atoms = n_atoms;
        }

        if (lines_read > 1)
        {
            int atnm = lines_read - 1;
            int atidx = lines_read - 2;
            int non_nei_read = -1;
            int nei_read = -1;
            int n_nei;
            int n_inei;
            int nei_offset;
            int inei_offset;

            /* set non-neighbor things */
            non_nei_read = sscanf(line_buffer, "%*d%*d%s%d%lf%lf%d%d%n", topo_data->atomname[atidx], &(topo_data->iac[atidx]), &(topo_data->masses[atidx]), &(topo_data->charges[atidx]), &(topo_data->cg_id[atidx]), &(n_nei), &nei_offset);
            if (non_nei_read != 6)
            {
                fprintf(stderr, "error: something went wrong reading SOLUTEATOM block in the following line:\n%s\n", line_buffer);
                exit(EXIT_FAILURE);
            }

            /* read exclusions in this line and the ones that follow */
            while (n_nei > 0)
            {
                int off;
                int excl;
                if (sscanf(line_buffer + nei_offset + 1, "%d%n", &excl, &off) == 1)
                {
                    nei_offset += off;
                    n_nei = n_nei - 1;
                    /* update topology */
                    topo_data->n_excl++;
                    topo_data->exclusions[topo_data->n_excl-1].ai = atnm;
                    topo_data->exclusions[topo_data->n_excl-1].aj = excl;
                }
                else
                {
                    if (top_read_next_line (fp, line_buffer) == -1)
                    {
                        fprintf(stderr, "error: unexpected end of topology inside SOLUTEATOM block\n");
                        exit(EXIT_FAILURE);
                    }
                    nei_offset = 0;
                }
            }
            /* get 1-4 interactions */
            /* read next line and do the same thing as with exclusions, 
             * but keeping track of the pairs */
            if (top_read_next_line (fp, line_buffer) == -1)
            {
                fprintf(stderr, "error: unexpected end of topology inside SOLUTEATOM block\n");
                exit(EXIT_FAILURE);
            }
            if (sscanf(line_buffer, "%d%n", &n_inei, &inei_offset) != 1)
            {
                fprintf(stderr, "error reading SOLUTEATOM block; unexpected format in line %s\n", line_buffer);
                exit(EXIT_FAILURE);
            }
            while (n_inei > 0)
            {
                int off;
                int pair;
                if (sscanf(line_buffer + inei_offset + 1, "%d%n", &pair, &off) == 1)
                {
                    inei_offset += off;
                    n_inei = n_inei - 1;
                    /* update topology */
                    topo_data->n_pairs++;
                    topo_data->pairs[topo_data->n_pairs-1].ai = atnm;
                    topo_data->pairs[topo_data->n_pairs-1].aj = pair;
                }
                else
                {
                    if (top_read_next_line (fp, line_buffer) == -1)
                    {
                        fprintf(stderr, "error: unexpected end of topology inside SOLUTEATOM block\n");
                        exit(EXIT_FAILURE);
                    }
                    inei_offset = 0;
                }
            }
            
            if (lines_read == n_atoms + 1)
            {
                /* finished parsing block */
                /* before returning, correct cgid to GROMACS standard */
                correct_cgid (topo_data);
                return;
            }
        }
    }
}

void parse_block_bondstretchtype (FILE* fp, topology_t* topo_data, forcefield_t* ff_data)
{
    int lines_read = 0;
    int n_bondtypes = -1;
    char line_buffer[TOP_LINE_BUFFER];

    while (1)
    {
        if (top_read_next_line (fp, line_buffer) != -1)
            lines_read++;

        else
        {
            fprintf(stderr, "unexpected end of file inside BONDSTRETCHTYPE block\n");
            exit(EXIT_FAILURE);
        }

        if (lines_read == 1)
        {
            if (sscanf(line_buffer, "%d", &n_bondtypes) != 1)
            {
                fprintf(stderr, "error reading number of bond types in BONDSTRETCHTYPE block\n");
                exit(EXIT_FAILURE);
            }
            ff_data->n_bondtypes = n_bondtypes;
        }

        if (lines_read > 1)
        {
            int idx = lines_read - 2;
            /* set code */
            ff_data->bondtypes[idx].code = idx + 1;
            /* set length and force constant */
            /* for bonds, GROMOS and GROMACS units are the same */
            if (sscanf(line_buffer, "%lf%*f%lf", &ff_data->bondtypes[idx].k, &ff_data->bondtypes[idx].l) != 2)
            {
                fprintf(stderr, "error: reading bondtype in BONDSTRETCHTYPE block. unexpected format in line \n%s", line_buffer);
                exit(EXIT_FAILURE);
            }

            if (lines_read == n_bondtypes + 1)
            {
                /* finished parsing block */
                return;
            }
        }
    }
}

void parse_block_bondh (FILE* fp, topology_t* topo_data, forcefield_t* ff_data)
{
    int lines_read = 0;
    int n_bondh = -1;
    char line_buffer[TOP_LINE_BUFFER];

    while (1)
    {
        if (top_read_next_line (fp, line_buffer) != -1)
            lines_read++;

        else
        {
            fprintf(stderr, "unexpected end of file inside BONDH block\n");
            exit(EXIT_FAILURE);
        }

        if (lines_read == 1)
        {
            if (sscanf(line_buffer, "%d", &n_bondh) != 1)
            {
                fprintf(stderr, "error reading number of bonds in BONDH block\n");
                exit(EXIT_FAILURE);
            }
            if (n_bondh == 0)
                return;
        }

        if (lines_read > 1)
        {
            topo_data->n_bonds++;
            int idx = topo_data->n_bonds - 1;

            if (sscanf(line_buffer, "%d%d%d", &topo_data->bonds[idx].ai, &topo_data->bonds[idx].aj, &topo_data->bonds[idx].code) != 3)
            {
                fprintf(stderr, "error: reading bond in BONDH block. unexpected format in line \n%s", line_buffer);
                exit(EXIT_FAILURE);
            }

            if (lines_read == n_bondh + 1)
            {
                /* finished parsing block */
                return;
            }
        }
    }
}

void parse_block_bond (FILE* fp, topology_t* topo_data, forcefield_t* ff_data)
{
    int lines_read = 0;
    int n_bond = -1;
    char line_buffer[TOP_LINE_BUFFER];

    while (1)
    {
        if (top_read_next_line (fp, line_buffer) != -1)
            lines_read++;

        else
        {
            fprintf(stderr, "unexpected end of file inside BOND block\n");
            exit(EXIT_FAILURE);
        }

        if (lines_read == 1)
        {
            if (sscanf(line_buffer, "%d", &n_bond) != 1)
            {
                fprintf(stderr, "error reading number of bonds in BOND block\n");
                exit(EXIT_FAILURE);
            }
            if (n_bond == 0)
                return;
        }

        if (lines_read > 1)
        {
            topo_data->n_bonds++;
            int idx = topo_data->n_bonds - 1;

            if (sscanf(line_buffer, "%d%d%d", &topo_data->bonds[idx].ai, &topo_data->bonds[idx].aj, &topo_data->bonds[idx].code) != 3)
            {
                fprintf(stderr, "error: reading bond in BONDH block. unexpected format in line \n%s", line_buffer);
                exit(EXIT_FAILURE);
            }

            if (lines_read == n_bond + 1)
            {
                /* finished parsing block */
                return;
            }
        }
    }
}

void parse_block_bondanglebendtype (FILE* fp, topology_t* topo_data, forcefield_t* ff_data)
{
    int lines_read = 0;
    int n_bondanglebendtypes = -1;
    char line_buffer[TOP_LINE_BUFFER];

    while (1)
    {
        if (top_read_next_line (fp, line_buffer) != -1)
            lines_read++;

        else
        {
            fprintf(stderr, "unexpected end of file inside BONDANGLEBENDTYPE block\n");
            exit(EXIT_FAILURE);
        }

        if (lines_read == 1)
        {
            if (sscanf(line_buffer, "%d", &n_bondanglebendtypes) != 1)
            {
                fprintf(stderr, "error reading number of bond angle types in BONDANGLEBENDTYPE block\n");
                exit(EXIT_FAILURE);
            }
            ff_data->n_angletypes = n_bondanglebendtypes;
        }

        if (lines_read > 1)
        {
            int idx = lines_read - 2;
            /* set code */
            ff_data->angletypes[idx].code = idx + 1;
            /* set angle and force constant */
            /* for angles, GROMOS and GROMACS units are the same */
            if (sscanf(line_buffer, "%lf%*f%lf", &ff_data->angletypes[idx].k, &ff_data->angletypes[idx].theta) != 2)
            {
                fprintf(stderr, "error: reading angletype in BONDANGLEBENDTYPE block. unexpected format in line \n%s", line_buffer);
                exit(EXIT_FAILURE);
            }

            if (lines_read == n_bondanglebendtypes + 1)
            {
                /* finished parsing block */
                return;
            }
        }
    }
}

void parse_block_bondangleh (FILE* fp, topology_t* topo_data, forcefield_t* ff_data)
{
    int lines_read = 0;
    int n_bondangleh = -1;
    char line_buffer[TOP_LINE_BUFFER];

    while (1)
    {
        if (top_read_next_line (fp, line_buffer) != -1)
            lines_read++;

        else
        {
            fprintf(stderr, "unexpected end of file inside BONDANGLEH block\n");
            exit(EXIT_FAILURE);
        }

        if (lines_read == 1)
        {
            if (sscanf(line_buffer, "%d", &n_bondangleh) != 1)
            {
                fprintf(stderr, "error reading number of angles in BONDANGLEH block\n");
                exit(EXIT_FAILURE);
            }
            if (n_bondangleh == 0)
                return;
        }

        if (lines_read > 1)
        {
            topo_data->n_angles++;
            int idx = topo_data->n_angles - 1;

            if (sscanf(line_buffer, "%d%d%d%d", &topo_data->angles[idx].ai, &topo_data->angles[idx].aj, &topo_data->angles[idx].ak, &topo_data->angles[idx].code) != 4)
            {
                fprintf(stderr, "error: reading angle in BONDANGLEH block. unexpected format in line \n%s", line_buffer);
                exit(EXIT_FAILURE);
            }

            if (lines_read == n_bondangleh + 1)
            {
                /* finished parsing block */
                return;
            }
        }
    }
}

void parse_block_bondangle (FILE* fp, topology_t* topo_data, forcefield_t* ff_data)
{
    int lines_read = 0;
    int n_bondangle = -1;
    char line_buffer[TOP_LINE_BUFFER];

    while (1)
    {
        if (top_read_next_line (fp, line_buffer) != -1)
            lines_read++;

        else
        {
            fprintf(stderr, "unexpected end of file inside BONDANGLE block\n");
            exit(EXIT_FAILURE);
        }

        if (lines_read == 1)
        {
            if (sscanf(line_buffer, "%d", &n_bondangle) != 1)
            {
                fprintf(stderr, "error reading number of angles in BONDANGLE block\n");
                exit(EXIT_FAILURE);
            }
            if (n_bondangle == 0)
                return;
        }

        if (lines_read > 1)
        {
            topo_data->n_angles++;
            int idx = topo_data->n_angles - 1;

            if (sscanf(line_buffer, "%d%d%d%d", &topo_data->angles[idx].ai, &topo_data->angles[idx].aj, &topo_data->angles[idx].ak, &topo_data->angles[idx].code) != 4)
            {
                fprintf(stderr, "error: reading angle in BONDANGLE block. unexpected format in line \n%s", line_buffer);
                exit(EXIT_FAILURE);
            }

            if (lines_read == n_bondangle + 1)
            {
                /* finished parsing block */
                return;
            }
        }
    }
}

void parse_block_impdihedraltype (FILE* fp, topology_t* topo_data, forcefield_t* ff_data)
{
    int lines_read = 0;
    int n_impdihedraltypes = -1;
    char line_buffer[TOP_LINE_BUFFER];

    while (1)
    {
        if (top_read_next_line (fp, line_buffer) != -1)
            lines_read++;

        else
        {
            fprintf(stderr, "unexpected end of file inside IMPDIHEDRALTYPE block\n");
            exit(EXIT_FAILURE);
        }

        if (lines_read == 1)
        {
            if (sscanf(line_buffer, "%d", &n_impdihedraltypes) != 1)
            {
                fprintf(stderr, "error reading number of improper dihedral types in IMPDIHEDRALTYPE block\n");
                exit(EXIT_FAILURE);
            }
            ff_data->n_idihtypes = n_impdihedraltypes;
        }

        if (lines_read > 1)
        {
            int idx = lines_read - 2;
            /* set code */
            ff_data->idihtypes[idx].code = idx + 1;
            /* set angle and force constant */
            /* for improper dihedrals, GROMOS and GROMACS units are NOT the same for the force constant */
            if (sscanf(line_buffer, "%lf%lf", &ff_data->idihtypes[idx].k, &ff_data->idihtypes[idx].phi) != 2)
            {
                fprintf(stderr, "error: reading impropertype in IMPDIHEDRALTYPE block. unexpected format in line \n%s", line_buffer);
                exit(EXIT_FAILURE);
            }
            /* convert force constant */
            ff_data->idihtypes[idx].k /= (DEG_TO_RAD)*(DEG_TO_RAD);

            if (lines_read == n_impdihedraltypes + 1)
            {
                /* finished parsing block */
                return;
            }
        }
    }
}

void parse_block_impdihedralh (FILE* fp, topology_t* topo_data, forcefield_t* ff_data)
{
    int lines_read = 0;
    int n_impdihedralh = -1;
    char line_buffer[TOP_LINE_BUFFER];

    while (1)
    {
        if (top_read_next_line (fp, line_buffer) != -1)
            lines_read++;

        else
        {
            fprintf(stderr, "unexpected end of file inside IMPDIHEDRALH block\n");
            exit(EXIT_FAILURE);
        }

        if (lines_read == 1)
        {
            if (sscanf(line_buffer, "%d", &n_impdihedralh) != 1)
            {
                fprintf(stderr, "error reading number of improper dihedrals in IMPDIHEDRALH block\n");
                exit(EXIT_FAILURE);
            }
            if (n_impdihedralh == 0)
                return;
        }

        if (lines_read > 1)
        {
            topo_data->n_impropers++;
            int idx = topo_data->n_impropers - 1;

            if (sscanf(line_buffer, "%d%d%d%d%d", &topo_data->impropers[idx].ai, &topo_data->impropers[idx].aj, &topo_data->impropers[idx].ak, &topo_data->impropers[idx].al, &topo_data->impropers[idx].code) != 5)
            {
                fprintf(stderr, "error: reading improper in block IMPDIHEDRALH. unexpected format in line \n%s", line_buffer);
                exit(EXIT_FAILURE);
            }

            if (lines_read == n_impdihedralh + 1)
            {
                /* finished parsing block */
                return;
            }
        }
    }
}

void parse_block_impdihedral (FILE* fp, topology_t* topo_data, forcefield_t* ff_data)
{
    int lines_read = 0;
    int n_impdihedral = -1;
    char line_buffer[TOP_LINE_BUFFER];

    while (1)
    {
        if (top_read_next_line (fp, line_buffer) != -1)
            lines_read++;

        else
        {
            fprintf(stderr, "unexpected end of file inside IMPDIHEDRAL block\n");
            exit(EXIT_FAILURE);
        }

        if (lines_read == 1)
        {
            if (sscanf(line_buffer, "%d", &n_impdihedral) != 1)
            {
                fprintf(stderr, "error reading number of improper dihedrals in IMPDIHEDRAL block\n");
                exit(EXIT_FAILURE);
            }
            if (n_impdihedral == 0)
                return;
        }

        if (lines_read > 1)
        {
            topo_data->n_impropers++;
            int idx = topo_data->n_impropers - 1;

            if (sscanf(line_buffer, "%d%d%d%d%d", &topo_data->impropers[idx].ai, &topo_data->impropers[idx].aj, &topo_data->impropers[idx].ak, &topo_data->impropers[idx].al, &topo_data->impropers[idx].code) != 5)
            {
                fprintf(stderr, "error: reading improper in block IMPDIHEDRAL. unexpected format in line \n%s", line_buffer);
                exit(EXIT_FAILURE);
            }

            if (lines_read == n_impdihedral + 1)
            {
                /* finished parsing block */
                return;
            }
        }
    }
}

void parse_block_dihedraltype (FILE* fp, topology_t* topo_data, forcefield_t* ff_data)
{
    int lines_read = 0;
    int n_dihedraltypes = -1;
    char line_buffer[TOP_LINE_BUFFER];

    while (1)
    {
        if (top_read_next_line (fp, line_buffer) != -1)
            lines_read++;

        else
        {
            fprintf(stderr, "unexpected end of file inside TORSDIHEDRALTYPE block\n");
            exit(EXIT_FAILURE);
        }

        if (lines_read == 1)
        {
            if (sscanf(line_buffer, "%d", &n_dihedraltypes) != 1)
            {
                fprintf(stderr, "error reading number of dihedral types in TORSDIHEDRALTYPE block\n");
                exit(EXIT_FAILURE);
            }
            ff_data->n_dihtypes = n_dihedraltypes;
        }

        if (lines_read > 1)
        {
            int idx = lines_read - 2;
            /* set code */
            ff_data->dihtypes[idx].code = idx + 1;
            /* set angle, force constant and multiplicity */
            /* for dihedrals, GROMOS and GROMACS units are the same */
            if (sscanf(line_buffer, "%lf%lf%d", &ff_data->dihtypes[idx].k, &ff_data->dihtypes[idx].phi, &ff_data->dihtypes[idx].mult) != 3)
            {
                fprintf(stderr, "error: reading dihedral type in TORSDIHEDRALTYPE block. unexpected format in line \n%s", line_buffer);
                exit(EXIT_FAILURE);
            }

            if (lines_read == n_dihedraltypes + 1)
            {
                /* finished parsing block */
                return;
            }
        }
    }
}

void parse_block_dihedralh (FILE* fp, topology_t* topo_data, forcefield_t* ff_data)
{
    int lines_read = 0;
    int n_dihedralh = -1;
    char line_buffer[TOP_LINE_BUFFER];

    while (1)
    {
        if (top_read_next_line (fp, line_buffer) != -1)
            lines_read++;

        else
        {
            fprintf(stderr, "unexpected end of file inside TORSDIHEDRALH block\n");
            exit(EXIT_FAILURE);
        }

        if (lines_read == 1)
        {
            if (sscanf(line_buffer, "%d", &n_dihedralh) != 1)
            {
                fprintf(stderr, "error reading number of dihedrals in TORSDIHEDRALH block\n");
                exit(EXIT_FAILURE);
            }
            if (n_dihedralh == 0)
                return;
        }

        if (lines_read > 1)
        {
            topo_data->n_dihedrals++;
            int idx = topo_data->n_dihedrals- 1;

            if (sscanf(line_buffer, "%d%d%d%d%d", &topo_data->dihedrals[idx].ai, &topo_data->dihedrals[idx].aj, &topo_data->dihedrals[idx].ak, &topo_data->dihedrals[idx].al, &topo_data->dihedrals[idx].code) != 5)
            {
                fprintf(stderr, "error: reading dihedral in block TORSDIHEDRALH. unexpected format in line \n%s", line_buffer);
                exit(EXIT_FAILURE);
            }

            if (lines_read == n_dihedralh + 1)
            {
                /* finished parsing block */
                return;
            }
        }
    }
}

void parse_block_dihedral (FILE* fp, topology_t* topo_data, forcefield_t* ff_data)
{
    int lines_read = 0;
    int n_dihedral = -1;
    char line_buffer[TOP_LINE_BUFFER];

    while (1)
    {
        if (top_read_next_line (fp, line_buffer) != -1)
            lines_read++;

        else
        {
            fprintf(stderr, "unexpected end of file inside TORSDIHEDRAL block\n");
            exit(EXIT_FAILURE);
        }

        if (lines_read == 1)
        {
            if (sscanf(line_buffer, "%d", &n_dihedral) != 1)
            {
                fprintf(stderr, "error reading number of dihedrals in TORSDIHEDRAL block\n");
                exit(EXIT_FAILURE);
            }
            if (n_dihedral == 0)
                return;
        }

        if (lines_read > 1)
        {
            topo_data->n_dihedrals++;
            int idx = topo_data->n_dihedrals - 1;

            if (sscanf(line_buffer, "%d%d%d%d%d", &topo_data->dihedrals[idx].ai, &topo_data->dihedrals[idx].aj, &topo_data->dihedrals[idx].ak, &topo_data->dihedrals[idx].al, &topo_data->dihedrals[idx].code) != 5)
            {
                fprintf(stderr, "error: reading dihedral in block TORSDIHEDRAL. unexpected format in line \n%s", line_buffer);
                exit(EXIT_FAILURE);
            }

            if (lines_read == n_dihedral + 1)
            {
                /* finished parsing block */
                return;
            }
        }
    }
}

/* master function for reading the GROMOS topology file */
void read_gromos_topology (char* fn, struct topology_t* topo_data, struct forcefield_t* ff_data)
{
    FILE* fp = fopen (fn, "r");
    char line_buffer[TOP_LINE_BUFFER];


    /* reset bonds, angles, etc. */
    topo_data->n_bonds = 0;
    topo_data->n_angles = 0;
    topo_data->n_dihedrals = 0;
    topo_data->n_impropers = 0;
    topo_data->n_pairs = 0;

    while (1)
    {
        /* if no line can be read, break */
        if (top_read_next_line (fp, line_buffer) == -1)
            break;

        /* if line is read, check for a particular block start */
        /* fork program according to results */

        if (is_block (line_buffer, "ATOMTYPENAME"))
            parse_block_atomtypename (fp, topo_data, ff_data);

        if (is_block (line_buffer, "RESNAME"))
            parse_block_resname (fp, topo_data, ff_data);

        if (is_block (line_buffer, "SOLUTEATOM"))
            parse_block_soluteatom (fp, topo_data, ff_data);

        if (is_block (line_buffer, "LJPARAMETERS"))
            parse_block_ljparameters (fp, topo_data, ff_data);

        if (is_block (line_buffer, "BONDSTRETCHTYPE"))
            parse_block_bondstretchtype (fp, topo_data, ff_data);
        else
        {
            if (is_block (line_buffer, "BONDH"))
                parse_block_bondh (fp, topo_data, ff_data);
            else
            {
                /* dummy block -- avoid entering BONDANGLE*** blocks */
                if ( (is_block (line_buffer, "BOND")) && !(is_block(line_buffer, "BONDA")))
                    parse_block_bond (fp, topo_data, ff_data);
            }
        }

        if (is_block (line_buffer, "BONDANGLEBENDTYPE"))
            parse_block_bondanglebendtype (fp, topo_data, ff_data);
        else
        {
            if (is_block (line_buffer, "BONDANGLEH"))
                parse_block_bondangleh (fp, topo_data, ff_data);
            else
            {
                if (is_block (line_buffer, "BONDANGLE"))
                    parse_block_bondangle (fp, topo_data, ff_data);
            }
        }

        if (is_block (line_buffer, "IMPDIHEDRALTYPE"))
            parse_block_impdihedraltype (fp, topo_data, ff_data);
        else
        {
            if (is_block (line_buffer, "IMPDIHEDRALH"))
                parse_block_impdihedralh (fp, topo_data, ff_data);
            else
            {
                if (is_block (line_buffer, "IMPDIHEDRAL"))
                    parse_block_impdihedral (fp, topo_data, ff_data);
            }
        }

        if (is_block (line_buffer, "TORSDIHEDRALTYPE"))
            parse_block_dihedraltype (fp, topo_data, ff_data);
        else
        {
            if (is_block (line_buffer, "DIHEDRALH"))
                parse_block_dihedralh (fp, topo_data, ff_data);
            else
            {
                if (is_block (line_buffer, "DIHEDRAL"))
                    parse_block_dihedral (fp, topo_data, ff_data);
            }
        }

    }

    fclose(fp);
    return;
}
#endif
