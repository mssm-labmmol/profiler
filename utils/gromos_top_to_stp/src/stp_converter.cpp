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

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>

#include "topo.h"
#include "forcefield.h"
#include "gromos_top_parser.h"
#include "stp_io.h"
#include "argparse.h"

#define MAX_TOPS 500

int main (int argc, char** argv)
{
    char* topos[MAX_TOPS];
    char* pairs_pattern, *dihedrals_pattern;
    int  ntopos;
    int  sel_atoms[N_ATOMS_MAX_SIZE];
    int  sel_pairs[N_PAIRS_MAX_SIZE];
    int  nsel = 0;
    int  nsel_pairs = 0;
    int  sel_dih[N_DIH_MAX_SIZE];
    int  nsel_dih;
    int  refdih;
    int  optflag = 0;

    /* flags for input options */
    int flag_tops  = multiopt_bool_flag (argc, argv, "-t");
    int flag_atoms = multiopt_bool_flag (argc, argv, "-a");
    int flag_dihs  = multiopt_bool_flag (argc, argv, "-d");
    int flag_ref   = multiopt_bool_flag (argc, argv, "-r");
    int flag_pairs = multiopt_bool_flag (argc, argv, "-p");
    char known_flags[][MULTIOPT_OPTION_BUFFER] = {"-t", "-a", "-p", "-d", "-r", "-inter"};

    multiopt_check_unknown (argc, argv, 6, known_flags);

    /* consistency check */
    if ( !flag_tops )
    {
        fprintf(stderr, "error: option -t is required\n");
        exit(EXIT_FAILURE);
    }
    if (flag_atoms && multiopt_getn(argc, argv, "-a") != 1)
    {
        fprintf(stderr, "error: -a requires exactly one type\n");
        exit(EXIT_FAILURE);
    }
    if (flag_dihs && multiopt_getn(argc, argv, "-d") != 1)
    {
        fprintf(stderr, "error: -d requires exactly one pattern\n");
        exit(EXIT_FAILURE);
    }
    if (flag_pairs && multiopt_getn(argc, argv, "-p") != 1)
    {
        fprintf(stderr, "error: -p requires exactly one pattern\n");
        exit(EXIT_FAILURE);
    }

    if (multiopt_bool_flag(argc, argv, "-a") ||
	multiopt_bool_flag(argc, argv, "-d") ||
	multiopt_bool_flag(argc, argv, "-r") ||
	multiopt_bool_flag(argc, argv, "-p") ||
	multiopt_bool_flag(argc, argv, "-inter")) {
	optflag = 1;
    } 

    /* get values of options */
    ntopos = multiopt_getn(argc, argv, "-t");
    if (ntopos > MAX_TOPS)
    {
        fprintf(stderr, "error: Number of topologies is larger than maximum allowed of %d.\n", MAX_TOPS);
        fprintf(stderr, "       Recompile the program defining a larger value for the MAX_TOPS macro.\n");
        exit(EXIT_FAILURE);
    }
    else if (ntopos < 1)
    {
        fprintf(stderr, "You have to supply at least one topology.\n");
        exit(EXIT_FAILURE);
    }
    else
    {
        fprintf(stderr, "Read %d topologies.\n", ntopos);
    }

    for (int i = 0; i < ntopos; i++)
    {
        topos[i] = multiopt2str(argc,argv,"-t",i);

        int namesize = strlen(topos[i]);

        if (strcmp(&(topos[i][namesize-4]), ".top"))
        {
            fprintf(stderr, "error: topology extension must be *.top\n");
            exit(EXIT_FAILURE);
        }
        else
        {
            fprintf(stdout, "topology %d: %s\n", i+1, topos[i]);
        }
    }

    if ( (flag_atoms && flag_pairs) ) {
      fprintf(stderr, "error: cannot have -a <pattern> and -p <pattern> at the same time\n");
      exit(EXIT_FAILURE);
    }
   
    if (flag_dihs)
        dihedrals_pattern = multiopt2str(argc, argv, "-d", 0); 
    else
        dihedrals_pattern = "NULL";

    if (flag_pairs)
      pairs_pattern = multiopt2str(argc, argv, "-p", 0);
    else
      pairs_pattern = "NULL";

    for (int i = 0; i < ntopos; i++)
    {
        struct topology_t topo_data;
        struct forcefield_t ff_data;
        char filename[512];
        FILE *fp;

        if (strrchr(topos[i], '/') != NULL)
            sprintf(filename, "%s", (char *) strrchr(topos[i],'/') + 1);
        else
            sprintf(filename, "%s", topos[i]);

        int n = strlen(filename);
       
        /* change extension to *.stp */
        filename[n-3] = 's';
        filename[n-2] = 't';

        fprintf(stderr, "Writing to %s \n", filename);
        fp = fopen(filename, "w");

        /* reads topology filling topo_data with topology data and ff_data with force-field data */
        read_gromos_topology (topos[i], &topo_data, &ff_data);

        /* select atoms */
        if (flag_atoms)
        {
            nsel = select_atoms_by_type (multiopt2str(argc, argv, "-a", 0), &topo_data, &ff_data,
                    sel_atoms);
        }
        else
        {
	  if (optflag) {
	    if (flag_pairs) {
	      nsel_pairs = select_pairs_by_pattern(pairs_pattern, &topo_data, &ff_data, sel_pairs);
	    } else {
	      nsel = select_atoms_by_inter (&topo_data, &ff_data, sel_atoms);
	      if (nsel == -1) {
		nsel_pairs = select_pairs_by_inter (&topo_data, &ff_data, sel_pairs);
	      }
	    }
	  }
        }

        /* select dihedrals */
        if (flag_dihs)
        {
            nsel_dih = select_dihedrals_by_pattern(dihedrals_pattern, &topo_data, &ff_data, sel_dih);
        }
        else
        {
	  if (optflag) {
            nsel_dih = select_dihedrals_by_inter(&topo_data, &ff_data, sel_dih);
	  }
        }

        /* ref */
        if (flag_ref)
        {
            refdih = multiopt2int(argc, argv, "-r", i);
        }
        else
        {
	  if (optflag) {
            refdih = select_ref_by_inter(&topo_data, &ff_data);
	  }
        }

        fprintf(stderr, "Finished writing to %s \n", filename);
	if (nsel_pairs == 0) {
	  write_stp_to_stream (fp, &topo_data, &ff_data, sel_atoms, nsel, sel_dih, nsel_dih, refdih, optflag);
	} else {
	  write_pairs_stp_to_stream (fp, &topo_data, &ff_data, sel_pairs, nsel_pairs, sel_dih, nsel_dih, refdih, optflag);
	}

        fclose(fp);
    }

    return 0;
}
