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

#ifndef _ARGPARSE_
#define _ARGPARSE_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MULTIOPT_OPTION_BUFFER 20

int multiopt_check_unknown (int argc,
        char** argv,
        int nflags,
        const char known_flags[][MULTIOPT_OPTION_BUFFER])
{
    for (int i = 0; i < argc; ++i)
    {
        if (argv[i][0] == '-')
        {
            int match = 0;
            for (int j = 0; j < nflags; j++)
            {
                if (!strcmp(argv[i],known_flags[j]))
                    match = 1;
            }
            if (!match)
            {
                fprintf(stderr, "error: flag %s is not known\n", argv[i]);
                exit(EXIT_FAILURE);
            }
        }
    }
    return 1;
}

/* Verifies if a flag is there. */
int multiopt_bool_flag (int argc,
       char** argv,
      const char* flag)
{
    for (int i = 0; i < argc; i++)
    {
        if (!strcmp(flag,argv[i]))
        {
            return 1;
        }
    }
    return 0;
} 

/* Get number of options following a flag. */
int multiopt_getn (int argc,
        char** argv,
        const char* flag)
{
    int n = -1;
    for (int i = 0; i < argc; i++)
    {
        if (!strcmp(flag,argv[i]))
        {
            n = 0;
        }
        else if (n > -1 && argv[i][0] == '-')
        {
            /* This identifies that a flag is coming. */
            break;
        }
        else if (n > -1)
        {
            n++;
        }
    }
    /* This will happen if the flag is not found. */
    if (n == -1)
    {
        fprintf(stderr, "warning: flag %s was not found \n", flag);
        exit(EXIT_FAILURE);
    }
    return n;
}

char* multiopt2str (int argc,
                 char** argv,
            const char* flag,
            int        entry)
{
    char *out;
    for (int i = 0; i < argc; ++i)
    {
        if (!strcmp(flag,argv[i]))
        {
            int n = strlen(argv[i+entry+1]);
            out = (char*) malloc ((n+1)*sizeof(char));
            strcpy(out, argv[i+1+entry]);
            return out;
        }
    }
    return NULL;
}

double multiopt2double (int argc,
                 char** argv,
            const char* flag,
            int        entry)
{
    char *str;
    double db;
    str = multiopt2str(argc,argv,flag,entry);
    db = atof(str);
    free(str);
    return db;
}

double multiopt2int (int argc,
                 char** argv,
            const char* flag,
            int        entry)
{
    char *str;
    int db;
    str = multiopt2str(argc,argv,flag,entry);
    db = atoi(str);
    free(str);
    return db;
}

#endif
