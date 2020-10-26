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

#ifndef _REGEX_MATCHER_
#define _REGEX_MATCHER_

#define REGEX_STRING_MAX 20

#include <regex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int regex_match (char* string, char* regex_pattern)
{
    regex_t regex;
    int comp_regex;
    
    comp_regex = regcomp(&regex, regex_pattern, 0);

    if (comp_regex) {
        fprintf(stderr, "Regex pattern %s looks wrong. Re-check it.\n", regex_pattern);
        exit(1);
    }

    comp_regex = regexec(&regex, string, 0, NULL, 0);

    regfree(&regex);
    return !comp_regex;
}

int compare_pair_with_pattern (char* ti, char* tj, char* pattern)
{
    char dord[REGEX_STRING_MAX];
    char rord[REGEX_STRING_MAX];

    sprintf(dord, "%s-%s", ti, tj);
    sprintf(rord, "%s-%s", tj, ti);

    if ( regex_match(dord, pattern) || regex_match(rord, pattern) )
        return 1;
    else
        return 0;
}

int compare_dih_with_pattern (char* ti, char* tj, char* tk, char* tl, char* pattern)
{
    char dord[REGEX_STRING_MAX];
    char rord[REGEX_STRING_MAX];

    sprintf(dord, "%s-%s-%s-%s", ti, tj, tk, tl);
    sprintf(rord, "%s-%s-%s-%s", tl, tk, tj, ti);

    if ( regex_match(dord, pattern) || regex_match(rord, pattern) )
        return 1;
    else
        return 0;
}

#endif
