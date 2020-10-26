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
