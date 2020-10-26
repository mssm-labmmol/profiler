#ifndef _FORCEFIELD_
#define _FORCEFIELD_

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define  FF_ATOMTYPE_BUFFER       50
#define  FF_N_ATOMTYPES_MAX_SIZE   500
#define  FF_N_PAIRTYPES_MAX_SIZE   125000
#define  FF_N_BONDS_MAX_SIZE   5000
#define  FF_N_ANGLES_MAX_SIZE  5000
#define  FF_N_DIH_MAX_SIZE     5000
#define  FF_N_IDIH_MAX_SIZE    5000

struct atomtype_t {
    int code;
    char type[FF_ATOMTYPE_BUFFER];
};

struct pairtype_pars_t {
    int iac, jac;
    double c6, c12, cs6, cs12;
    /* NOTE: Charges are not needed because they can be retrieved from IAC and
     * JAC in the context of the forcefield_t structure. */
};

struct bond_pars_t {
    int code;
    double l,k;
};

struct angle_pars_t {
    int code;
    double theta, k;
};

struct dih_pars_t {
    int code;
    int mult;
    double phi,k;
};

struct idih_pars_t {
    int code;
    double phi,k;
};

struct forcefield_t {

    int n_atomtypes, n_pairtypes, n_bondtypes, n_angletypes, n_dihtypes, n_idihtypes;
    struct atomtype_t atomtypes[FF_N_ATOMTYPES_MAX_SIZE];
    struct pairtype_pars_t pairtypes[FF_N_PAIRTYPES_MAX_SIZE];
    struct bond_pars_t bondtypes[FF_N_BONDS_MAX_SIZE];
    struct angle_pars_t angletypes[FF_N_ANGLES_MAX_SIZE];
    struct dih_pars_t dihtypes[FF_N_DIH_MAX_SIZE];
    struct idih_pars_t idihtypes[FF_N_IDIH_MAX_SIZE];
};

void print_atomtypes (struct forcefield_t* ff_data)
{
    for (int i = 0; i < ff_data->n_atomtypes; i++)
    {
        fprintf(stderr, "atomtype iac %d name %s\n", ff_data->atomtypes[i].code, ff_data->atomtypes[i].type);
    }
}

void get_type_for_iac (struct forcefield_t* ff_data, int iac, char* out_type)
{
    for (int i = 0; i < ff_data->n_atomtypes; i++)
    {
        if (ff_data->atomtypes[i].code == iac)
        {
            strcpy(out_type, ff_data->atomtypes[i].type);
            return;
        }
    }

    fprintf(stderr, "error: atomtype with IAC %d was not found.\n", iac);
    exit(EXIT_FAILURE);
}

void get_pairtype_parameters_for_code (struct forcefield_t* ff_data, int iac, int jac, struct pairtype_pars_t* out_pairtype)
{
    for (int i = 0; i < ff_data->n_pairtypes; i++)
    {
        if ( ((ff_data->pairtypes[i].iac == iac) && (ff_data->pairtypes[i].jac == jac)) || 
                ((ff_data->pairtypes[i].iac == jac) && (ff_data->pairtypes[i].jac == iac)) )
        {
            out_pairtype->iac = ff_data->pairtypes[i].iac;
            out_pairtype->jac = ff_data->pairtypes[i].jac;
            out_pairtype->c6 = ff_data->pairtypes[i].c6;
            out_pairtype->c12 = ff_data->pairtypes[i].c12;
            out_pairtype->cs6 = ff_data->pairtypes[i].cs6;
            out_pairtype->cs12 = ff_data->pairtypes[i].cs12;
            return;
        }
    }
    fprintf(stderr, "error: pair type with codes %d,%d was not found.\n", iac, jac);
    exit(EXIT_FAILURE);
}

void get_bond_parameters_for_code (struct forcefield_t* ff_data, int code, struct bond_pars_t* out_bond)
{
    for (int i = 0; i < ff_data->n_bondtypes; i++)
    {
        if (ff_data->bondtypes[i].code == code)
        {
            out_bond->k = ff_data->bondtypes[i].k;
            out_bond->l = ff_data->bondtypes[i].l;
            return;
        }
    }

    fprintf(stderr, "error: bond type with code %d was not found.\n", code);
    exit(EXIT_FAILURE);
}

void get_angle_parameters_for_code (struct forcefield_t* ff_data, int code, struct angle_pars_t* out_angle)
{
    for (int i = 0; i < ff_data->n_angletypes; i++)
    {
        if (ff_data->angletypes[i].code == code)
        {
            out_angle->k = ff_data->angletypes[i].k;
            out_angle->theta = ff_data->angletypes[i].theta;
            return;
        }
    }

    fprintf(stderr, "error: angle type with code %d was not found.\n", code);
    exit(EXIT_FAILURE);
}

void get_dihedral_parameters_for_code (struct forcefield_t* ff_data, int code, struct dih_pars_t* out_dih)
{
    for (int i = 0; i < ff_data->n_dihtypes; i++)
    {
        if (ff_data->dihtypes[i].code == code)
        {
            out_dih->k = ff_data->dihtypes[i].k;
            out_dih->phi = ff_data->dihtypes[i].phi;
            out_dih->mult = ff_data->dihtypes[i].mult;
            return;
        }
    }

    fprintf(stderr, "error: dih type with code %d was not found.\n", code);
    exit(EXIT_FAILURE);
}

void get_improper_parameters_for_code (struct forcefield_t* ff_data, int code, struct idih_pars_t* out_idih)
{
    for (int i = 0; i < ff_data->n_idihtypes; i++)
    {
        if (ff_data->idihtypes[i].code == code)
        {
            out_idih->k = ff_data->idihtypes[i].k;
            out_idih->phi = ff_data->idihtypes[i].phi;
            return;
        }
    }

    fprintf(stderr, "error: idih type with code %d was not found.\n", code);
    exit(EXIT_FAILURE);
}
#endif
