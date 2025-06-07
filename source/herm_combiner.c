#include "ext.h"
#include "ext_obex.h"
#include <stdlib.h>

typedef struct _herm_comb {
    t_object ob;
    long n;    // matrix dimension
    void *out; // outlet pointer
} t_herm_comb;

/* Global class pointer */
static t_class *herm_comb_class = NULL;
void *herm_comb_new(t_symbol *s, long argc, t_atom *argv);
void herm_comb_free(t_herm_comb *x);
void herm_comb_assist(t_herm_comb *x, void *b, long m, long a, char *s);
void herm_comb_list(t_herm_comb *x, t_symbol *s, long argc, t_atom *argv);

void ext_main(void *r) {
    t_class *c;
    c = class_new("herm_combiner", (method)herm_comb_new, (method)herm_comb_free, sizeof(t_herm_comb), 0L, A_GIMME, 0);
    class_addmethod(c, (method)herm_comb_list, "list", A_GIMME, 0);
    class_addmethod(c, (method)herm_comb_assist, "assist", A_CANT, 0);
    class_register(CLASS_BOX, c);
    herm_comb_class = c;
}

void *herm_comb_new(t_symbol *s, long argc, t_atom *argv) {
    t_herm_comb *x = (t_herm_comb *)object_alloc(herm_comb_class);
    if (x) {
        x->n = 3; // default dimension
        if(argc > 0 && atom_gettype(argv) == A_LONG)
            x->n = atom_getlong(argv);
        x->out = outlet_new(x, NULL);
    }
    return (x);
}

void herm_comb_free(t_herm_comb *x) { }

void herm_comb_assist(t_herm_comb *x, void *b, long m, long a, char *s) {
    if(m == 1)
        sprintf(s, "Input: List of numbers for matrices & coefficients");
    else
        sprintf(s, "Output: Combined Hermitian matrix as list");
}

void herm_comb_list(t_herm_comb *x, t_symbol *s, long argc, t_atom *argv) {
    // For demonstration, we assume two matrices are concatenated, each with n*n numbers,
    // followed by two coefficients.
    long n = x->n;
    if(argc != 2*n*n + 2) {
        object_post((t_object *)x, "Expected %ld numbers", 2*n*n+2);
        return;
    }
    double *M1 = (double *)sysmem_newptr(n*n * sizeof(double));
    double *M2 = (double *)sysmem_newptr(n*n * sizeof(double));
    double c1 = atom_getfloat(argv + 2*n*n);
    double c2 = atom_getfloat(argv + 2*n*n + 1);
    for(long i = 0; i < n*n; i++){
        M1[i] = atom_getfloat(argv + i);
        M2[i] = atom_getfloat(argv + n*n + i);
    }
    t_atom *out_list = (t_atom *)sysmem_newptr(n*n * sizeof(t_atom));
    for(long i = 0; i < n*n; i++){
        double val = c1 * M1[i] + c2 * M2[i];
        atom_setfloat(out_list + i, val);
    }
    outlet_list(x->out, gensym("list"), n*n, out_list);
    sysmem_freeptr(out_list);
    sysmem_freeptr(M1);
    sysmem_freeptr(M2);
}

#include "ext_obex.h"#include <stdlib.h>typedef struct _herm_comb {    t_object ob;    long n;    // matrix dimension    void *out; // outlet pointer} t_herm_comb;void *herm_comb_new(t_symbol *s, long argc, t_atom *argv);void herm_comb_free(t_herm_comb *x);void herm_comb_assist(t_herm_comb *x, void *b, long m, long a, char *s);void herm_comb_list(t_herm_comb *x, t_symbol *s, long argc, t_atom *argv);void ext_main(void *r) {    t_class *c;    c = class_new("herm_combiner", (method)herm_comb_new, (method)herm_comb_free, sizeof(t_herm_comb), 0L, A_GIMME, 0);    class_addmethod(c, (method)herm_comb_list, "list", A_GIMME, 0);    class_addmethod(c, (method)herm_comb_assist, "assist", A_CANT, 0);    class_register(CLASS_BOX, c);}void *herm_comb_new(t_symbol *s, long argc, t_atom *argv) {    t_herm_comb *x = (t_herm_comb *)object_alloc(/* class pointer */);    if (x) {        x->n = 3; // default dimension        if(argc > 0 && atom_gettype(argv) == A_LONG)            x->n = atom_getlong(argv);        x->out = outlet_new(x, NULL);    }    return (x);}void herm_comb_free(t_herm_comb *x) { }void herm_comb_assist(t_herm_comb *x, void *b, long m, long a, char *s) {    if(m == 1)        sprintf(s, "Input: List of numbers for matrices & coefficients");    else        sprintf(s, "Output: Combined Hermitian matrix as list");}void herm_comb_list(t_herm_comb *x, t_symbol *s, long argc, t_atom *argv) {    // For demonstration, we assume two matrices are concatenated, each with n*n numbers,    // followed by two coefficients.    long n = x->n;    if(argc != 2*n*n + 2) {        object_post((t_object *)x, "Expected %ld numbers", 2*n*n+2);        return;    }    double *M1 = (double *)sysmem_newptr(n*n * sizeof(double));    double *M2 = (double *)sysmem_newptr(n*n * sizeof(double));    double c1 = atom_getfloat(argv + 2*n*n);    double c2 = atom_getfloat(argv + 2*n*n + 1);    for(long i = 0; i < n*n; i++){        M1[i] = atom_getfloat(argv + i);        M2[i] = atom_getfloat(argv + n*n + i);    }    t_atom *out_list = (t_atom *)sysmem_newptr(n*n * sizeof(t_atom));    for(long i = 0; i < n*n; i++){        double val = c1 * M1[i] + c2 * M2[i];        atom_setfloat(out_list + i, val);    }    outlet_list(x->out, gensym("list"), n*n, out_list);    sysmem_freeptr(out_list);    sysmem_freeptr(M1);    sysmem_freeptr(M2);}