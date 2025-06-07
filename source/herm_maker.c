/* herm_maker.c – Hermitian Matrix Maker external for Max/MSP
 *
 * This external takes an input list representing an upper–triangular matrix
 * (n×n, provided in row–major order; the entries below the diagonal are assumed
 * to be zero and the diagonal entries are real) and outputs the full Hermitian matrix.
 *
 * The computation is:
 *   For i == j:
 *      H[i][j] = 2 * U[i][j]
 *   For i != j:
 *      H[i][j] = U[i][j] + U[j][i]
 *
 * The output is sent as a flat list.
 *
 * Compile with Xcode using the Max SDK.
 */

#include "ext.h"
#include "ext_obex.h"
#include <stdlib.h>
#include <stdio.h>

typedef struct _herm_maker {
    t_object ob;
    long n;        // Matrix dimension.
    void *out;     // Outlet pointer.
} t_herm_maker;

/* Global class pointer */
static t_class *herm_maker_class = NULL;

/* Function prototypes */
void *herm_maker_new(t_symbol *s, long argc, t_atom *argv);
void herm_maker_free(t_herm_maker *x);
void herm_maker_assist(t_herm_maker *x, void *b, long m, long a, char *s);
void herm_maker_list(t_herm_maker *x, t_symbol *s, long argc, t_atom *argv);

/* Main entry point, called by Max at load time */
void ext_main(void *r)
{
    t_class *c = class_new("herm_maker", (method)herm_maker_new, (method)herm_maker_free,
                             sizeof(t_herm_maker), 0L, A_GIMME, 0);
    
    class_addmethod(c, (method)herm_maker_list, "list", A_GIMME, 0);
    class_addmethod(c, (method)herm_maker_assist, "assist", A_CANT, 0);
    
    class_register(CLASS_BOX, c);
    herm_maker_class = c;  // Save the class pointer for use in new
}

/* Create a new instance of the object */
void *herm_maker_new(t_symbol *s, long argc, t_atom *argv)
{
    t_herm_maker *x = (t_herm_maker *)object_alloc(herm_maker_class);
    if (x) {
        // Set default matrix dimension; default is 3.
        x->n = 3;
        if (argc > 0) {
            if (atom_gettype(argv) == A_LONG)
                x->n = atom_getlong(argv);
            else if (atom_gettype(argv) == A_FLOAT)
                x->n = (long)atom_getfloat(argv);
        }
        if (x->n <= 0)
            x->n = 3;
        // Create an outlet.
        x->out = outlet_new(x, NULL);
    }
    return x;
}

/* Free the object (nothing to free in this simple example) */
void herm_maker_free(t_herm_maker *x)
{
    ; // No dynamic memory allocated in the object structure
}

/* Provide assistance messages for inlets/outlets */
void herm_maker_assist(t_herm_maker *x, void *b, long m, long a, char *s)
{
    if (m == 1)
        sprintf(s, "Input: List representing an upper-triangular matrix (%ld numbers)", x->n * x->n);
    else
        sprintf(s, "Output: Hermitian matrix as flat list (%ld numbers)", x->n * x->n);
}

/* Process the input list message and output the Hermitian matrix */
void herm_maker_list(t_herm_maker *x, t_symbol *s, long argc, t_atom *argv)
{
    long n = x->n;
    long total = n * n;
    if (argc != total) {
        object_post((t_object *)x, "Expected %ld numbers, received %ld", total, argc);
        return;
    }
    
    // Allocate an output array for the Hermitian matrix.
    t_atom *out_list = (t_atom *)sysmem_newptr(total * sizeof(t_atom));
    if (!out_list) {
        object_error((t_object *)x, "Memory allocation failed");
        return;
    }
    
    // Compute the full Hermitian matrix:
    // For diagonal entries: H[i][i] = 2 * U[i][i].
    // For off-diagonals: H[i][j] = U[i][j] + U[j][i].
    for (long i = 0; i < n; i++) {
        for (long j = 0; j < n; j++) {
            long index = i * n + j;
            double val;
            if (i == j) {
                val = 2.0 * atom_getfloat(argv + index);
            } else {
                double v1 = atom_getfloat(argv + index);
                double v2 = atom_getfloat(argv + (j * n + i));
                val = v1 + v2;
            }
            atom_setfloat(out_list + index, val);
        }
    }
    
    // Output the resulting Hermitian matrix.
    outlet_list(x->out, gensym("list"), total, out_list);
    sysmem_freeptr(out_list);
}
 * * This external takes an input list representing an upper–triangular matrix * (n×n, provided in row–major order; the entries below the diagonal are assumed * to be zero and the diagonal entries are real) and outputs the full Hermitian matrix. * * The computation is: *   For i == j: *      H[i][j] = 2 * U[i][j] *   For i != j: *      H[i][j] = U[i][j] + U[j][i] * * The output is sent as a flat list. * * Compile with Xcode using the Max SDK. */#include "ext.h"#include "ext_obex.h"#include <stdlib.h>#include <stdio.h>typedef struct _herm_maker {    t_object ob;    long n;        // Matrix dimension.    void *out;     // Outlet pointer.} t_herm_maker;/* Global class pointer */static t_class *herm_maker_class = NULL;/* Function prototypes */void *herm_maker_new(t_symbol *s, long argc, t_atom *argv);void herm_maker_free(t_herm_maker *x);void herm_maker_assist(t_herm_maker *x, void *b, long m, long a, char *s);void herm_maker_list(t_herm_maker *x, t_symbol *s, long argc, t_atom *argv);/* Main entry point, called by Max at load time */void ext_main(void *r){    t_class *c = class_new("herm_maker", (method)herm_maker_new, (method)herm_maker_free,                             sizeof(t_herm_maker), 0L, A_GIMME, 0);        class_addmethod(c, (method)herm_maker_list, "list", A_GIMME, 0);    class_addmethod(c, (method)herm_maker_assist, "assist", A_CANT, 0);        class_register(CLASS_BOX, c);    herm_maker_class = c;  // Save the class pointer for use in new}/* Create a new instance of the object */void *herm_maker_new(t_symbol *s, long argc, t_atom *argv){    t_herm_maker *x = (t_herm_maker *)object_alloc(herm_maker_class);    if (x) {        // Set default matrix dimension; default is 3.        x->n = 3;        if (argc > 0) {            if (atom_gettype(argv) == A_LONG)                x->n = atom_getlong(argv);            else if (atom_gettype(argv) == A_FLOAT)                x->n = (long)atom_getfloat(argv);        }        if (x->n <= 0)            x->n = 3;        // Create an outlet.        x->out = outlet_new(x, NULL);    }    return x;}/* Free the object (nothing to free in this simple example) */void herm_maker_free(t_herm_maker *x){    ; // No dynamic memory allocated in the object structure}/* Provide assistance messages for inlets/outlets */void herm_maker_assist(t_herm_maker *x, void *b, long m, long a, char *s){    if (m == 1)        sprintf(s, "Input: List representing an upper-triangular matrix (%ld numbers)", x->n * x->n);    else        sprintf(s, "Output: Hermitian matrix as flat list (%ld numbers)", x->n * x->n);}/* Process the input list message and output the Hermitian matrix */void herm_maker_list(t_herm_maker *x, t_symbol *s, long argc, t_atom *argv){    long n = x->n;    long total = n * n;    if (argc != total) {        object_post((t_object *)x, "Expected %ld numbers, received %ld", total, argc);        return;    }        // Allocate an output array for the Hermitian matrix.    t_atom *out_list = (t_atom *)sysmem_newptr(total * sizeof(t_atom));    if (!out_list) {        object_error((t_object *)x, "Memory allocation failed");        return;    }        // Compute the full Hermitian matrix:    // For diagonal entries: H[i][i] = 2 * U[i][i].    // For off-diagonals: H[i][j] = U[i][j] + U[j][i].    for (long i = 0; i < n; i++) {        for (long j = 0; j < n; j++) {            long index = i * n + j;            double val;            if (i == j) {                val = 2.0 * atom_getfloat(argv + index);            } else {                double v1 = atom_getfloat(argv + index);                double v2 = atom_getfloat(argv + (j * n + i));                val = v1 + v2;            }            atom_setfloat(out_list + index, val);        }    }        // Output the resulting Hermitian matrix.    outlet_list(x->out, gensym("list"), total, out_list);    sysmem_freeptr(out_list);}