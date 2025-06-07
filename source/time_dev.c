#include "ext.h"
#include "ext_obex.h"
#include <stdlib.h>
#include <math.h>
#include <complex.h>

typedef struct _time_dev {
    t_object ob;
    long n;      // number of eigenstates
    long tsteps; // number of time steps
    double tmin; // minimum time
    double tmax; // maximum time
    void *out;   // outlet pointer
} t_time_dev;

/* Global class pointer */
static t_class *time_dev_class = NULL;
void *time_dev_new(t_symbol *s, long argc, t_atom *argv);
void time_dev_free(t_time_dev *x);
void time_dev_assist(t_time_dev *x, void *b, long m, long a, char *s);
void time_dev_bang(t_time_dev *x);

void ext_main(void *r) {
    t_class *c;
    c = class_new("time_dev", (method)time_dev_new, (method)time_dev_free, sizeof(t_time_dev), 0L, A_GIMME, 0);
    class_addmethod(c, (method)time_dev_bang, "bang", 0);
    class_addmethod(c, (method)time_dev_assist, "assist", A_CANT, 0);
    class_register(CLASS_BOX, c);
    time_dev_class = c;
}

void *time_dev_new(t_symbol *s, long argc, t_atom *argv) {
    t_time_dev *x = (t_time_dev *)object_alloc(time_dev_class);
    if(x){
        x->n = 3;        // default: 3 eigenstates
        x->tsteps = 100; // default time resolution
        x->tmin = 0.0;
        x->tmax = 1.0;
        if(argc >= 1)
            x->n = atom_getlong(argv);
        if(argc >= 2)
            x->tsteps = atom_getlong(argv+1);
        if(argc >= 3)
            x->tmin = atom_getfloat(argv+2);
        if(argc >= 4)
            x->tmax = atom_getfloat(argv+3);
        x->out = outlet_new(x, NULL);
    }
    return (x);
}

void time_dev_free(t_time_dev *x) { }

void time_dev_assist(t_time_dev *x, void *b, long m, long a, char *s) {
    if(m==1)
        sprintf(s, "Bang to compute time development");
    else
        sprintf(s, "Output: List of magnitude & phase pairs per time step");
}

void time_dev_bang(t_time_dev *x) {
    long n = x->n;
    long tsteps = x->tsteps;
    double tmin = x->tmin;
    double tmax = x->tmax;
    double dt = (tmax - tmin) / (tsteps - 1);
    
    // For demonstration, use dummy eigenvalues and coefficients.
    double *eigenvalues = (double *)sysmem_newptr(n * sizeof(double));
    double complex *coeff = (double complex *)sysmem_newptr(n * sizeof(double complex));
    for(long k = 0; k < n; k++){
        eigenvalues[k] = 1.0 + k;      // dummy eigenvalues
        coeff[k] = 1.0 + 0.0*I;        // dummy coefficients
    }
    // For each time step, compute: Ψ(t) = Σ_k coeff[k] * exp(-i * eigenvalues[k] * t)
    t_atom *out_list = (t_atom *)sysmem_newptr(tsteps * 2 * sizeof(t_atom)); // each time step: mag and phase
    for(long t = 0; t < tsteps; t++){
        double time = tmin + t * dt;
        double complex psi = 0.0 + 0.0*I;
        for(long k = 0; k < n; k++){
            psi += coeff[k] * cexp(-I * eigenvalues[k] * time);
        }
        double mag = cabs(psi);
        double phase = carg(psi);
        atom_setfloat(out_list + t*2, mag);
        atom_setfloat(out_list + t*2 + 1, phase);
    }
    outlet_list(x->out, gensym("list"), tsteps*2, out_list);
    sysmem_freeptr(out_list);
    sysmem_freeptr(eigenvalues);
    sysmem_freeptr(coeff);
}

#include "ext_obex.h"#include <stdlib.h>#include <math.h>#include <complex.h>typedef struct _time_dev {    t_object ob;    long n;      // number of eigenstates    long tsteps; // number of time steps    double tmin; // minimum time    double tmax; // maximum time    void *out;   // outlet pointer} t_time_dev;void *time_dev_new(t_symbol *s, long argc, t_atom *argv);void time_dev_free(t_time_dev *x);void time_dev_assist(t_time_dev *x, void *b, long m, long a, char *s);void time_dev_bang(t_time_dev *x);void ext_main(void *r) {    t_class *c;    c = class_new("time_dev", (method)time_dev_new, (method)time_dev_free, sizeof(t_time_dev), 0L, A_GIMME, 0);    class_addmethod(c, (method)time_dev_bang, "bang", 0);    class_addmethod(c, (method)time_dev_assist, "assist", A_CANT, 0);    class_register(CLASS_BOX, c);}void *time_dev_new(t_symbol *s, long argc, t_atom *argv) {    t_time_dev *x = (t_time_dev *)object_alloc(/* class pointer */);    if(x){        x->n = 3;        // default: 3 eigenstates        x->tsteps = 100; // default time resolution        x->tmin = 0.0;        x->tmax = 1.0;        if(argc >= 1)            x->n = atom_getlong(argv);        if(argc >= 2)            x->tsteps = atom_getlong(argv+1);        if(argc >= 3)            x->tmin = atom_getfloat(argv+2);        if(argc >= 4)            x->tmax = atom_getfloat(argv+3);        x->out = outlet_new(x, NULL);    }    return (x);}void time_dev_free(t_time_dev *x) { }void time_dev_assist(t_time_dev *x, void *b, long m, long a, char *s) {    if(m==1)        sprintf(s, "Bang to compute time development");    else        sprintf(s, "Output: List of magnitude & phase pairs per time step");}void time_dev_bang(t_time_dev *x) {    long n = x->n;    long tsteps = x->tsteps;    double tmin = x->tmin;    double tmax = x->tmax;    double dt = (tmax - tmin) / (tsteps - 1);        // For demonstration, use dummy eigenvalues and coefficients.    double *eigenvalues = (double *)sysmem_newptr(n * sizeof(double));    double complex *coeff = (double complex *)sysmem_newptr(n * sizeof(double complex));    for(long k = 0; k < n; k++){        eigenvalues[k] = 1.0 + k;      // dummy eigenvalues        coeff[k] = 1.0 + 0.0*I;        // dummy coefficients    }    // For each time step, compute: Ψ(t) = Σ_k coeff[k] * exp(-i * eigenvalues[k] * t)    t_atom *out_list = (t_atom *)sysmem_newptr(tsteps * 2 * sizeof(t_atom)); // each time step: mag and phase    for(long t = 0; t < tsteps; t++){        double time = tmin + t * dt;        double complex psi = 0.0 + 0.0*I;        for(long k = 0; k < n; k++){            psi += coeff[k] * cexp(-I * eigenvalues[k] * time);        }        double mag = cabs(psi);        double phase = carg(psi);        atom_setfloat(out_list + t*2, mag);        atom_setfloat(out_list + t*2 + 1, phase);    }    outlet_list(x->out, gensym("list"), tsteps*2, out_list);    sysmem_freeptr(out_list);    sysmem_freeptr(eigenvalues);    sysmem_freeptr(coeff);}