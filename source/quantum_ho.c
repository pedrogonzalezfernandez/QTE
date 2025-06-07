/* quantum_ho.c – Quantum Harmonic Oscillator external for Max/MSP
 *
 * This external computes the Hamiltonian matrix H = 0.5*(P^2 + Q^2)
 * where:
 *    - PImpulse = diag(0, 1, ..., n-1)
 *    - F is the Fourier matrix: F[k][l] = (1/sqrt(n)) * exp(2πi*k*l/n)
 *    - Finv is the conjugate transpose of F
 *    - P = Finv * (PImpulse * F)
 *    - Q = diag( -((n-1)*a/2) + i )  (with i = row index)
 *    - Finally, H = 0.5*(P^2 + Q^2) (with Q^2 added only on the diagonal)
 *
 * The result is output as a flat list of 2*n*n numbers (each matrix entry’s
 * real and imaginary parts in sequence).
 *
 * Compile with Xcode using the Max SDK.
 */

#include "ext.h"
#include "ext_obex.h"
#include <math.h>
#include <stdlib.h>
#include <complex.h>

// Define our object structure.
typedef struct _quantum_ho {
    t_object ob;
    long n;      // Matrix dimension.
    double a;    // Potential parameter.
    void *out;   // Outlet pointer.
} t_quantum_ho;

// Global class pointer.
static t_class *quantum_ho_class = NULL;

/* Function prototypes */
void *quantum_ho_new(t_symbol *s, long argc, t_atom *argv);
void quantum_ho_free(t_quantum_ho *x);
void quantum_ho_assist(t_quantum_ho *x, void *b, long m, long a, char *s);
void quantum_ho_bang(t_quantum_ho *x);

// Helper: Allocate an n x n matrix of double complex numbers.
static double complex **alloc_complex_matrix(long n) {
    double complex **matrix = (double complex **)malloc(n * sizeof(double complex *));
    if (!matrix)
        return NULL;
    for (long i = 0; i < n; i++) {
        matrix[i] = (double complex *)calloc(n, sizeof(double complex));
        if (!matrix[i]) {
            for (long j = 0; j < i; j++)
                free(matrix[j]);
            free(matrix);
            return NULL;
        }
    }
    return matrix;
}

// Helper: Free an n x n complex matrix.
static void free_complex_matrix(double complex **matrix, long n) {
    if (matrix) {
        for (long i = 0; i < n; i++)
            free(matrix[i]);
        free(matrix);
    }
}

// Helper: Compute the Fourier matrix F (n x n) where:
// F[k][l] = (1/sqrt(n)) * exp(2πi*k*l/n)
static void compute_fourier_matrix(double complex **F, long n) {
    double norm = 1.0 / sqrt((double)n);
    for (long k = 0; k < n; k++) {
        for (long l = 0; l < n; l++) {
            double angle = 2.0 * M_PI * k * l / n;
            F[k][l] = norm * (cos(angle) + I * sin(angle));
        }
    }
}

// Helper: Compute the conjugate transpose of F: Finv[i][j] = conj(F[j][i]).
static void compute_conjugate_transpose(double complex **F, double complex **Finv, long n) {
    for (long i = 0; i < n; i++) {
        for (long j = 0; j < n; j++) {
            Finv[i][j] = conj(F[j][i]);
        }
    }
}

// Helper: Multiply a diagonal matrix D (vector of length n) with matrix F.
static void multiply_diag_matrix(double *D, double complex **F, double complex **result, long n) {
    for (long i = 0; i < n; i++) {
        for (long j = 0; j < n; j++) {
            result[i][j] = D[i] * F[i][j];
        }
    }
}

// Helper: Multiply two n x n complex matrices A and B.
static void multiply_complex_matrices(double complex **A, double complex **B, double complex **result, long n) {
    for (long i = 0; i < n; i++) {
        for (long j = 0; j < n; j++) {
            double complex sum = 0.0 + 0.0 * I;
            for (long k = 0; k < n; k++)
                sum += A[i][k] * B[k][j];
            result[i][j] = sum;
        }
    }
}

// Helper: Round a double to 5 decimal places.
static double round5(double x) {
    return round(x * 100000.0) / 100000.0;
}

// Compute the Hamiltonian matrix H = 0.5*(P^2 + Q^2).
static double complex **compute_hamiltonian(t_quantum_ho *x) {
    long n = x->n;
    double a = x->a;
    double complex **H = alloc_complex_matrix(n);
    if (!H)
        return NULL;
    
    double complex **F = alloc_complex_matrix(n);
    double complex **Finv = alloc_complex_matrix(n);
    double complex **P_temp = alloc_complex_matrix(n); // PImpulse * F
    double complex **P = alloc_complex_matrix(n);      // Finv * (PImpulse * F)
    double complex **P2 = alloc_complex_matrix(n);     // P^2
    if (!F || !Finv || !P_temp || !P || !P2) {
        free_complex_matrix(H, n);
        free_complex_matrix(F, n);
        free_complex_matrix(Finv, n);
        free_complex_matrix(P_temp, n);
        free_complex_matrix(P, n);
        free_complex_matrix(P2, n);
        return NULL;
    }
    
    compute_fourier_matrix(F, n);
    compute_conjugate_transpose(F, Finv, n);
    
    // Create PImpulse as a diagonal vector [0, 1, 2, ..., n-1].
    double *PImpulse = (double *)malloc(n * sizeof(double));
    for (long i = 0; i < n; i++)
        PImpulse[i] = (double)i;
    
    // Compute P_temp = PImpulse * F.
    multiply_diag_matrix(PImpulse, F, P_temp, n);
    
    // Compute P = Finv * P_temp.
    multiply_complex_matrices(Finv, P_temp, P, n);
    
    // Compute P^2 = P * P.
    multiply_complex_matrices(P, P, P2, n);
    
    // Create Q diagonal: Q[i] = -((n-1)*a/2) + i.
    double *Q_diag = (double *)malloc(n * sizeof(double));
    for (long i = 0; i < n; i++)
        Q_diag[i] = -((n - 1) * a / 2.0) + i;
    
    // Compute H = 0.5*(P^2 + Q^2) (Q^2 is only added on the diagonal).
    for (long i = 0; i < n; i++) {
        for (long j = 0; j < n; j++) {
            double complex qterm = 0.0 + 0.0 * I;
            if (i == j) {
                double q2 = Q_diag[i] * Q_diag[i];
                qterm = q2 + 0.0 * I;
            }
            H[i][j] = 0.5 * (P2[i][j] + qterm);
            // Round both real and imaginary parts.
            double real_part = round5(creal(H[i][j]));
            double imag_part = round5(cimag(H[i][j]));
            H[i][j] = real_part + I * imag_part;
        }
    }
    
    free(PImpulse);
    free(Q_diag);
    free_complex_matrix(F, n);
    free_complex_matrix(Finv, n);
    free_complex_matrix(P_temp, n);
    free_complex_matrix(P, n);
    free_complex_matrix(P2, n);
    
    return H;
}

/* ---------------------- Max External Methods ------------------------- */

// The main entry point called by Max at load time.
void ext_main(void *r)
{
    t_class *c = class_new("quantum_ho", (method)quantum_ho_new, (method)quantum_ho_free,
                             sizeof(t_quantum_ho), 0L, A_GIMME, 0);
    
    class_addmethod(c, (method)quantum_ho_bang, "bang", 0);
    class_addmethod(c, (method)quantum_ho_assist, "assist", A_CANT, 0);
    class_register(CLASS_BOX, c);
    quantum_ho_class = c;
}

// Create a new instance of the object.
void *quantum_ho_new(t_symbol *s, long argc, t_atom *argv)
{
    t_quantum_ho *x = (t_quantum_ho *)object_alloc(quantum_ho_class);
    if (x) {
        // Set default values.
        x->n = 8;      // default dimension.
        x->a = 1.0;    // default potential parameter.
        if (argc >= 1) {
            if (atom_gettype(argv) == A_LONG)
                x->n = atom_getlong(argv);
            else if (atom_gettype(argv) == A_FLOAT)
                x->n = (long)atom_getfloat(argv);
        }
        if (argc >= 2)
            x->a = atom_getfloat(argv + 1);
        x->out = outlet_new(x, NULL);
    }
    return x;
}

// Free the object (nothing to free in this example).
void quantum_ho_free(t_quantum_ho *x)
{
    ;
}

// Provide assistance messages for inlets/outlets.
void quantum_ho_assist(t_quantum_ho *x, void *b, long m, long a, char *s)
{
    if(m == 1)
        sprintf(s, "Bang to compute Hamiltonian matrix");
    else
        sprintf(s, "Output: Flattened Hamiltonian matrix as list (real, imag pairs)");
}

// Method invoked on a bang: compute and output the Hamiltonian matrix.
void quantum_ho_bang(t_quantum_ho *x)
{
    long n = x->n;
    double complex **H = compute_hamiltonian(x);
    if (!H) {
        object_error((t_object *)x, "Memory allocation failed during Hamiltonian computation");
        return;
    }
    long list_size = 2 * n * n;
    t_atom *out_list = (t_atom *)sysmem_newptr(list_size * sizeof(t_atom));
    for (long i = 0; i < n; i++) {
        for (long j = 0; j < n; j++) {
            long index = 2 * (i * n + j);
            atom_setfloat(out_list + index, creal(H[i][j]));
            atom_setfloat(out_list + index + 1, cimag(H[i][j]));
        }
    }
    outlet_list(x->out, gensym("list"), list_size, out_list);
    sysmem_freeptr(out_list);
    free_complex_matrix(H, n);
}
