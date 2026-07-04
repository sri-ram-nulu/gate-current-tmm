/*
This code makes up for the 2nd part of the project. This code requires two input files in it's directory.
    1. input_parameters.txt:
        Contains the gate voltage, oxide thickness, channel thickness, m_gate/me, m_oxide/m_e, m_channel/m_e, no_of_slices_btw_Emin_Emax
        Emin and Emax are the Energy ranges to sweep and calculate T(E).

    2. Eminmax :
        Contains Emin and Emax values to define range (Emin, Emax) to compute T(E)

    The above 2 files can be used by any part of the project.
    It can be further modified with prior information with the other members of the team.
    
    3. poisson_sio2.txt:
        It contains the electrostatic potential profile accross the oxide.
        This is output textfile generated after running the poisson.c file.

The above three files are essential to run this code.
This program outputs "transmission_data.txt" containing Total Energy of Electron and the corresponding T(E).

This code also had standard tests (1, 2, 3) like for rectangular barriers and triangular barriers. The plots can be seen in the T_E_test_{i} folder.
Observations made: 
    1. Transmission Probability is never greater than 1.
    2. Transmission Probability is same for test case 2 & 3, which demonstrates that probability of an electron moving 
        from left to the right has the same probability to move from right to them left in the same environment.

Code Author: Sri Ram Nulu
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_SCAN 1e7
#define h_bar 1.0545718e-34 // Reduced planck's constant in J-s
#define m_e 9.10938356e-31 // Mass of electron in kg
#define eV 1.60218e-19  // eV to Joules conversion

#define electron  (1.6e-19)
const double Kt = (electron * 25.9e-3);
const double Nc =2.8e25 ;
const double Nv =1.04e25;
const double Eg = (1.12*electron);
const double Esisio2 = 3.1*electron;
const double epsi = 11.7*8.85e-12;
const double epsio2 = 3.9*8.85e-12;
const double intrinsic2 = 2.25e32;

// Complex number type
typedef struct {
    double r;  // real part
    double i;  // imaginary part
} Complex;

// 2x2 complex matrix
typedef struct {
    Complex m[2][2];
} cmat2;

// Complex arithmetic
Complex cadd(Complex a, Complex b) {
    return (Complex){a.r + b.r, a.i + b.i};
}

Complex cmul(Complex a, Complex b) {
    return (Complex){a.r*b.r - a.i*b.i, a.r*b.i + a.i*b.r};
}

Complex cdiv(Complex a, Complex b) {
    double denom = b.r*b.r + b.i*b.i;
    return (Complex){(a.r*b.r + a.i*b.i)/denom, (a.i*b.r - a.r*b.i)/denom};
}

Complex cexpj(double x) {
    return (Complex){cos(x), sin(x)};
}

double c_abs2(Complex a) {
    return a.r*a.r + a.i*a.i;
}

double c_abs(Complex a) {
    return sqrt(c_abs2(a));
}

// Matrix operations
cmat2 mmul(cmat2 A, cmat2 B) {
    cmat2 C;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            Complex sum = {0, 0};
            for (int k = 0; k < 2; k++) {
                sum = cadd(sum, cmul(A.m[i][k], B.m[k][j]));
            }
            C.m[i][j] = sum;
        }
    }
    return C;
}

Complex cdet(cmat2 A) {
    return cadd(cmul(A.m[0][0], A.m[1][1]),
               (Complex){-cmul(A.m[0][1], A.m[1][0]).r, -cmul(A.m[0][1], A.m[1][0]).i});
}

// Physics functions
Complex wave_number(double E, double V, double m_reg) {
    double arg = 2 * m_reg * (E - V) / (h_bar * h_bar);
    if (arg >= 0) {
        return (Complex){sqrt(arg), 0};
    } else {
        return (Complex){0, sqrt(-arg)};
    }
}

cmat2 M0(Complex k, double del_x) {
    cmat2 M;
    Complex e1 = cexpj(k.r * del_x);  // if k is real
    if (k.i != 0) e1 = (Complex){exp(-k.i * del_x), 0}; // if k is imaginary

    Complex e2 = cdiv((Complex){1, 0}, e1);

    M.m[0][0] = e1;              M.m[0][1] = (Complex){0, 0};
    M.m[1][0] = (Complex){0, 0}; M.m[1][1] = e2;
    return M;
}

cmat2 Ms(Complex k2, Complex k1, double mp, double mm) {
    Complex ratio = cdiv((Complex){k1.r * mp, k1.i * mp},
                        (Complex){k2.r * mm, k2.i * mm});

    Complex one = {1, 0};
    Complex m_one = {-1, 0};
    Complex half = {0.5, 0};

    Complex ratio1 = cadd(one, ratio);              // (1 + ratio)
    Complex ratio2 = cadd(one, cmul(m_one, ratio)); // (1 - ratio)

    cmat2 M = {{
        {cmul(half, ratio1), cmul(half, ratio2)},
        {cmul(half, ratio2), cmul(half, ratio1)}
    }};

    return M;
}

cmat2 overall_transfer_matrix(int N, double E, double* potentials, double* lengths, double* m_stars, Complex* k_vals){

    // Calculate wave numbers
    for (int j = 0; j < N + 2; j++) {
        k_vals[j] = wave_number(E, potentials[j], m_stars[j]);
    }

    // Initialize M as identity matrix
    cmat2 M = {{{{1, 0}, {0, 0}}, {{0, 0}, {1, 0}}}};

    // First interface
    M = mmul(Ms(k_vals[1], k_vals[0], m_stars[1], m_stars[0]), M);

    // Loop through regions
    for (int i = 1; i < N + 1; i++) {
        M = mmul(M0(k_vals[i], lengths[i-1]), M);
        if (i < N) {
            M = mmul(Ms(k_vals[i+1], k_vals[i], m_stars[i+1], m_stars[i]), M);
        } else {
            M = mmul(Ms(k_vals[N+1], k_vals[i], m_stars[N+1], m_stars[i]), M);
        }
    }

    return M;
}

void reflection_transmission(int N, cmat2 M, Complex* k_vals, double* m_stars, Complex* r, Complex* t, double* R, double* T) {
    *r = cdiv(M.m[1][0], M.m[1][1]);
    *t = cdiv(cdet(M), M.m[1][1]);

    double ini_k_by_m = c_abs(k_vals[0]) / m_stars[0];
    double fin_k_by_m = c_abs(k_vals[N+1]) / m_stars[N+1];
    *T = c_abs2(*t) * fin_k_by_m / ini_k_by_m;
    *R = c_abs2(*r);
}



int main(){
    /* Test cases:
        0: Uses Poisson Data (Default)
        1: Uses Rectangular Barrier (Testing Only)
        2: Uses Rising Triangular Barrier
        3: Uses Falling Triangular Barrier (Mirror of Rising)
    */

    int test = 0;

    char * potentials_file = "poisson_sio2.txt";
    char * Eminmax_file = "Eminmax";
    char * input_parameters_file = "input_parameters.txt";

    if(test == 1){ //Rectangular Barrier
        potentials_file = "./T_E_test_1/T_E_test_1_potentials.txt";
        Eminmax_file = "./T_E_test_1/T_E_test_1_Eminmax.txt";
        input_parameters_file = "./T_E_test_1/T_E_test_1_input_parameters.txt";
    }

    if(test == 2){ // Triangle 1
        potentials_file = "./T_E_test_2/T_E_test_2_potentials.txt";
        Eminmax_file = "./T_E_test_2/T_E_test_2_Eminmax.txt";
        input_parameters_file = "./T_E_test_2/T_E_test_2_input_parameters.txt";
    }

    if(test == 3){ // Triangle 2 (Mirrored)
        potentials_file = "./T_E_test_3/T_E_test_3_potentials.txt";
        Eminmax_file = "./T_E_test_3/T_E_test_3_Eminmax.txt";
        input_parameters_file = "./T_E_test_3/T_E_test_3_input_parameters.txt";
    }


    // 1. Declaring Parameters
    double barrier_width; // Oxide thickness

    double m_left; // Channel
    double m_barrier; // Oxide
    double m_right; // Gate/metal

    double Emin, Emax; // For Energy Sweep
    int steps_for_E;

    int N; // No. of slices

	FILE *minmax;
	minmax = fopen(Eminmax_file,"r");
	fscanf(minmax,"%lf %lf",&Emin,&Emax);
    //printf("%lf %lf\n", Emin, Emax);
	fclose(minmax);


    double gate_voltage; // not used
    double channel_width; // not used
    
    FILE *file_params;
    file_params = fopen(input_parameters_file, "r");
    
    if (file_params == NULL) { printf("Error: Could not open parameters.txt\n"); return 1;}
    int param_count = fscanf(file_params,"%lf,  %lf, %lf, %lf, %lf, %lf, %d", &gate_voltage, &barrier_width, &channel_width, &m_right, &m_barrier, &m_left, &steps_for_E);
    m_left *= m_e; m_barrier *= m_e; m_right *= m_e; // masses in kg

    fclose(file_params);


    // 2. Fetching the potentials in the oxide
    double *potential_E = NULL;

    int pot_count = 0;
    int index;
    double val1, val2;
    FILE *file = fopen(potentials_file, "r");
    if (!file) { perror(potentials_file); return 1;}
    while (pot_count < MAX_SCAN && fscanf(file, "%d %le %le", &index, &val1, &val2) == 3) {
        double *tmp = realloc(potential_E, (pot_count + 1) * sizeof(double));
        if (!tmp) {perror("realloc"); free(potential_E); fclose(file); return 1; }
        potential_E = tmp;
        potential_E[pot_count++] = val2;
    }
    N = pot_count-2;

    // 3. Updating mstars, lengths, wave numbers
    double *mstars     = malloc((N+2) * sizeof(double));
    double *lengths    = malloc(N * sizeof(double));
    Complex *kvals     = malloc((N+2) * sizeof(Complex));

    mstars[0] = m_left;
    for (int i=0;i<N;i++) mstars[i+1] = m_barrier;
    mstars[N+1] = m_right;

    double slice_width = barrier_width / N;
    for (int i=0;i<N;i++) lengths[i] = slice_width;


    // 4. Generating the Transmission Data file
    FILE *fp = fopen("transmission_data.txt", "w");
    if (!fp) {
        printf("Error: could not open file\n");
        return 1;
    }

    // Energy sweep
    //int steps_for_E = 3000;
    //double Emin = 0.1, Emax = 8.0;
    for (int idx=0; idx<steps_for_E; idx++) {
        double E_eV = Emin + (Emax-Emin)*idx/(steps_for_E-1);
        double E = E_eV * eV;

        cmat2 M;
        M = overall_transfer_matrix(N, E, potential_E, lengths, mstars, kvals);

        Complex r, t;
        double R, T;
        reflection_transmission(N, M, kvals, mstars, &r, &t, &R, &T);

        fprintf(fp, "%f\t%g\n", E_eV, T);
    }

    fclose(fp);
    printf("Results written to transmission_data.txt\n");

    free(potential_E);
    free(mstars);
    free(lengths);
    free(kvals);

    return 0;
}

