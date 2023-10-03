#include <stdio.h>
#include <stdlib.h>
#include <windows.h>
#include <time.h>
#include "C:\Users\Sandro\Desktop\uni\TO\codes\nlopt\nlopt\src\api\nlopt.h"

#define SIZE (46)
#define LIB_NAME "One_Fractions_new.dll"
#define CALC_NAME "CALCULATE_S"
#define CALC_PRINT_NAME "CALCULATE_S_WITHSAVE"
#define false (0)
#define true (1)
#define bool char
#define dxr (1e-6)
#define dxw (4e-6)
#define dxd (4e-6)
#define alpha (1e-8)

typedef void (*SUB)(double* R, double* W, double* D, int* LAY_M, double* S_PAR);
typedef void (*SUB_WS)(double* R, double* W, double* D, int* LAY_M);

double obj_f(unsigned n, const double* vars, double* grad, void* data);
void constr_f(unsigned m, double *result, unsigned n, const double* x, double* grad, void* f_data);
int init_libraries(HINSTANCE hDLL);
void init_ideal();
void init_init(double* vars, int type);
void init_bounds(double* lb, double* ub);

SUB calc_s;
SUB_WS calc_s_ws;
const double A = 5.075e-3;
int LAY_M = 4;
double* S_PAR;
// double ideal[] = {0.002176, 0.004369, 0.009583, 0.023721, 0.069401, 0.088145,
//                  0.112899, 0.145750, 0.189398, 0.247112, 0.322319, 0.417455,
//                  0.531767, 0.658455, 0.783347, 0.888349, 0.959585, 0.994013,
//                  0.994986, 0.969038, 0.940275, 0.921247, 0.916528, 0.926231,
//                  0.947418, 0.973853, 0.995152, 0.997085, 0.965412, 0.917989,
//                  0.851448, 0.770967, 0.683565, 0.596071, 0.513684, 0.439498,
//                  0.374760, 0.319441, 0.272788, 0.233733, 0.201145, 0.122289,
//                  0.068552, 0.042943, 0.029354, 0.021492};
// double ideal[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.707106,
//                   0.99999, 0.99999, 0.99999, 0.99999, 0.99999, 0.99999,
//                   0.99999, 0.99999, 0.99999, 0.99999, 0.99999, 0.99999,
//                   0.99999, 0.99999, 0.99999, 0.707106, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//                   0, 0, 0, 0, 0};
double ideal[SIZE];
int num_obj = 0;
double cur_best = -1;
long long start, end;

int main() {
    HINSTANCE hDLL = LoadLibraryA(LIB_NAME);
    int load_lib_res = init_libraries(hDLL);
    if(load_lib_res == 1) {
        printf("error loading %s. error code: %lu", LIB_NAME, GetLastError());
        return 1;
    }
    else if (load_lib_res == 2) {
        printf("error loading %s. error code: %lu", CALC_NAME, GetLastError());
        FreeLibrary(hDLL);
        return 1;
    }
    else if (load_lib_res == 3) {
        printf("error loading %s. error code: %lu", CALC_PRINT_NAME, GetLastError());
        FreeLibrary(hDLL);
        return 1;
    }  
    init_ideal();
    S_PAR = malloc(sizeof(double) * SIZE);
    if(S_PAR == NULL) {
        printf("could not allocate memory for S_PAR, exiting!");
        FreeLibrary(hDLL);
        return 1;
    }
    unsigned n = 3*LAY_M-1;
    
    double vars[n];
    init_init(vars, 1);

    double lb[n], ub[n];
    init_bounds(lb, ub);
    lb[0] = lb[1] = 1e-6;
    ub[0] = ub[1] = A/2;

    nlopt_opt nlo = nlopt_create(NLOPT_LD_MMA, n);
    nlopt_result nlo_res;
    nlo_res = nlopt_set_lower_bounds(nlo, lb);
    if(nlo_res != NLOPT_SUCCESS) {
        printf("error %d setting lower bounds", nlo_res);
        FreeLibrary(hDLL);
        free(S_PAR);
        return 1;
    }
    nlo_res = nlopt_set_upper_bounds(nlo, ub);
    if(nlo_res != NLOPT_SUCCESS) {
        printf("error %d setting upper bounds", nlo_res);
        FreeLibrary(hDLL);
        free(S_PAR);
        return 1;
    }
    unsigned m = 3*LAY_M-1;
    const double tol = 1e-14;
    nlo_res = nlopt_set_max_objective(nlo, obj_f, NULL);
    if(nlo_res != NLOPT_SUCCESS) {
        printf("could not set objective\n");
        FreeLibrary(hDLL);
        free(S_PAR);
        return 1;
    }
    nlo_res = nlopt_add_inequality_mconstraint(nlo, m, constr_f, NULL, &tol);
    if(nlo_res != NLOPT_SUCCESS) {
        printf("could not add inequality constraints\n");
        FreeLibrary(hDLL);
        free(S_PAR);
        return 1;
    }

    double cost = -1;
    printf("starting optimization.\n\n");
    
    start = clock();
    if (nlopt_optimize(nlo, vars, &cost) < 0) {
        printf("nlopt failed!\n");
    } else {
        printf("found maximum at:\n");
        for(int i = 0; i < LAY_M; i++) {
            printf("r%d=%f w%d=%f d%d=%f\n", i+1, vars[i], i+1, vars[i+LAY_M], i+1, vars[i+2*LAY_M-1]);
        }
    }
    end = clock();
    double timeTaken = (end - start) / CLOCKS_PER_SEC;
    printf("elapsed time: %f\n", timeTaken);
    printf("function evaluations took: %d\n\n", num_obj);
    printf("S_PAR:\n");
    for(int i = 0; i < SIZE; i++) {
        printf("%f, ", S_PAR[i]);
    }
    printf("\ncost was %f", cost);
    free(S_PAR);
    nlopt_destroy(nlo);
    FreeLibrary(hDLL);
    return 0;
}

double obj_f(unsigned n, const double* vars, double* grad, void* data) {
    double* r = (double*)vars;
    double* w = (double*)r + LAY_M;
    double d[LAY_M];
    d[0] = 5e-3;
    for(int i = 1; i < LAY_M; i++) {
        d[i] = *((double*)vars + 2*LAY_M + i-1);
    }

    calc_s(r, w, d, &LAY_M, S_PAR);
    
    double res = 0;
    double temp;
    for(int i = 0; i < SIZE; i++) {
        temp = S_PAR[i] - ideal[i];
        res += temp*temp;
    }
    
    res = 1 / (res + alpha);
    if(res > cur_best) {
        for(int i = 0; i < LAY_M; i++) {
            printf("r=%f w=%f d=%f\n", r[i], w[i], d[i]);
        }
        printf("\nnumber of evals: %d\n", num_obj);
        printf("objective: %f\n", res);
        printf("time since the beginning: %f seconds\n\n", ((double)(clock() - start)) / CLOCKS_PER_SEC);
        cur_best = res;
    }
    if(num_obj%20 == 0) {
        printf("\nnumber of evals: %d\n", num_obj);
        printf("time since the beginning: %f seconds\n\n", ((double)(clock() - start)) / CLOCKS_PER_SEC);
    }
    num_obj++;

    if(grad) {
        double temp[n];
        memcpy(temp, vars, n * sizeof(double));
        for(int i = 0; i < LAY_M; i++) {
            temp[i] += dxr;
            double newObj = obj_f(n, temp, NULL, NULL);
            grad[i] = (newObj-res) / dxr;
            temp[i] = vars[i];
        }
        for(int i = LAY_M; i < 2*LAY_M; i++) {
            temp[i] += dxw;
            double newObj = obj_f(n, temp, NULL, NULL);
            grad[i] = (newObj-res) / dxw;
            temp[i] = vars[i];
        }
        for(int i = 2*LAY_M; i < n; i++) {
            temp[i] += dxd;
            double newObj = obj_f(n, temp, NULL, NULL);
            grad[i] = (newObj-res) / dxd;
            temp[i] = vars[i];
        }
    }
    return res;
}

/*
imports functions from specified DLL (Dinamically Linked Library),
returns 0 if import was successful,
        1 if encountered problem loading library LIB_NAME
        2 if encountered problem loading function CALC_NAME from library
        3 if encountered problem loading function CALC_PRINT_NAME from library
*/
int init_libraries(HINSTANCE hDLL) {
    hDLL = LoadLibraryA(LIB_NAME);
    if (hDLL == NULL) return 1;

    calc_s = (SUB)GetProcAddress(hDLL, CALC_NAME);
    if(calc_s == NULL) return 2;
    
    calc_s_ws = (SUB_WS)GetProcAddress(hDLL, CALC_PRINT_NAME);
    if(calc_s_ws == NULL) return 3;

    return 0;
}

void init_ideal() {
    memset(ideal, 0, SIZE * sizeof(double));
    ideal[14] = ideal[33] = 0.707106;
    for (int i = 15; i < 33 ; i++) {
        ideal[i] = 1.1;
    }
}

void init_init(double* vars, int type) {
    switch(type) {
        case 0:
            for(int i = 0; i < LAY_M; i++) {
                vars[i] = 0.4e-3; 
            }
            for(int i = LAY_M; i < 2*LAY_M; i++) {
                vars[i] = A/2; 
            }
            for(int i = 2*LAY_M; i < 3*LAY_M-1; i++) {
                vars[i] = 5e-3; 
            }
            break;
        case 1:
            for(int i = 0; i < LAY_M; i++) {
                // vars[i] = 0.4e-3 - (i%2)*0.1e-3;
            }
            for(int i = LAY_M; i < 2*LAY_M; i++) {
                // vars[i] = A/2 + (i%2)*A/4 - ((i+1)%2)*A/4;
            }
            for(int i = 2*LAY_M; i < 3*LAY_M-1; i++) {
                // vars[i] = 4e-3;
            }
            vars[0] = vars[3] = (0.7733e-3)/2;
            vars[1] = vars[2] = (0.7498e-3)/2;
            vars[4] = vars[7] = A/2 - 1.0987e-3;
            vars[5] = vars[6] = A/2 + 0.021e-3;
            vars[8] = vars[10] = 4.611e-3;
            vars[9] = 5.075e-3;
            break;
    }
    
}

void init_bounds(double* lb, double* ub) {
    for(int i = 0; i < LAY_M; i++) {
        lb[i] = 1e-5;
        ub[i] = A/2 - 1e-5;
    }
    for(int i = LAY_M; i < 2*LAY_M; i++) {
        lb[i] = 1e-5;
        ub[i] = A-1e-5;
    }
    for(int i = 2*LAY_M; i < 3*LAY_M-1; i++) {
        lb[i] = 1e-5;
        ub[i] = 10e-3;
    }
}

void constr_f(unsigned m, double *result, unsigned n, const double* x, double* grad, void* f_data) {
    double* r = (double*)x;
    double* w = (double*)r + LAY_M;
    double* d = (double*)r + 2*LAY_M;

    for(int i = 0; i < LAY_M-1; i++) {
        result[i] = (r[i] + r[i+1]) - d[i];
        if(grad) {
            memset(grad + (i*n), 0, n);
            grad[i*n + i] = 1;
            grad[i*n + i+1] = 1;
            grad[i*n + 2*LAY_M + i] = -1;
        }
    }
    for(int i = 0; i < LAY_M; i++) {
        result[i + LAY_M-1] = r[i]-w[i];
        if(grad) {
            memset(grad + ((i + LAY_M-1)*n), 0, n);
            grad[(i + LAY_M-1)*n + i] = 1;
            grad[(i + LAY_M-1)*n + LAY_M + i] = -1;
        }
    }
    for(int i = 0; i < LAY_M; i++) {
        result[i+2*LAY_M-1] = r[i]+w[i]-A;
        if(grad) {
            memset(grad+((i + 2*LAY_M-1)*n), 0, n);
            grad[(i + 2*LAY_M-1)*n + i] = 1;
            grad[(i + 2*LAY_M-1)*n + LAY_M + i] = 1;
        }
    }

}