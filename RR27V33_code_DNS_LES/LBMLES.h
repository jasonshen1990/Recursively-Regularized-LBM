#ifndef LBMLES
#define LBMLES

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <sys/time.h>

#include "parameters.h"

// Device constants declaration
extern __constant__ int d_cix[NPOP];
extern __constant__ int d_ciy[NPOP];
extern __constant__ int d_ciz[NPOP];
extern __constant__ double d_tp[NPOP];
extern __constant__ double d_tau;
extern __constant__ double d_visc;
extern __constant__ double d_ustar, d_ystar;
extern __constant__ double d_force_in_x;
extern __device__ double d_rand_phase;
extern __device__ int d_istep;
extern __constant__ double d_pi, d_pi2;
extern __constant__ double d_RT, d_cs;
extern __constant__ double d_sigma;

extern __constant__ double d_H100[NPOP], d_H010[NPOP], d_H001[NPOP];
extern __constant__ double d_H110[NPOP], d_H101[NPOP], d_H011[NPOP], d_H200[NPOP], d_H020[NPOP], d_H002[NPOP];
extern __constant__ double d_H210[NPOP], d_H201[NPOP], d_H021[NPOP], d_H120[NPOP], d_H102[NPOP], d_H012[NPOP], d_H111[NPOP];
// extern __constant__ double d_H220[NPOP], d_H202[NPOP], d_H022[NPOP], d_H211[NPOP], d_H121[NPOP], d_H112[NPOP];
// extern __constant__ double d_H221[NPOP], d_H212[NPOP], d_H122[NPOP], d_H222[NPOP];

// Global variables declaration
extern double visc, tau;
extern double u0, ystar;
extern double force_in_x;
extern int continue_step;
extern double pi, pi2;
extern double RT, cs;
extern double sigma;

// Device pointers declaration
extern double *d_f, *d_rho, *d_ux, *d_uy, *d_uz;
extern double *d_force_realx, *d_force_realy, *d_force_realz;
extern double *d_f_temp;
extern double *d_taut;    // LES

// Host pointers declaration
extern double *h_rho, *h_ux, *h_uy, *h_uz;
extern double *h_f;
extern double *h_taut;    // LES

// Function declarations
void initialize();

__global__ void init_arrays(double *f, double *rho, double *ux, double *uy, double *uz, 
                            double *force_realx, double *force_realy, double *force_realz);

__global__ void LES_WALE(double *ux, double *uy, double *uz, double *taut);

__global__ void collision_RR(double *f, double *force_realx, double *force_realy, double *force_realz, 
                             double *rho, double *ux, double *uy, double *uz, double *taut);

__global__ void streaming(double *f, double *f_temp);
__global__ void BB(double *f, double *f_temp);

__global__ void RRBC(double *f, double *force_realx, double *force_realy, double *force_realz, 
                     double *rho, double *ux, double *uy, double *uz, double *taut);
__global__ void FDBC(double *f, double *force_realx, double *force_realy, double *force_realz, 
                     double *rho, double *ux, double *uy, double *uz, double *taut);

__global__ void macrovar(double *f, double *force_realx, double *force_realy, double *force_realz, 
                        double *rho, double *ux, double *uy, double *uz);
__global__ void macrovar_BN(double *f, double *force_realx, double *force_realy, double *force_realz, 
                        double *rho, double *ux, double *uy, double *uz);
__global__ void forcing(double *force_realx, double *force_realy, double *force_realz);
__global__ void forcing_GP(double *force_realx, double *force_realy, double *force_realz);

void diag_flow(int istep);
void output_flow(int istep);
void output_profile(int istep);
void output_f(int istep);

// Error checking macro
#define CHECK_CUDA_ERROR(call) \
do { \
    cudaError_t err = call; \
    if (err != cudaSuccess) { \
        fprintf(stderr, "CUDA error in %s:%d: %s\n", __FILE__, __LINE__, \
                cudaGetErrorString(err)); \
        exit(EXIT_FAILURE); \
    } \
} while(0)

#endif 