#ifndef RB3D_H
#define RB3D_H

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
extern __constant__ double d_grav0;
extern __constant__ double d_beta;
extern __constant__ double d_t0;
extern __constant__ double d_visc;
extern __constant__ double d_ustar, d_ystar;
extern __constant__ double d_force_in_x;
extern __device__ double d_rand_phase;
extern __device__ int d_istep;
extern __constant__ double d_pi, d_pi2;

// Global variables declaration
// extern double rayl, prand;
extern double visc, tau, diff, tauc, tauci, grav0, beta;
extern double u0, ystar;
extern double force_in_x;
extern double t0;
extern int continue_step;
extern double pi, pi2;

// Device pointers declaration
extern double *d_f, *d_rho, *d_ux, *d_uy, *d_uz;
extern double *d_force_realx, *d_force_realy, *d_force_realz;
extern double *d_f_temp; // *d_g_temp;
extern double *d_uyT;
extern double *d_Tx, *d_Ty, *d_Tz, *d_ted;
extern double *d_Ked;

// Host pointers declaration
extern double *h_rho, *h_ux, *h_uy, *h_uz, *h_phi;
extern double *h_uyT;
extern double *h_f, *h_g;
extern double *h_ted, *h_Ked;

// Function declarations
void initialize();

__global__ void init_arrays(double *f, double *rho, double *ux, double *uy, double *uz, 
                            double *force_realx, double *force_realy, double *force_realz);

__global__ void collision_BGK(double *f, double *force_realx, double *force_realy, double *force_realz, 
                             double *rho, double *ux, double *uy, double *uz);

__global__ void streaming(double *f, double *f_temp);
__global__ void streaming2(double *f, double *f_temp);

__global__ void macrovar(double *f, double *force_realx, double *force_realy, double *force_realz, 
                        double *rho, double *ux, double *uy, double *uz);
__global__ void forcing(double *force_realx, double *force_realy, double *force_realz);
__global__ void forcing_GP(double *force_realx, double *force_realy, double *force_realz);

void diag_flow(int istep);
void output_flow(int istep);
void output_profile(int istep);
void output_fg(int istep);

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