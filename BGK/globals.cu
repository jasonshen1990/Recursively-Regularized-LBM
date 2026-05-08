#include "rb3d.h"
#include "parameters.h"

__constant__ int d_cix[NPOP];
__constant__ int d_ciy[NPOP];
__constant__ int d_ciz[NPOP];
__constant__ double d_tp[NPOP];
__constant__ double d_tau, d_tauc, d_grav0, d_beta, d_tHot, d_tCold, d_t0;
__constant__ double d_diff, d_visc;
__constant__ double d_ustar, d_ystar;
__constant__ double d_force_in_x;
__constant__ double d_pi, d_pi2;
__device__ int d_istep;
__device__ double d_rand_phase;

double pi = acos(-1.0);
double pi2 = 2 * pi;

double Re = 180.0;
double visc = 0.0036;
double u0 = 2.0 * Re * visc / LY; 
double tau = 3.0 * visc + 0.5;
double ystar = visc / u0;

double grav0 = GRAVITY;
double beta = BETA;
double t0 = 0.5; // used to be the mean temperature* (tCold + tHot);

// double force_in_x = 8.0 * visc * u0 / LY / LY;
double force_in_x = 2.0 * u0 * u0 / LY;

int continue_step = CONTINUE_STEP;

// Device pointers definition
double *d_f = nullptr, *d_rho = nullptr, *d_ux = nullptr, *d_uy = nullptr, *d_uz = nullptr;
double *d_force_realx = nullptr, *d_force_realy = nullptr, *d_force_realz = nullptr;
double *d_f_temp = nullptr, *d_g_temp = nullptr;
double *d_uyT = nullptr;
double *d_Tx = nullptr, *d_Ty = nullptr, *d_Tz = nullptr, *d_ted = nullptr;
double *d_Ked = nullptr;

// Host pointers definition
double *h_rho = nullptr, *h_ux = nullptr, *h_uy = nullptr, *h_uz = nullptr, *h_phi = nullptr;
double *h_uyT = nullptr;
double *h_f = nullptr, *h_g = nullptr;
double *h_ted = nullptr, *h_Ked = nullptr;