#include "LBMLES.h"
#include "parameters.h"

__constant__ int d_cix[NPOP];
__constant__ int d_ciy[NPOP];
__constant__ int d_ciz[NPOP];
__constant__ double d_tp[NPOP];
__constant__ double d_tau;
__constant__ double d_visc;
__constant__ double d_ustar, d_ystar;
__constant__ double d_force_in_x;
__constant__ double d_pi, d_pi2;
__constant__ double d_RT, d_cs;
__constant__ double d_sigma;
__device__ int d_istep;
__device__ double d_rand_phase;

double pi = acos(-1.0);
double pi2 = 2 * pi;
double RT = 1.0 / 3.0;
double cs = sqrt(RT);

//RR_Hermite coefficients
__constant__ double d_H100[NPOP], d_H010[NPOP], d_H001[NPOP];
__constant__ double d_H110[NPOP], d_H101[NPOP], d_H011[NPOP], d_H200[NPOP], d_H020[NPOP], d_H002[NPOP];
__constant__ double d_H210[NPOP], d_H201[NPOP], d_H021[NPOP], d_H120[NPOP], d_H102[NPOP], d_H012[NPOP], d_H111[NPOP];
// __constant__ double d_H220[NPOP], d_H202[NPOP], d_H022[NPOP], d_H211[NPOP], d_H121[NPOP], d_H112[NPOP];
// __constant__ double d_H221[NPOP], d_H212[NPOP], d_H122[NPOP], d_H222[NPOP];

double sigma = 1;

double Re = 180.0;
double visc = 0.000622;
double u0 = 2.0 * Re * visc / LY; 
double tau = 3.0 * visc + 0.5;
double ystar = visc / u0;

double force_in_x = 2.0 * u0 * u0 / LY;

int continue_step = CONTINUE_STEP;

// Device pointers definition
double *d_f = nullptr, *d_rho = nullptr, *d_ux = nullptr, *d_uy = nullptr, *d_uz = nullptr;
double *d_force_realx = nullptr, *d_force_realy = nullptr, *d_force_realz = nullptr;
double *d_f_temp = nullptr;
double *d_taut = nullptr; //LES turbulent viscosity  taut = nu_t / (c_s^2 * dt) = 3 * nu_t

// Host pointers definition
double *h_rho = nullptr, *h_ux = nullptr, *h_uy = nullptr, *h_uz = nullptr;
double *h_f = nullptr;
double *h_taut = nullptr; //LES turbulent viscosity  taut = nu_t / (c_s^2 * dt) = 3 * nu_t