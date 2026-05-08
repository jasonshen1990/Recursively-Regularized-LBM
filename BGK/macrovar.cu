#include "rb3d.h"
#include <math.h>

__global__ void macrovar(double *f, double *force_realx, double *force_realy, double *force_realz, 
                        double *rho, double *ux, double *uy, double *uz) {
    int ix = blockIdx.x * blockDim.x + threadIdx.x;
    int iy = blockIdx.y * blockDim.y + threadIdx.y;
    int iz = blockIdx.z * blockDim.z + threadIdx.z;
    if (ix >= LX || iy >= LY || iz >= LZ) return;

    int idx = (iz * LY + iy) * LX + ix;
    double rho_local = 0.0;
    double ux_local = 0.0;
    double uy_local = 0.0;
    double uz_local = 0.0;

    for (int ip = 0; ip < NPOP; ip++) {
        double f_val = f[ip * LXYZ + idx];

        rho_local += f_val;
        ux_local  += d_cix[ip] * f_val;
        uy_local  += d_ciy[ip] * f_val;
        uz_local  += d_ciz[ip] * f_val;
    }

    rho[idx] = rho_local;
    ux[idx]  = ux_local + 0.5 * force_realx[idx]; 
    uy[idx]  = uy_local + 0.5 * force_realy[idx];
    uz[idx]  = uz_local + 0.5 * force_realz[idx];

}

__global__ void forcing(double *force_realx, double *force_realy, double *force_realz) {
    int ix = blockIdx.x * blockDim.x + threadIdx.x;
    int iy = blockIdx.y * blockDim.y + threadIdx.y;
    int iz = blockIdx.z * blockDim.z + threadIdx.z;
    if (ix >= LX || iy >= LY || iz >= LZ) return;

    int idx = (iz * LY + iy) * LX + ix;

    force_realx[idx] = d_force_in_x;
    force_realy[idx] = 0.0;
    force_realz[idx] = 0.0;  
}

__global__ void forcing_GP(double *force_realx, double *force_realy, double *force_realz) {
    int ix = blockIdx.x * blockDim.x + threadIdx.x;
    int iy = blockIdx.y * blockDim.y + threadIdx.y;
    int iz = blockIdx.z * blockDim.z + threadIdx.z;
    if (ix >= LX || iy >= LY || iz >= LZ) return;

    int idx = (iz * LY + iy) * LX + ix;

    double Tpd = 2000.0;
    //int Tpdp = 1500.0;
    // double LYh = LY / 2;
    double LYq = LY / 4;
    int h0 = 2;
    double beta9 = 3.0;
    double gamma9 = 2.0;
    // double phase9 = 0.0;   //0.25x   0.15x
    double Amp0 = ( 100.0 * beta9 / LY ) * sin(d_pi2 * d_istep / Tpd);

    force_realx[idx] = d_force_in_x;
    force_realy[idx] = 0.0;
    force_realz[idx] = 0.0;

    if (iy > h0 && iy < LYq + h0) {
        force_realx[idx] = d_force_in_x * (1.0 - Amp0 * LX / beta9 * sin(d_pi2 * ( iy - h0 ) / LYq) *
                           sin(beta9 * d_pi2 * ix /LX) * cos(gamma9 * d_pi2 * iz / LZ));

        force_realy[idx] = d_force_in_x * Amp0 * 0.5 * LYq * (1 - cos(d_pi2 * ( iy - h0 ) / LYq)) *
                           cos(beta9 * d_pi2 * ix /LX) * cos(gamma9 * d_pi2 * iz / LZ);

        force_realz[idx] = d_force_in_x * Amp0 * 0.5 * LZ / gamma9 * sin(d_pi2 * ( iy - h0 ) / LYq) *
                           cos(beta9 * d_pi2 * ix /LX) * sin(gamma9 * d_pi2 * iz / LZ); 
    } else if(iy > 0.75 * LY - h0 && iy < LY - h0) {
        force_realx[idx] = d_force_in_x * (1.0 + Amp0 * LX / beta9 * sin(d_pi2 * (-LY + iy + h0 + LYq ) / LYq) *
                           sin(beta9 * d_pi2 * (ix /LX + d_rand_phase)) * cos(gamma9 * d_pi2 * iz / LZ));

        force_realy[idx] = -d_force_in_x * Amp0 * 0.5 * LYq * (1 - cos(d_pi2 * (-LY + iy + h0 + LYq ) / LYq)) *
                           cos(beta9 * d_pi2 * (ix /LX + d_rand_phase)) * cos(gamma9 * d_pi2 * iz / LZ);

        force_realz[idx] = -d_force_in_x * Amp0 * 0.5 * LZ / gamma9 * sin(d_pi2 * (-LY + iy + h0 + LYq ) / LYq) *
                           cos(beta9 * d_pi2 * (ix /LX + d_rand_phase)) * sin(gamma9 * d_pi2 * iz / LZ);
    }

}