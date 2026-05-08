#include "rb3d.h"

__global__ void collision_BGK(double *f, double *force_realx, double *force_realy, double *force_realz, 
                             double *rho, double *ux, double *uy, double *uz) {
    int ix = blockIdx.x * blockDim.x + threadIdx.x;
    int iy = blockIdx.y * blockDim.y + threadIdx.y;
    int iz = blockIdx.z * blockDim.z + threadIdx.z;
    if (ix >= LX || iy >= LY || iz >= LZ) return;

    int idx = (iz * LY + iy) * LX + ix;
    double fx9  = force_realx[idx];
    double fy9  = force_realy[idx];
    double fz9  = force_realz[idx];

    double rho9 = rho[idx];
    double u9   = ux[idx];
    double v9   = uy[idx];
    double w9   = uz[idx];

    double G3 = u9 * fx9 + v9 * fy9 + w9 * fz9;

    for (int ip = 0; ip < NPOP; ip++) {
        double RT = 1.0 / 3.0; 
        double eu = (d_cix[ip] * u9 + d_ciy[ip] * v9 + d_ciz[ip] * w9) / RT;
        double uv = (u9 * u9 + v9 * v9 + w9 * w9) / RT;
        double feq = d_tp[ip] * (rho9 + eu + 0.5 * (eu * eu - uv));
        double G1 = d_cix[ip] * fx9 + d_ciy[ip] * fy9 + d_ciz[ip] * fz9;

        double fchange = -(f[ip * LXYZ + idx] - feq) / d_tau;
        double feq_rho0 = d_tp[ip] * (1.0 + eu + 0.5 * (eu * eu - uv));
        fchange += (1.0 - 0.5 / d_tau) * (G1-G3) / RT * feq_rho0;
        f[ip * LXYZ + idx] += fchange;
    }
}