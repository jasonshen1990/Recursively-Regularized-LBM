#include "LBMLES.h"

__global__ void streaming(double *f, double *f_temp) {
    int ix = blockIdx.x * blockDim.x + threadIdx.x;
    int iy = blockIdx.y * blockDim.y + threadIdx.y;
    int iz = blockIdx.z * blockDim.z + threadIdx.z;
    if (ix >= LX || iy >= LY || iz >= LZ) return;

    int idx = (iz * LY + iy) * LX + ix;

    f_temp[idx] = f[idx];

    for (int ip = 0; ip < NPOP; ip++) {
        int imove = ix + d_cix[ip];
        int jmove = iy + d_ciy[ip];
        int kmove = iz + d_ciz[ip];
        imove = (imove + LX) % LX;
        jmove = (jmove + LY) % LY;
        kmove = (kmove + LZ) % LZ;
        int idxnew = (kmove * LY + jmove) * LX + imove;
        f_temp[ip * LXYZ + idxnew] = f[ip * LXYZ + idx];
    }

}

__global__ void BB(double *f, double *f_temp) {
    int ix = blockIdx.x * blockDim.x + threadIdx.x;
    int iy = blockIdx.y * blockDim.y + threadIdx.y;
    int iz = blockIdx.z * blockDim.z + threadIdx.z;
    if (ix >= LX || iy >= LY || iz >= LZ) return;

    int idx = (iz * LY + iy) * LX + ix;

    if (iy == 0) {  // Bottom wall
        f_temp[2  * LXYZ + idx] = f[4  * LXYZ + idx];
        f_temp[7  * LXYZ + idx] = f[9  * LXYZ + idx];
        f_temp[8  * LXYZ + idx] = f[10 * LXYZ + idx];
        f_temp[12 * LXYZ + idx] = f[18 * LXYZ + idx];
        f_temp[19 * LXYZ + idx] = f[25 * LXYZ + idx];
        f_temp[20 * LXYZ + idx] = f[26 * LXYZ + idx];
        f_temp[16 * LXYZ + idx] = f[14 * LXYZ + idx];
        f_temp[23 * LXYZ + idx] = f[21 * LXYZ + idx];
        f_temp[24 * LXYZ + idx] = f[22 * LXYZ + idx];
    }
    if (iy == LY-1) {  // Top wall
        f_temp[4  * LXYZ + idx] = f[2  * LXYZ + idx];
        f_temp[9  * LXYZ + idx] = f[7  * LXYZ + idx];
        f_temp[10 * LXYZ + idx] = f[8  * LXYZ + idx];
        f_temp[18 * LXYZ + idx] = f[12 * LXYZ + idx];
        f_temp[25 * LXYZ + idx] = f[19 * LXYZ + idx];
        f_temp[26 * LXYZ + idx] = f[20 * LXYZ + idx];
        f_temp[14 * LXYZ + idx] = f[16 * LXYZ + idx];
        f_temp[21 * LXYZ + idx] = f[23 * LXYZ + idx];
        f_temp[22 * LXYZ + idx] = f[24 * LXYZ + idx];
    }
}