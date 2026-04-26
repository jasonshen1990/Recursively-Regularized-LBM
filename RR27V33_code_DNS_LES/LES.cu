#include "LBMLES.h"
#include <math.h>

__global__ void LES_WALE(double *ux, double *uy, double *uz, double *taut) {
    int ix = blockIdx.x * blockDim.x + threadIdx.x;
    int iy = blockIdx.y * blockDim.y + threadIdx.y;
    int iz = blockIdx.z * blockDim.z + threadIdx.z;
    if (ix >= LX || iy >= LY || iz >= LZ) return;

    int idx = (iz * LY + iy) * LX + ix;
    
    double cw = 0.5;

    double dudy, dvdy, dwdy;
    double dudx, dvdx, dwdx;
    double dudz, dvdz, dwdz;

    // wall-normal direction
    if ( iy == 0 ) {
        int idxp =  (iz * LY + iy + 1) * LX + ix;
        dudy = (3.0 * ux[idx] + ux[idxp] ) / 3.0 ;
        dvdy = (3.0 * uy[idx] + uy[idxp] ) / 3.0 ;
        dwdy = (3.0 * uz[idx] + uz[idxp] ) / 3.0 ;
    } else if ( iy == LY - 1 ){
        int idxm =  (iz * LY + iy - 1) * LX + ix;
        dudy = ( -3.0 * ux[idx] - ux[idxm] ) / 3.0 ;
        dvdy = ( -3.0 * uy[idx] - uy[idxm] ) / 3.0 ;
        dwdy = ( -3.0 * uz[idx] - uz[idxm] ) / 3.0 ;
    } else {
        int idxp =  (iz * LY + iy + 1) * LX + ix;
        int idxm =  (iz * LY + iy - 1) * LX + ix;
        dudy = 0.5 * ( ux[idxp] - ux[idxm] );
        dvdy = 0.5 * ( uy[idxp] - uy[idxm] );
        dwdy = 0.5 * ( uz[idxp] - uz[idxm] );
    }
    
    // streamwise direction: central difference
    if ( ix == 0 ) {
        int idxp =  (iz * LY + iy) * LX + ix + 1;
        int idxm =  (iz * LY + iy) * LX + LX - 1;
        dudx = 0.5 * ( ux[idxp] - ux[idxm] );
        dvdx = 0.5 * ( uy[idxp] - uy[idxm] );
        dwdx = 0.5 * ( uz[idxp] - uz[idxm] );
    } else if ( ix == LX - 1 ){
        int idxp =  (iz * LY + iy) * LX + 0;
        int idxm =  (iz * LY + iy) * LX + ix - 1;
        dudx = 0.5 * ( ux[idxp] - ux[idxm] );
        dvdx = 0.5 * ( uy[idxp] - uy[idxm] );
        dwdx = 0.5 * ( uz[idxp] - uz[idxm] );
    } else {
        int idxp =  (iz * LY + iy) * LX + ix + 1;
        int idxm =  (iz * LY + iy) * LX + ix - 1;
        dudx = 0.5 * ( ux[idxp] - ux[idxm] );
        dvdx = 0.5 * ( uy[idxp] - uy[idxm] );
        dwdx = 0.5 * ( uz[idxp] - uz[idxm] );
    }

    // spanwise direction: central difference
    if ( iz == 0 ) {
        int idxp = (( iz + 1 ) * LY + iy) * LX + ix;
        int idxm = (( LZ - 1 ) * LY + iy) * LX + ix;
        dudz = 0.5 * ( ux[idxp] - ux[idxm] );
        dvdz = 0.5 * ( uy[idxp] - uy[idxm] );
        dwdz = 0.5 * ( uz[idxp] - uz[idxm] );
    } else if ( iz == LZ - 1 ){
        int idxp = ( 0 * LY + iy) * LX + ix;
        int idxm = (( iz - 1 ) * LY + iy) * LX + ix;
        dudz = 0.5 * ( ux[idxp] - ux[idxm] );
        dvdz = 0.5 * ( uy[idxp] - uy[idxm] );
        dwdz = 0.5 * ( uz[idxp] - uz[idxm] );
    } else {
        int idxp = (( iz + 1 ) * LY + iy) * LX + ix;
        int idxm = (( iz - 1 ) * LY + iy) * LX + ix;
        dudz = 0.5 * ( ux[idxp] - ux[idxm] );
        dvdz = 0.5 * ( uy[idxp] - uy[idxm] );
        dwdz = 0.5 * ( uz[idxp] - uz[idxm] );
    }

    double sijsq = (dudy + dvdx)*(dudy + dvdx) + (dudz + dwdx)*(dudz + dwdx) + (dvdz + dwdy)*(dvdz + dwdy);
    sijsq = 0.5*sijsq + dudx * dudx + dvdy * dvdy + dwdz * dwdz;

    double s11d = 2.* dudx * dudx - dvdy * dvdy - dwdz * dwdz + dudy * dvdx + dudz * dwdx - 2.* dvdz * dwdy;
    s11d = s11d / 3.;
    double s12d = dudx * (dudy + dvdx) + dudy * dvdy + dvdx * dvdy + dvdz * dwdx + dudz * dwdy;
    s12d = 0.5 * s12d;
    double s13d = dudx * (dudz + dwdx) + dudy * dvdz + dvdx * dwdy + dudz * dwdz + dwdx * dwdz;
    s13d = 0.5 * s13d;
    double s22d = 2.* dvdy * dvdy - dudx * dudx - dwdz * dwdz + dudy * dvdx + dvdz * dwdy - 2.* dudz * dwdx;
    s22d = s22d / 3.;
    double s23d = dvdy * (dvdz + dwdy) + dudz * dvdx + dudy * dwdx + dvdz * dwdz + dwdy * dwdz;
    s23d = 0.5 * s23d;
    double s33d = 2.* dwdz * dwdz - dudx * dudx - dvdy * dvdy + dudz * dwdx + dvdz * dwdy - 2.* dudy * dvdx;
    s33d = s33d / 3.;

    double sijdsq = s11d * s11d + s22d * s22d + s33d * s33d + 2.0 *(s12d * s12d + s13d * s13d + s23d * s23d);

    taut[idx] = 3.0 * cw * cw * pow(sijdsq, 1.5) / ( pow(sijsq, 2.5) + pow(sijdsq, 1.25) );
}
