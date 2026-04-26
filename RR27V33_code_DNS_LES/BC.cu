#include "LBMLES.h"
#include <math.h>

__global__ void RRBC(double *f, double *force_realx, double *force_realy, double *force_realz, 
                      double *rho, double *ux, double *uy, double *uz, double *taut) {
    int ix = blockIdx.x * blockDim.x + threadIdx.x;
    int iy = blockIdx.y * blockDim.y + threadIdx.y;
    int iz = blockIdx.z * blockDim.z + threadIdx.z;
    if (ix >= LX || iy >= LY || iz >= LZ) return;
    // if (ix < 1 || iy < 1 || iz < 1) return;       // 0: buffer layer
    // if (ix > LX-2 || iy > LY-2 || iz > LZ-2) return;
    if(iy != LY - 1 && iy != 0)return;

    int idx = (iz * LY + iy) * LX + ix;

    double Fbar[NPOP];
    double feqRR[NPOP];
    double forcefeq[NPOP];

    double fx9 = force_realx[idx];
    double fy9 = force_realy[idx];
    double fz9 = force_realz[idx];

    double rho9 = rho[idx];
    double u9   = ux[idx];
    double v9   = uy[idx];
    double w9   = uz[idx];
    double G3   = u9 * fx9 + v9 * fy9 + w9 * fz9;

    double uu = u9 / d_cs;
    double vv = v9 / d_cs;
    double ww = w9 / d_cs;
    double a110 = uu * vv, a101 = uu * ww, a011 = vv * ww, a200 = uu * uu, a020 = vv * vv, a002 = ww * ww;
    double a_110 = rho9 * uu * vv, a_101 = rho9 * uu * ww, a_011 = rho9 * vv * ww;
    double a_200 = rho9 * uu * uu, a_020 = rho9 * vv * vv, a_002 = rho9 * ww * ww;
    double a210 = a200 * vv, a201 = a200 * ww, a120 = a020 * uu, a021 = a020 * ww, a012 = a002 * vv, a102 = a002 * uu, a111 = a110 * ww;
    // double a220 = a210 * vv, a202 = a201 * ww, a022 = a021 * ww, a211 = a210 * ww, a121 = a120 * ww, a112 = a102 * vv;
    // double a221 = a220 * ww, a212 = a202 * vv, a122 = a022 * uu, a222 = a221 * ww;
    double ane110 = 0.0, ane101 = 0.0, ane011 = 0.0, ane200 = 0.0, ane020 = 0.0, ane002 = 0.0;

    for (int ip = 0; ip < NPOP; ip++) {
        double G1 = d_cix[ip] * fx9 + d_ciy[ip] * fy9 + d_ciz[ip] * fz9;
        Fbar[ip] = 0.5 * (G1 - G3) / d_RT;
        // feqRR & external body force
        double fir = d_H100[ip] * uu + d_H010[ip] * vv + d_H001[ip] * ww;
        double sec = 2.0 * (d_H110[ip] * a110 + d_H101[ip] * a101 + d_H011[ip] * a011)
                    + d_H200[ip] * a200 + d_H020[ip] * a020 + d_H002[ip] * a002;
        double thir = d_H210[ip] * a210 + d_H201[ip] * a201 + d_H021[ip] * a021 + d_H120[ip] * a120 +
                      d_H012[ip] * a012 + d_H102[ip] * a102 + 2.0 * d_H111[ip] * a111;
        // double four = d_H220[ip] * a220 + d_H202[ip] * a202 + d_H022[ip] * a022 + 
        //               2.0 * (d_H211[ip] * a211 + d_H121[ip] * a121 + d_H112[ip] * a112);
        // double fifsix = (d_H221[ip] * a221 + d_H212[ip] * a212 + d_H122[ip] * a122) / 4.0 + (d_H222[ip] * a222) / 8.0;
        // feqRR[ip] = d_tp[ip] * rho9 * (1.0 + fir + 0.5 * (sec + thir) + four / 4.0 + fifsix);
        // forcefeq[ip] = d_tp[ip] * rho9 * (1.0 + fir + 0.5 * (sec + thir) + four / 4.0);
        feqRR[ip] = d_tp[ip] * rho9 * (1.0 + fir + 0.5 * (sec + thir));
        forcefeq[ip] = d_tp[ip] * rho9 * (1.0 + fir + 0.5 * (sec + thir));
        // aneq
        ane110 += d_H110[ip] * (f[ip * LXYZ + idx] + Fbar[ip] * forcefeq[ip]);
        ane101 += d_H101[ip] * (f[ip * LXYZ + idx] + Fbar[ip] * forcefeq[ip]);
        ane011 += d_H011[ip] * (f[ip * LXYZ + idx] + Fbar[ip] * forcefeq[ip]);
        ane200 += d_H200[ip] * (f[ip * LXYZ + idx] + Fbar[ip] * forcefeq[ip]);
        ane020 += d_H020[ip] * (f[ip * LXYZ + idx] + Fbar[ip] * forcefeq[ip]);
        ane002 += d_H002[ip] * (f[ip * LXYZ + idx] + Fbar[ip] * forcefeq[ip]);
    }
        ane110 = ane110 - a_110;
        ane101 = ane101 - a_101;
        ane011 = ane011 - a_011;
        ane200 = ane200 - a_200;
        ane020 = ane020 - a_020;
        ane002 = ane002 - a_002;

        double ane210 = ane200 * vv + 2. * uu * ane110;      
        double ane201 = ane200 * ww + 2. * uu * ane101;    
        double ane120 = ane020 * uu + 2. * vv * ane110;
        double ane021 = ane020 * ww + 2. * vv * ane011;     
        double ane102 = ane002 * uu + 2. * ww * ane101;
        double ane012 = ane002 * vv + 2. * ww * ane011;
        double ane111 = ane110 * ww + ane101 * vv + ane011 * uu; // 3rd
        // double ane220 = ane210 * vv + uu*uu * ane020 + 2. * uu*vv * ane110;
        // double ane202 = ane201 * ww + uu*uu * ane002 + 2. * uu*ww * ane101;
        // double ane022 = ane021 * ww + vv*vv * ane002 + 2. * vv*ww * ane011;
        // double ane211 = ane210 * ww + uu*uu * ane011 + 2. * uu*vv * ane101;
        // double ane121 = ane120 * ww + 2. * uu*vv * ane011 + vv*vv * ane101;
        // double ane112 = ane111 * ww + uu*vv * ane002 + uu*ww * ane011 + vv*ww * ane101;   //4th   
        // double ane221 = ane220 * ww + 2. * a_200*vv * ane011 + 2. * a_020*uu * ane101;
        // double ane212 = ane211 * ww + 2. * a_110*ww * ane101 + a_200 * (vv*ane002 + ww*ane011);
        // double ane122 = ane121 * ww + 2. * a_110*ww * ane011 + a_020 * (uu*ane002 + ww*ane101);
        // double ane222 = ane221 * ww + a_200 * a_020 * ane002 + 2. * a_200 * a_011 * ane011+ 2.0 * a_101 * a_020 * ane101; // 5th & 6th
    if (iy == 0) {
        for (int ip = 0; ip < NPOP; ip++) {
            if(d_ciy[ip] > 0) {
                // fneqRR
                double secneq  = 2.0 * (d_H110[ip] * ane110 + d_H101[ip] * ane101 + d_H011[ip] * ane011)
                               + d_H200[ip] * ane200 + d_H020[ip] * ane020 + d_H002[ip] * ane002;
                double thirneq = d_H210[ip] * ane210 + d_H201[ip] * ane201 + d_H021[ip] * ane021 + d_H120[ip] * ane120 +
                                 d_H012[ip] * ane012 + d_H102[ip] * ane102 + 2.0 * d_H111[ip] * ane111;
                // double fourneq = d_H220[ip] * ane220 + d_H202[ip] * ane202 + d_H022[ip] * ane022 + 
                //           2.0 * (d_H211[ip] * ane211 + d_H121[ip] * ane121 + d_H112[ip] * ane112);
                // double fifsixneq = (d_H221[ip] * ane221 + d_H212[ip] * ane212 + d_H122[ip] * ane122) / 4.0 + (d_H222[ip] * ane222) / 8.0;
                double fneqRR = d_tp[ip] * (0.5 * (secneq + thirneq));
                // double fneqRR = d_tp[ip] * (0.5 * (secneq + thirneq) + fourneq / 4.0 + fifsixneq);

                f[ip * LXYZ + idx] = feqRR[ip] + fneqRR - Fbar[ip] * forcefeq[ip];
            }
        }
    }
    
    if (iy == LY - 1) {
        for (int ip = 0; ip < NPOP; ip++) {
            if(d_ciy[ip] < 0) {
                // fneqRR
                double secneq  = 2.0 * (d_H110[ip] * ane110 + d_H101[ip] * ane101 + d_H011[ip] * ane011)
                               + d_H200[ip] * ane200 + d_H020[ip] * ane020 + d_H002[ip] * ane002;
                double thirneq = d_H210[ip] * ane210 + d_H201[ip] * ane201 + d_H021[ip] * ane021 + d_H120[ip] * ane120 +
                                 d_H012[ip] * ane012 + d_H102[ip] * ane102 + 2.0 * d_H111[ip] * ane111;
                // double fourneq = d_H220[ip] * ane220 + d_H202[ip] * ane202 + d_H022[ip] * ane022 + 
                //           2.0 * (d_H211[ip] * ane211 + d_H121[ip] * ane121 + d_H112[ip] * ane112);
                // double fifsixneq = (d_H221[ip] * ane221 + d_H212[ip] * ane212 + d_H122[ip] * ane122) / 4.0 + (d_H222[ip] * ane222) / 8.0;
                double fneqRR = d_tp[ip] * (0.5 * (secneq + thirneq));
                // double fneqRR = d_tp[ip] * (0.5 * (secneq + thirneq) + fourneq / 4.0 + fifsixneq);

                f[ip * LXYZ + idx] = feqRR[ip] + fneqRR - Fbar[ip] * forcefeq[ip];
            }
        }
    }
}

// ============================== FDBC ===========================================
__global__ void FDBC(double *f, double *force_realx, double *force_realy, double *force_realz, 
                      double *rho, double *ux, double *uy, double *uz, double *taut) {
    int ix = blockIdx.x * blockDim.x + threadIdx.x;
    int iy = blockIdx.y * blockDim.y + threadIdx.y;
    int iz = blockIdx.z * blockDim.z + threadIdx.z;
    if (ix >= LX || iy >= LY || iz >= LZ) return;
    if(iy != LY - 1 && iy != 0)return;

    int idx = (iz * LY + iy) * LX + ix;

    double tau_eff;

    #if USE_LES
        tau_eff = d_tau + taut[idx];
    #else
        tau_eff = d_tau;
    #endif

    // compute velocity at the first node(y=0, y=1.5, y=2.5 -> y=0.5)
    if ( iy == 0 ) {
        int idxp1 =  (iz * LY + iy + 1) * LX + ix;
        int idxp2 =  (iz * LY + iy + 2) * LX + ix;
        ux[idx] = 2.0 / 3.0 * ux[idxp1] - 1.0 / 5.0 * ux[idxp2];
        uy[idx] = 2.0 / 3.0 * uy[idxp1] - 1.0 / 5.0 * uy[idxp2];
        uz[idx] = 2.0 / 3.0 * uz[idxp1] - 1.0 / 5.0 * uz[idxp2];
    } else if ( iy == LY - 1) {
        int idxm1 =  (iz * LY + iy - 1) * LX + ix;
        int idxm2 =  (iz * LY + iy - 2) * LX + ix;
        ux[idx] = 2.0 / 3.0 * ux[idxm1] - 1.0 / 5.0 * ux[idxm2];
        uy[idx] = 2.0 / 3.0 * uy[idxm1] - 1.0 / 5.0 * uy[idxm2];
        uz[idx] = 2.0 / 3.0 * uz[idxm1] - 1.0 / 5.0 * uz[idxm2];
    }

    // compute shear stress use second order center difference (y=0, y=1.5 -> y=0.5)
    double dudy, dvdy, dwdy;
    double dudx, dvdx, dwdx;
    double dudz, dvdz, dwdz;

    if ( iy == 0 ) {
        int idyp1 =  (iz * LY + iy + 1) * LX + ix;
        int idyp2 =  (iz * LY + iy + 2) * LX + ix;
        dudy = ( - 3.0 * ux[idx] + 4.0 * ux[idyp1] - ux[idyp2] ) / 2.0;
        dvdy = ( - 3.0 * uy[idx] + 4.0 * uy[idyp1] - uy[idyp2] ) / 2.0;
        dwdy = ( - 3.0 * uz[idx] + 4.0 * uz[idyp1] - uz[idyp2] ) / 2.0;
    } else if ( iy == LY - 1 ){
        int idym1 =  (iz * LY + iy - 1) * LX + ix;
        int idym2 =  (iz * LY + iy - 2) * LX + ix;
        dudy = ( 3.0 * ux[idx] - 4.0 * ux[idym1] + ux[idym2] ) / 2.0 ;
        dvdy = ( 3.0 * uy[idx] - 4.0 * uy[idym1] + uy[idym2] ) / 2.0 ;
        dwdy = ( 3.0 * uz[idx] - 4.0 * uz[idym1] + uz[idym2] ) / 2.0 ;
    } 
    
    int ixp = (ix + 1) % LX; 
    int ixm = (ix - 1 + LX) % LX;
    int izp = (iz + 1) % LZ; 
    int izm = (iz - 1 + LZ) % LZ;
    
    dudx = 0.5 * (ux[(iz * LY + iy) * LX + ixp] - ux[(iz * LY + iy) * LX + ixm]);
    dvdx = 0.5 * (uy[(iz * LY + iy) * LX + ixp] - uy[(iz * LY + iy) * LX + ixm]);
    dwdx = 0.5 * (uz[(iz * LY + iy) * LX + ixp] - uz[(iz * LY + iy) * LX + ixm]);
    
    dudz = 0.5 * (ux[(izp * LY + iy) * LX + ix] - ux[(izm * LY + iy) * LX + ix]);
    dvdz = 0.5 * (uy[(izp * LY + iy) * LX + ix] - uy[(izm * LY + iy) * LX + ix]);
    dwdz = 0.5 * (uz[(izp * LY + iy) * LX + ix] - uz[(izm * LY + iy) * LX + ix]);

    double S110 = 0.5 * (dudy + dvdx);
    double S101 = 0.5 * (dudz + dwdx);
    double S011 = 0.5 * (dvdz + dwdy);
    double S200 = dudx;
    double S020 = dvdy;
    double S002 = dwdz;

    // compute \rho
    double Fbar[NPOP];
    double feqRR[NPOP];
    double forcefeq[NPOP];
    double fx9 = force_realx[idx];
    double fy9 = force_realy[idx];
    double fz9 = force_realz[idx];

    double rho9 = 1.0;           // target rho = 1.0
    double u9   = ux[idx];
    double v9   = uy[idx];
    double w9   = uz[idx];
    double G3   = u9 * fx9 + v9 * fy9 + w9 * fz9;

    double uu = u9 / d_cs;
    double vv = v9 / d_cs;
    double ww = w9 / d_cs;
    double a110 = uu * vv, a101 = uu * ww, a011 = vv * ww, a200 = uu * uu, a020 = vv * vv, a002 = ww * ww;
    double a210 = a200 * vv, a201 = a200 * ww, a120 = a020 * uu, a021 = a020 * ww, a012 = a002 * vv, a102 = a002 * uu, a111 = a110 * ww;

    double ane110 = - rho9 * tau_eff * S110;
    double ane101 = - rho9 * tau_eff * S101;
    double ane011 = - rho9 * tau_eff * S011;
    double ane200 = - rho9 * tau_eff * S200;
    double ane020 = - rho9 * tau_eff * S020;
    double ane002 = - rho9 * tau_eff * S002;

    double ane210 = ane200 * vv + 2. * uu * ane110;      
    double ane201 = ane200 * ww + 2. * uu * ane101;    
    double ane120 = ane020 * uu + 2. * vv * ane110;
    double ane021 = ane020 * ww + 2. * vv * ane011;     
    double ane102 = ane002 * uu + 2. * ww * ane101;
    double ane012 = ane002 * vv + 2. * ww * ane011;
    double ane111 = ane110 * ww + ane101 * vv + ane011 * uu; // 3rd
    
    double sum_f_ps = 0.0;
    double sum_g_in = 0.0;

    for (int ip = 0; ip < NPOP; ip++) {
        bool known = ((iy == 0) ? (d_ciy[ip] < 0) : (d_ciy[ip] > 0));
        // known direction
        if(known) {
       // 1. f^ps
            double f_ps_temp = f[ip * LXYZ + idx];
            sum_f_ps += f_ps_temp;
        // 2. g^in = g^eq(1,u)+g^neq(1,S) 
            double G1 = d_cix[ip] * fx9 + d_ciy[ip] * fy9 + d_ciz[ip] * fz9;
            Fbar[ip] = 0.5 * (G1 - G3) / d_RT;
            // geqRR & external body force
            double fir = d_H100[ip] * uu + d_H010[ip] * vv + d_H001[ip] * ww;
            double sec = 2.0 * (d_H110[ip] * a110 + d_H101[ip] * a101 + d_H011[ip] * a011)
                        + d_H200[ip] * a200 + d_H020[ip] * a020 + d_H002[ip] * a002;
            double thir = d_H210[ip] * a210 + d_H201[ip] * a201 + d_H021[ip] * a021 + d_H120[ip] * a120 +
                          d_H012[ip] * a012 + d_H102[ip] * a102 + 2.0 * d_H111[ip] * a111;
            feqRR[ip] = d_tp[ip] * rho9 * (1.0 + fir + 0.5 * (sec + thir));
            forcefeq[ip] = d_tp[ip] * rho9 * (1.0 + fir + 0.5 * (sec + thir));
            // gneqRR
            double secneq  = 2.0 * (d_H110[ip] * ane110 + d_H101[ip] * ane101 + d_H011[ip] * ane011)
                           + d_H200[ip] * ane200 + d_H020[ip] * ane020 + d_H002[ip] * ane002;
            double thirneq = d_H210[ip] * ane210 + d_H201[ip] * ane201 + d_H021[ip] * ane021 + d_H120[ip] * ane120 +
                             d_H012[ip] * ane012 + d_H102[ip] * ane102 + 2.0 * d_H111[ip] * ane111;
            double fneqFD = d_tp[ip] * (secneq + thirneq);
            // g^in = g^eq + g^neq
            double g_in = feqRR[ip] + fneqFD - Fbar[ip] * forcefeq[ip];
            sum_g_in += g_in;
        }
    }

    // 3. compute rho9: sum(f^ps) / sum(g^in)
        double rho_local = sum_f_ps / sum_g_in;
        // printf("sum_f_ps=%e, sum_g_in=%e, rho=%e\n",sum_f_ps, sum_g_in, rho_local);

        double aneq110 = - rho_local * tau_eff * S110;
        double aneq101 = - rho_local * tau_eff * S101;
        double aneq011 = - rho_local * tau_eff * S011;
        double aneq200 = - rho_local * tau_eff * S200;
        double aneq020 = - rho_local * tau_eff * S020;
        double aneq002 = - rho_local * tau_eff * S002;

        double aneq210 = aneq200 * vv + 2. * uu * aneq110;      
        double aneq201 = aneq200 * ww + 2. * uu * aneq101;    
        double aneq120 = aneq020 * uu + 2. * vv * aneq110;
        double aneq021 = aneq020 * ww + 2. * vv * aneq011;     
        double aneq102 = aneq002 * uu + 2. * ww * aneq101;
        double aneq012 = aneq002 * vv + 2. * ww * aneq011;
        double aneq111 = aneq110 * ww + aneq101 * vv + aneq011 * uu; // 3rd

    // reconstruct distribution
    for (int ip = 0; ip < NPOP; ip++) {
        double G1 = d_cix[ip] * fx9 + d_ciy[ip] * fy9 + d_ciz[ip] * fz9;
        Fbar[ip] = 0.5 * (G1 - G3) / d_RT;
        // feqRR & external body force
        double fir = d_H100[ip] * uu + d_H010[ip] * vv + d_H001[ip] * ww;
        double sec = 2.0 * (d_H110[ip] * a110 + d_H101[ip] * a101 + d_H011[ip] * a011)
                    + d_H200[ip] * a200 + d_H020[ip] * a020 + d_H002[ip] * a002;
        double thir = d_H210[ip] * a210 + d_H201[ip] * a201 + d_H021[ip] * a021 + d_H120[ip] * a120 +
                      d_H012[ip] * a012 + d_H102[ip] * a102 + 2.0 * d_H111[ip] * a111;
        feqRR[ip] = d_tp[ip] * rho_local * (1.0 + fir + 0.5 * (sec + thir));
        forcefeq[ip] = d_tp[ip] * rho_local * (1.0 + fir + 0.5 * (sec + thir));
        // fneqRR
        double secneq  = 2.0 * (d_H110[ip] * aneq110 + d_H101[ip] * aneq101 + d_H011[ip] * aneq011)
                       + d_H200[ip] * aneq200 + d_H020[ip] * aneq020 + d_H002[ip] * aneq002;
        double thirneq = d_H210[ip] * aneq210 + d_H201[ip] * aneq201 + d_H021[ip] * aneq021 + d_H120[ip] * aneq120 +
                         d_H012[ip] * aneq012 + d_H102[ip] * aneq102 + 2.0 * d_H111[ip] * aneq111;
        double fneqFD = d_tp[ip] * (secneq + thirneq);

        f[ip * LXYZ + idx] = feqRR[ip] + fneqFD - Fbar[ip] * forcefeq[ip];
    }
}