#include "LBMLES.h"
#include <math.h>

__global__ void collision_RR(double *f, double *force_realx, double *force_realy, double *force_realz, 
                             double *rho, double *ux, double *uy, double *uz, double *taut) {
    int ix = blockIdx.x * blockDim.x + threadIdx.x;
    int iy = blockIdx.y * blockDim.y + threadIdx.y;
    int iz = blockIdx.z * blockDim.z + threadIdx.z;
    if (ix >= LX || iy >= LY || iz >= LZ) return;

    int idx = (iz * LY + iy) * LX + ix;

    double Fbar[NPOP];
    double feqRR[NPOP];
    double forcefeq[NPOP];

    double fx9  = force_realx[idx];
    double fy9  = force_realy[idx];
    double fz9  = force_realz[idx];
    double rho9 = rho[idx];
    double u9   = ux[idx];
    double v9   = uy[idx];
    double w9   = uz[idx];
    double G3   = u9 * fx9 + v9 * fy9 + w9 * fz9;
    double tau_eff;

    #if USE_LES
        tau_eff = d_tau + taut[idx];
    #else
        tau_eff = d_tau;
    #endif

    if(iy == LY/2-1 && ix == 5 && iz == 5 && d_istep % 1000 == 0){
        #if USE_LES
            printf("LES: taut[idx]=%.6e, tau_eff=%.6e\n", taut[idx], tau_eff);
        #else
            printf("DNS: tau_eff=%.6e\n", tau_eff);
        #endif
    }

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
        // forcefeq[ip] = d_tp[ip]  * rho9 * (1.0 + fir + 0.5 * (sec + thir) + four / 4.0);

        feqRR[ip] = d_tp[ip] * rho9 * (1.0 + fir + 0.5 * (sec + thir));
        forcefeq[ip] = d_tp[ip]  * rho9 * (1.0 + fir + 0.5 * (sec + thir));
    
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
        // double ane222 = ane221 * ww + a_200 * a_020 * ane002 + 2. * a_200 * a_011 * ane011+ 2.0 * a_101 * a_020 * ane101; // 5th & 6th // 5th & 6th

    for (int ip = 0; ip < NPOP; ip++) {
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

        f[ip * LXYZ + idx] = feqRR[ip] + (1.0 - 1.0 / tau_eff ) * fneqRR + Fbar[ip] * forcefeq[ip];
    }
}