#include "LBMLES.h"
#include <math.h>

void initialize() {
    // Allocate device memory
    CHECK_CUDA_ERROR(cudaMalloc(&d_f,      NPOP * LXYZ * sizeof(double)));
    CHECK_CUDA_ERROR(cudaMalloc(&d_f_temp, NPOP * LXYZ * sizeof(double)));

    CHECK_CUDA_ERROR(cudaMalloc(&d_rho,           LXYZ * sizeof(double)));
    CHECK_CUDA_ERROR(cudaMalloc(&d_ux,            LXYZ * sizeof(double)));
    CHECK_CUDA_ERROR(cudaMalloc(&d_uy,            LXYZ * sizeof(double)));
    CHECK_CUDA_ERROR(cudaMalloc(&d_uz,            LXYZ * sizeof(double)));

    CHECK_CUDA_ERROR(cudaMalloc(&d_taut,          LXYZ * sizeof(double)));    // LES

    CHECK_CUDA_ERROR(cudaMalloc(&d_force_realx,   LXYZ * sizeof(double)));
    CHECK_CUDA_ERROR(cudaMalloc(&d_force_realy,   LXYZ * sizeof(double)));
    CHECK_CUDA_ERROR(cudaMalloc(&d_force_realz,   LXYZ * sizeof(double)));


    // Allocate host memory with error checking
    h_rho = (double*)malloc(LXYZ * sizeof(double));
    if (h_rho == NULL) {
        fprintf(stderr, "Failed to allocate host memory for h_rho\n");
        exit(EXIT_FAILURE);
    }
    
    h_ux = (double*)malloc(LXYZ * sizeof(double));
    if (h_ux == NULL) {
        fprintf(stderr, "Failed to allocate host memory for h_ux\n");
        free(h_rho);
        exit(EXIT_FAILURE);
    }
    
    h_uy = (double*)malloc(LXYZ * sizeof(double));
    if (h_uy == NULL) {
        fprintf(stderr, "Failed to allocate host memory for h_uy\n");
        free(h_rho); free(h_ux);
        exit(EXIT_FAILURE);
    }
    
    h_uz = (double*)malloc(LXYZ * sizeof(double));
    if (h_uz == NULL) {
        fprintf(stderr, "Failed to allocate host memory for h_uz\n");
        free(h_rho); free(h_ux); free(h_uy);
        exit(EXIT_FAILURE);
    }
    
    
    h_f = (double*)malloc(NPOP * LXYZ * sizeof(double));
    if (h_f == NULL) {
        fprintf(stderr, "Failed to allocate host memory for h_f\n");
        free(h_rho); free(h_ux); free(h_uy); free(h_uz);
        exit(EXIT_FAILURE);
    }
    

    h_taut = (double*)malloc(LXYZ * sizeof(double));
    if (h_taut == NULL) {
        fprintf(stderr, "Failed to allocate host memory for h_phi\n");
        free(h_rho); free(h_ux); free(h_uy); free(h_uz); free(h_f);
        exit(EXIT_FAILURE);
    }

    // Initialize constant memory
    //                 0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26 
    int h_cix[NPOP] = {0,  1,  0, -1,  0,  0,  0,  1, -1, -1,  1,  1,  0, -1,  0,  1,  0, -1,  0,  1, -1, -1,  1,  1, -1, -1,  1};
    int h_ciy[NPOP] = {0,  0,  1,  0, -1,  0,  0,  1,  1, -1, -1,  0,  1,  0, -1,  0,  1,  0, -1,  1,  1, -1, -1,  1,  1, -1, -1};
    int h_ciz[NPOP] = {0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1,  1,  1,  1, -1, -1, -1, -1,  1,  1,  1,  1, -1, -1, -1, -1};
     
    double h_tp[NPOP];
    for (int i = 0; i < NPOP; i++) {
        if      (i == 0)  h_tp[i] = 8.0/27.0;
        else if (i <= 6)  h_tp[i] = 2.0/27.0;
        else if (i <= 18) h_tp[i] = 1.0/54.0;
        else              h_tp[i] = 1.0/216.0;
    }

    double RT = 1.0 / 3.0;
    double cs = sqrt(RT);
    double cs3 = RT * cs;
    double cs4 = RT * RT;
    double cs5 = RT * RT * cs;
    double cs6 = RT * RT * RT;

    double h_H100[NPOP], h_H010[NPOP], h_H001[NPOP];
    double h_H110[NPOP], h_H101[NPOP], h_H011[NPOP], h_H200[NPOP], h_H020[NPOP], h_H002[NPOP];
    double h_H210[NPOP], h_H201[NPOP], h_H021[NPOP], h_H120[NPOP], h_H102[NPOP], h_H012[NPOP], h_H111[NPOP];
    // double h_H220[NPOP], h_H202[NPOP], h_H022[NPOP], h_H211[NPOP], h_H121[NPOP], h_H112[NPOP];
    // double h_H221[NPOP], h_H212[NPOP], h_H122[NPOP], h_H222[NPOP];
    // double H000 = 1.0;
    for (int i = 0; i < NPOP; i++) {
        h_H100[i] = h_cix[i] / cs;
        h_H010[i] = h_ciy[i] / cs;
        h_H001[i] = h_ciz[i] / cs;
        // second order
        h_H110[i] = h_cix[i] * h_ciy[i] / RT;
        h_H101[i] = h_cix[i] * h_ciz[i] / RT;
        h_H011[i] = h_ciy[i] * h_ciz[i] / RT;
        h_H200[i] = h_cix[i] * h_cix[i] / RT - 1;
        h_H020[i] = h_ciy[i] * h_ciy[i] / RT - 1;
        h_H002[i] = h_ciz[i] * h_ciz[i] / RT - 1;
        // third oeder
        h_H210[i] = h_cix[i] * h_cix[i] * h_ciy[i] / cs3 - h_ciy[i] / cs;
        h_H201[i] = h_cix[i] * h_cix[i] * h_ciz[i] / cs3 - h_ciz[i] / cs;
        h_H021[i] = h_ciy[i] * h_ciy[i] * h_ciz[i] / cs3 - h_ciz[i] / cs;
        h_H120[i] = h_cix[i] * h_ciy[i] * h_ciy[i] / cs3 - h_cix[i] / cs;
        h_H102[i] = h_cix[i] * h_ciz[i] * h_ciz[i] / cs3 - h_cix[i] / cs;
        h_H012[i] = h_ciy[i] * h_ciz[i] * h_ciz[i] / cs3 - h_ciy[i] / cs;
        h_H111[i] = h_cix[i] * h_ciy[i] * h_ciz[i] / cs3;
        // forth order
        // h_H220[i] = h_cix[i] * h_cix[i] * h_ciy[i] * h_ciy[i] / cs4 - 
        //             h_cix[i] * h_cix[i] / RT - h_ciy[i] * h_ciy[i] / RT + 1;
        // h_H202[i] = h_cix[i] * h_cix[i] * h_ciz[i] * h_ciz[i] / cs4 - 
        //             h_cix[i] * h_cix[i] / RT - h_ciz[i] * h_ciz[i] / RT + 1;
        // h_H022[i] = h_ciy[i] * h_ciy[i] * h_ciz[i] * h_ciz[i] / cs4 - 
        //             h_ciy[i] * h_ciy[i] / RT - h_ciz[i] * h_ciz[i] / RT + 1;
        // h_H211[i] = h_cix[i] * h_cix[i] * h_ciy[i] * h_ciz[i] / cs4 - h_ciy[i] * h_ciz[i] / RT;
        // h_H121[i] = h_cix[i] * h_ciy[i] * h_ciy[i] * h_ciz[i] / cs4 - h_cix[i] * h_ciz[i] / RT;
        // h_H112[i] = h_cix[i] * h_ciy[i] * h_ciz[i] * h_ciz[i] / cs4 - h_cix[i] * h_ciy[i] / RT;
        // // fifth & sixth order
        // h_H221[i] = h_cix[i] * h_cix[i] * h_ciy[i] * h_ciy[i] * h_ciz[i]/ cs5 - 
        //             h_cix[i] * h_cix[i] * h_ciz[i] / cs3 - h_ciy[i] * h_ciy[i] * h_ciz[i]/ cs3 + h_ciz[i] / cs;
        // h_H212[i] = h_cix[i] * h_cix[i] * h_ciy[i] * h_ciz[i] * h_ciz[i]/ cs5 - 
        //             h_cix[i] * h_cix[i] * h_ciy[i] / cs3 - h_ciy[i] * h_ciz[i] * h_ciz[i]/ cs3 + h_ciy[i] / cs;
        // h_H122[i] = h_cix[i] * h_ciy[i] * h_ciy[i] * h_ciz[i] * h_ciz[i]/ cs5 - 
        //             h_cix[i] * h_ciy[i] * h_ciy[i] / cs3 - h_cix[i] * h_ciz[i] * h_ciz[i]/ cs3 + h_cix[i] / cs;
        // h_H222[i] = h_cix[i] * h_cix[i] * h_ciy[i] * h_ciy[i] * h_ciz[i] * h_ciz[i]/ cs6 - 
        //             h_cix[i] * h_cix[i] * h_ciy[i] * h_ciy[i] / cs4 - h_cix[i] * h_cix[i] * h_ciz[i] * h_ciz[i]/ cs4 + 
        //             h_ciy[i] * h_ciy[i] * h_ciz[i] * h_ciz[i] / cs4 + h_cix[i] * h_cix[i] / RT + h_ciy[i] * h_ciy[i] / RT + h_ciz[i] * h_ciz[i] / RT - 1;
    }

    // Copy parameters to constant memory
    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_cix,   h_cix,  NPOP * sizeof(int)));
    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_ciy,   h_ciy,  NPOP * sizeof(int)));
    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_ciz,   h_ciz,  NPOP * sizeof(int)));
    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_tp,    h_tp,   NPOP * sizeof(double)));
 
    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_H100,  h_H100, NPOP * sizeof(double)));
    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_H010,  h_H010, NPOP * sizeof(double)));
    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_H001,  h_H001, NPOP * sizeof(double)));  // 1st
    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_H200,  h_H200, NPOP * sizeof(double)));
    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_H020,  h_H020, NPOP * sizeof(double)));
    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_H002,  h_H002, NPOP * sizeof(double)));
    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_H110,  h_H110, NPOP * sizeof(double)));
    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_H101,  h_H101, NPOP * sizeof(double)));
    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_H011,  h_H011, NPOP * sizeof(double)));  // 2nd
    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_H210,  h_H210, NPOP * sizeof(double)));
    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_H201,  h_H201, NPOP * sizeof(double)));
    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_H021,  h_H021, NPOP * sizeof(double)));
    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_H120,  h_H120, NPOP * sizeof(double)));
    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_H102,  h_H102, NPOP * sizeof(double)));
    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_H012,  h_H012, NPOP * sizeof(double)));
    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_H111,  h_H111, NPOP * sizeof(double)));  // 3rd
    // CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_H220,  h_H220, NPOP * sizeof(double)));
    // CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_H202,  h_H202, NPOP * sizeof(double)));
    // CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_H022,  h_H022, NPOP * sizeof(double)));
    // CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_H211,  h_H211, NPOP * sizeof(double)));
    // CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_H121,  h_H121, NPOP * sizeof(double)));
    // CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_H112,  h_H112, NPOP * sizeof(double)));  // 4th
    // CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_H221,  h_H221, NPOP * sizeof(double)));
    // CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_H212,  h_H212, NPOP * sizeof(double)));
    // CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_H122,  h_H122, NPOP * sizeof(double)));  // 5th
    // CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_H222,  h_H222, NPOP * sizeof(double)));  // 6th

    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_tau,    &tau,   sizeof(double)));
    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_ustar,  &u0,    sizeof(double)));
    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_ystar,  &ystar, sizeof(double)));

    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_force_in_x,  &force_in_x,        sizeof(double)));

    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_pi,    &pi,          sizeof(double)));
    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_pi2,   &pi2,         sizeof(double)));
    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_RT,    &RT,          sizeof(double)));
    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_cs,    &cs,          sizeof(double)));
    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_sigma, &sigma,       sizeof(double)));

    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_visc,  &visc,        sizeof(double)));

    // Initialize device arrays
    dim3 grid((LX + BLOCK_X - 1) / BLOCK_X, (LY + BLOCK_Y - 1) / BLOCK_Y, (LZ + BLOCK_Z - 1) / BLOCK_Z);
    dim3 block(BLOCK_X, BLOCK_Y, BLOCK_Z);

    if (continue_step == 0) {
        init_arrays<<<grid, block>>>(d_f, d_rho, d_ux, d_uy, d_uz,
                                     d_force_realx, d_force_realy, d_force_realz);
        CHECK_CUDA_ERROR(cudaGetLastError());
        CHECK_CUDA_ERROR(cudaDeviceSynchronize());
    } else {
        char dirname[256];
        snprintf(dirname, sizeof(dirname), "data");

        char data_filename[512];
        snprintf(data_filename, sizeof(data_filename), "%s/%09d.dat", dirname, continue_step);

        FILE *f_data = fopen(data_filename, "rb");
        if (f_data == NULL) {
            fprintf(stderr, "Cannot open file: %s\n", data_filename);
            exit(EXIT_FAILURE);
        }

        size_t read_elements;
        read_elements = fread(h_rho, sizeof(double), LXYZ, f_data);
        if (read_elements != LXYZ) {
            fprintf(stderr, "Error reading h_rho data: read %zu of %d elements\n", 
                    read_elements, LXYZ);
            fclose(f_data);
            exit(EXIT_FAILURE);
        }

        read_elements = fread(h_ux, sizeof(double), LXYZ, f_data);
        if (read_elements != LXYZ) {
            fprintf(stderr, "Error reading h_ux data: read %zu of %d elements\n", 
                    read_elements, LXYZ);
            fclose(f_data);
            exit(EXIT_FAILURE);
        }

        read_elements = fread(h_uy, sizeof(double), LXYZ, f_data);
        if (read_elements != LXYZ) {
            fprintf(stderr, "Error reading h_uy data: read %zu of %d elements\n", 
                    read_elements, LXYZ);
            fclose(f_data);
            exit(EXIT_FAILURE);
        }

        read_elements = fread(h_uz, sizeof(double), LXYZ, f_data);
        if (read_elements != LXYZ) {
            fprintf(stderr, "Error reading h_uz data: read %zu of %d elements\n", 
                    read_elements, LXYZ);
            fclose(f_data);
            exit(EXIT_FAILURE);
        }

        fclose(f_data);

        CHECK_CUDA_ERROR(cudaMemcpy(d_rho, h_rho, LXYZ * sizeof(double), cudaMemcpyHostToDevice));
        CHECK_CUDA_ERROR(cudaMemcpy(d_ux,  h_ux,  LXYZ * sizeof(double), cudaMemcpyHostToDevice));
        CHECK_CUDA_ERROR(cudaMemcpy(d_uy,  h_uy,  LXYZ * sizeof(double), cudaMemcpyHostToDevice));
        CHECK_CUDA_ERROR(cudaMemcpy(d_uz,  h_uz,  LXYZ * sizeof(double), cudaMemcpyHostToDevice));

        char f_filename[512];
        snprintf(f_filename, sizeof(f_filename), "%s/f%09d.dat", dirname, continue_step);

        FILE *f_f = fopen(f_filename, "rb");
        if (f_f == NULL) {
            fprintf(stderr, "Cannot open file: %s\n", f_filename);
            exit(EXIT_FAILURE);
        }

        read_elements = fread(h_f, sizeof(double), NPOP * LXYZ, f_f);
        if (read_elements != NPOP * LXYZ) {
            fprintf(stderr, "Error reading h_f data: read %zu of %d elements\n", 
                    read_elements, NPOP * LXYZ);
            fclose(f_f);
            exit(EXIT_FAILURE);
        }

        fclose(f_f);
        CHECK_CUDA_ERROR(cudaMemcpy(d_f, h_f, NPOP * LXYZ * sizeof(double), cudaMemcpyHostToDevice));

        forcing<<<grid, block>>>(d_force_realx, d_force_realy, d_force_realz);
        CHECK_CUDA_ERROR(cudaGetLastError());
        CHECK_CUDA_ERROR(cudaDeviceSynchronize());
    }
}

__global__ void init_arrays(double *f, double *rho, double *ux, double *uy, double *uz, 
                            double *force_realx, double *force_realy, double *force_realz) {
    int ix = blockIdx.x * blockDim.x + threadIdx.x;
    int iy = blockIdx.y * blockDim.y + threadIdx.y;
    int iz = blockIdx.z * blockDim.z + threadIdx.z;
    if (ix >= LX || iy >= LY || iz >= LZ) return;

    int idx = (iz * LY + iy) * LX + ix;

    // Initialize macroscopic variables
    rho[idx] = 1.0;

    double A9 = 0.1*sqrt(2.0);
    double alpha = 1.0;
    double beta9 = 1.0;
    double cc = 60.0;
    double ccc1 = - LX / d_pi2 / alpha / d_ystar * A9 * d_ustar / cc / cc; 

    //mean flow
    uy[idx] = 0.0;
    uz[idx] = 0.0;
    double z9 = d_pi2 * (iz + 0.5) / LZ;
    double x9 = d_pi2 * (ix + 0.5) / LX;
    if (iy <= LY/2) {
        double yplus = (iy + 0.5) / d_ystar;
        if (yplus <= 10.8) {
            ux[idx] = d_ustar * yplus;
        } else {
            ux[idx] = log10(yplus) / 0.41 + 5.0;
            ux[idx] *= d_ustar;
        }
        double ccc9 = exp(-yplus/cc);
        double uu9 = ccc1 * yplus * sin(alpha * x9 + beta9 * z9);
        ux[idx] += uu9;
    } else if (iy > LY/2) {
        double yplus = (LY - iy - 0.5)/ d_ystar;
        if (yplus <= 10.8) {
            ux[idx] = d_ustar * yplus;
        } else {
            ux[idx] = log10(yplus) / 0.41 + 5.0;
            ux[idx] *= d_ustar;
        }
        double ccc9 = exp(-yplus/cc);
        double uu9 = ccc1 * yplus * sin(alpha * x9 + beta9 * z9);
        ux[idx] += uu9;
    }

    // Initialize force field
    force_realx[idx] = d_force_in_x;
    force_realy[idx] = 0.0;
    force_realz[idx] = 0.0;

    // Initialize distribution functions
    double usq = 0.0;  // u^2 + v^2 + w^2 = 0 initially
    for (int ip = 0; ip < NPOP; ip++) {
        double feq = d_tp[ip] * rho[idx] * (1.0 - 1.5 * usq);
        f[ip * LXYZ + idx] = feq - 0.5 * force_realx[idx] * d_cix[ip] / (1.0/3.0) * feq;
    }
}