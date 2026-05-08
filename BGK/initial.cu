#include "rb3d.h"

void initialize() {
    // Allocate device memory
    CHECK_CUDA_ERROR(cudaMalloc(&d_f,      NPOP * LXYZ * sizeof(double)));
    CHECK_CUDA_ERROR(cudaMalloc(&d_f_temp, NPOP * LXYZ * sizeof(double)));

    CHECK_CUDA_ERROR(cudaMalloc(&d_rho,           LXYZ * sizeof(double)));
    CHECK_CUDA_ERROR(cudaMalloc(&d_ux,            LXYZ * sizeof(double)));
    CHECK_CUDA_ERROR(cudaMalloc(&d_uy,            LXYZ * sizeof(double)));
    CHECK_CUDA_ERROR(cudaMalloc(&d_uz,            LXYZ * sizeof(double)));

    CHECK_CUDA_ERROR(cudaMalloc(&d_force_realx,   LXYZ * sizeof(double)));
    CHECK_CUDA_ERROR(cudaMalloc(&d_force_realy,   LXYZ * sizeof(double)));
    CHECK_CUDA_ERROR(cudaMalloc(&d_force_realz,   LXYZ * sizeof(double)));

    // Allocate host memory with error checking
    h_rho = (double*)malloc(LXYZ * sizeof(double));
    h_ux = (double*)malloc(LXYZ * sizeof(double));
    h_uy = (double*)malloc(LXYZ * sizeof(double));
    h_uz = (double*)malloc(LXYZ * sizeof(double));

    h_f = (double*)malloc(NPOP * LXYZ * sizeof(double));

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

    // Copy parameters to constant memory
    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_cix,   h_cix, NPOP * sizeof(int)));
    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_ciy,   h_ciy, NPOP * sizeof(int)));
    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_ciz,   h_ciz, NPOP * sizeof(int)));
    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_tp,    h_tp,  NPOP * sizeof(double)));

    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_tau,    &tau,         sizeof(double)));
    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_ustar,  &u0,    sizeof(double)));
    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_ystar,  &ystar, sizeof(double)));
    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_force_in_x,  &force_in_x,   sizeof(double)));
    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_pi,    &pi,          sizeof(double)));
    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_pi2,   &pi2,         sizeof(double)));

    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_grav0, &grav0,       sizeof(double)));
    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_beta,  &beta,        sizeof(double)));
    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_t0,    &t0,          sizeof(double)));
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

        size_t read_elements;
        read_elements = fread(h_rho, sizeof(double), LXYZ, f_data);
        if (read_elements != LXYZ) {
            fprintf(stderr, "Error reading h_rho data: read %zu of %d elements\n", 
                    read_elements, LXYZ);
            fclose(f_data);
            exit(EXIT_FAILURE);
        }

        read_elements = fread(h_ux, sizeof(double), LXYZ, f_data);
        read_elements = fread(h_uy, sizeof(double), LXYZ, f_data);
        read_elements = fread(h_uz, sizeof(double), LXYZ, f_data);
        fclose(f_data);

        CHECK_CUDA_ERROR(cudaMemcpy(d_rho, h_rho, LXYZ * sizeof(double), cudaMemcpyHostToDevice));
        CHECK_CUDA_ERROR(cudaMemcpy(d_ux,  h_ux,  LXYZ * sizeof(double), cudaMemcpyHostToDevice));
        CHECK_CUDA_ERROR(cudaMemcpy(d_uy,  h_uy,  LXYZ * sizeof(double), cudaMemcpyHostToDevice));
        CHECK_CUDA_ERROR(cudaMemcpy(d_uz,  h_uz,  LXYZ * sizeof(double), cudaMemcpyHostToDevice));

        char fg_filename[512];
        snprintf(fg_filename, sizeof(fg_filename), "%s/fg%09d.dat", dirname, continue_step);

        FILE *f_fg = fopen(fg_filename, "rb");
        if (f_fg == NULL) {
            fprintf(stderr, "Cannot open file: %s\n", fg_filename);
            exit(EXIT_FAILURE);
        }

        read_elements = fread(h_f, sizeof(double), NPOP * LXYZ, f_fg);
        if (read_elements != NPOP * LXYZ) {
            fprintf(stderr, "Error reading h_f data: read %zu of %d elements\n", 
                    read_elements, NPOP * LXYZ);
            fclose(f_fg);
            exit(EXIT_FAILURE);
        }

        fclose(f_fg);

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
    rho[idx] = 0.0;
    ux[idx]  = 0.0;
    uy[idx]  = 0.0;
    uz[idx]  = 0.0;

    // double A9 = 0.0;
    double A9 = 0.1*sqrt(2.0);
    double alpha = 1.0;
    double beta9 = 1.0;
    // double gamma1 = 2.0;
    // double gamma2 = 1.0;
    // double phi1 = 0.111;
    // double phi2 = 0.555;
    double cc = 60.0;
    // double ystar = visc / u0;

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
            // if(iy == LY/2-1 && ix == 5 && iz == 5){
            //     printf( "Check_ini_ux,%.3e\n", ux[idx]);
            // }
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
    force_realy[idx] = 0.0; //d_grav0 * d_beta * (phi[idx] - d_t0);
    force_realz[idx] = 0.0;

    // Initialize distribution functions
    double usq = 0.0;  // u^2 + v^2 + w^2 = 0 initially
    for (int ip = 0; ip < NPOP; ip++) {
        double feq = d_tp[ip] * rho[idx] * (1.0 - 1.5 * usq);
        f[ip * LXYZ + idx] = feq - 0.5 * force_realx[idx] * d_cix[ip] / (1.0/3.0) * feq;

    }
}