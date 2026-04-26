#include "LBMLES.h"

void diag_flow(int istep) {
    char dirname[256];
    snprintf(dirname, sizeof(dirname), "data");


    struct stat st = {0};
    if (stat(dirname, &st) == -1) {
        #ifdef _WIN32
            _mkdir(dirname);
        #else
            mkdir(dirname, 0700);
        #endif
        printf("Created directory: %s\n", dirname);
    }

    // Copy data from device to host
    CHECK_CUDA_ERROR(cudaMemcpy(h_rho,  d_rho,  LXYZ * sizeof(double), cudaMemcpyDeviceToHost));
    CHECK_CUDA_ERROR(cudaMemcpy(h_ux,  d_ux,  LXYZ * sizeof(double), cudaMemcpyDeviceToHost));
    CHECK_CUDA_ERROR(cudaMemcpy(h_uy,  d_uy,  LXYZ * sizeof(double), cudaMemcpyDeviceToHost));
    CHECK_CUDA_ERROR(cudaMemcpy(h_uz,  d_uz,  LXYZ * sizeof(double), cudaMemcpyDeviceToHost));

    // Calculate statistics
    double umean = 0.0, vmean = 0.0, wmean = 0.0, rhomean = 0.0;
    for (int i = 0; i < LXYZ; i++) {
        umean += h_ux[i];
        vmean += h_uy[i];
        wmean += h_uz[i];
        rhomean += h_rho[i];
    }
    umean /= LXYZ;
    vmean /= LXYZ;
    wmean /= LXYZ;
    rhomean /= LXYZ;

    // Output to console
    printf("Step %-7d umean = %.17e, vmean = %.17e, wmean = %.17e, rhomean = %14.6e\n", 
            istep, umean/u0, vmean, wmean, rhomean);
    fflush(stdout);

    // Output to statistics file
    char stat_filename[512];
    snprintf(stat_filename, sizeof(stat_filename), "%s/statistics.dat", dirname);
    FILE *f_stat = fopen(stat_filename, "a");
    if (f_stat == NULL) {
        fprintf(stderr, "Error opening file: %s\n", stat_filename);
        exit(EXIT_FAILURE);
    }

    int write_status = fprintf(f_stat, "%-7d %14.6e %14.6e %14.6e %14.6e\n", 
                             istep, umean, vmean, wmean, rhomean);
    if (write_status < 0) {
        fprintf(stderr, "Error writing to statistics file\n");
        fclose(f_stat);
        exit(EXIT_FAILURE);
    }

    if (fflush(f_stat) != 0) {
        fprintf(stderr, "Error flushing statistics file\n");
        fclose(f_stat);
        exit(EXIT_FAILURE);
    }

    fclose(f_stat);
}

void output_flow(int istep) {
    char dirname[256];
    // snprintf(dirname, sizeof(dirname), "Ra%.1ePr%.2f", rayl, prand);
    snprintf(dirname, sizeof(dirname), "data");

    struct stat st = {0};
    if (stat(dirname, &st) == -1) {
        #ifdef _WIN32
            _mkdir(dirname);
        #else
            mkdir(dirname, 0700);
        #endif
        printf("Created directory: %s\n", dirname);
    }

    // Copy data from device to host
    CHECK_CUDA_ERROR(cudaMemcpy(h_rho, d_rho, LXYZ * sizeof(double), cudaMemcpyDeviceToHost));
    CHECK_CUDA_ERROR(cudaMemcpy(h_ux,  d_ux,  LXYZ * sizeof(double), cudaMemcpyDeviceToHost));
    CHECK_CUDA_ERROR(cudaMemcpy(h_uy,  d_uy,  LXYZ * sizeof(double), cudaMemcpyDeviceToHost));
    CHECK_CUDA_ERROR(cudaMemcpy(h_uz,  d_uz,  LXYZ * sizeof(double), cudaMemcpyDeviceToHost));

    // Output flow field data to binary file
    char data_filename[512];
    snprintf(data_filename, sizeof(data_filename), "%s/%09d.dat", dirname, istep);

    FILE *f_data = fopen(data_filename, "wb");
    if (f_data == NULL) {
        fprintf(stderr, "Cannot open file: %s\n", data_filename);
        exit(EXIT_FAILURE);
    }

    size_t written_elements;
    written_elements = fwrite(h_rho, sizeof(double), LXYZ, f_data);  //0
    if (written_elements != LXYZ) {
        fprintf(stderr, "Error writing h_rho data: wrote %zu of %d elements\n", 
                written_elements, LXYZ);
        fclose(f_data);
        exit(EXIT_FAILURE);
    }
    
    written_elements = fwrite(h_ux, sizeof(double), LXYZ, f_data);  //1
    if (written_elements != LXYZ) {
        fprintf(stderr, "Error writing h_ux data: wrote %zu of %d elements\n", 
                written_elements, LXYZ);
        fclose(f_data);
        exit(EXIT_FAILURE);
    }
    
    written_elements = fwrite(h_uy, sizeof(double), LXYZ, f_data);  //2
    if (written_elements != LXYZ) {
        fprintf(stderr, "Error writing h_uy data: wrote %zu of %d elements\n", 
                written_elements, LXYZ);
        fclose(f_data);
        exit(EXIT_FAILURE);
    }
    
    written_elements = fwrite(h_uz, sizeof(double), LXYZ, f_data);  //3
    if (written_elements != LXYZ) {
        fprintf(stderr, "Error writing h_uz data: wrote %zu of %d elements\n", 
                written_elements, LXYZ);
        fclose(f_data);
        exit(EXIT_FAILURE);
    }

    
    if (fflush(f_data) != 0) {
        fprintf(stderr, "Error flushing output file\n");
        fclose(f_data);
        exit(EXIT_FAILURE);
    }
    
    fclose(f_data);
}

void output_profile(int istep) {
    char dirname[256];
    snprintf(dirname, sizeof(dirname), "data");

    struct stat st = {0};
    if (stat(dirname, &st) == -1) {
        #ifdef _WIN32
            _mkdir(dirname);
        #else
            mkdir(dirname, 0700);
        #endif
        printf("Created directory: %s\n", dirname);
    }

    char filename[512]; 

    snprintf(filename, sizeof(filename), "%s/Profiles.txt", dirname);

    FILE *f_profile = fopen(filename, "a");
    if (f_profile == NULL) {
        fprintf(stderr, "Error opening file: %s\n", filename);
        exit(EXIT_FAILURE);
    }

    // Copy data from device to host
    CHECK_CUDA_ERROR(cudaMemcpy(h_rho, d_rho, LXYZ * sizeof(double), cudaMemcpyDeviceToHost));
    CHECK_CUDA_ERROR(cudaMemcpy(h_ux,  d_ux,  LXYZ * sizeof(double), cudaMemcpyDeviceToHost));
    CHECK_CUDA_ERROR(cudaMemcpy(h_uy,  d_uy,  LXYZ * sizeof(double), cudaMemcpyDeviceToHost));
    CHECK_CUDA_ERROR(cudaMemcpy(h_uz,  d_uz,  LXYZ * sizeof(double), cudaMemcpyDeviceToHost));
    CHECK_CUDA_ERROR(cudaMemcpy(h_taut,  d_taut,  LXYZ * sizeof(double), cudaMemcpyDeviceToHost));

    // Write the current step to the file
    int write_status = fprintf(f_profile, "%d\n", istep);
    if (write_status < 0) {
        fprintf(stderr, "Error writing step number to profile file\n");
        fclose(f_profile);
        exit(EXIT_FAILURE);
    }

    // Loop over y-planes
    for (int iy = 0; iy < LY; iy++) {
        double sum_ux = 0.0, sum_uy = 0.0, sum_uz = 0.0, sum_taut = 0.0;
        double s_uu=0, s_vv=0, s_ww=0, s_uv=0, s_uw=0, s_vw=0;

        int count = 0;

        for (int ix = 0; ix < LX; ix++) {
            for (int iz = 0; iz < LZ; iz++) {
                int idx = (iz * LY + iy) * LX + ix;
                sum_ux    += h_ux[idx];
                sum_uy    += h_uy[idx];
                sum_uz    += h_uz[idx];
                sum_taut    += h_taut[idx];
                count++;
            }
        }

        double avg_ux  = sum_ux  / count;
        double avg_uy  = sum_uy  / count;
        double avg_uz  = sum_uz  / count;
        double avg_taut  = sum_taut  / count;

        for (int ix = 0; ix < LX; ix++) {
            for (int iz = 0; iz < LZ; iz++) {
                int idx = (iz * LY + iy) * LX + ix;
                double ur = h_ux[idx] - avg_ux, vr = h_uy[idx] - avg_uy, wr = h_uz[idx] - avg_uz;
                s_uu += ur*ur; s_vv += vr*vr; s_ww += wr*wr;
                s_uv += ur*vr; s_uw += ur*wr; s_vw += vr*wr;
            }
        }

        fprintf(f_profile, "%d, %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e\n",
                iy + 0.5, (iy + 0.5)/ystar, avg_ux, avg_uy, avg_uz, s_uu/count, s_vv/count, s_ww/count, s_uv/count, s_uw/count, s_vw/count, avg_taut);

    }

    if (fflush(f_profile) != 0) {
        fprintf(stderr, "Error flushing profile file\n");
        fclose(f_profile);
        exit(EXIT_FAILURE);
    }

    fclose(f_profile);
}

void output_f(int istep) {
    char dirname[256];
    snprintf(dirname, sizeof(dirname), "data");

    char f_filename[512];
    snprintf(f_filename, sizeof(f_filename), "%s/f%09d.dat", dirname, istep);

    CHECK_CUDA_ERROR(cudaMemcpy(h_f, d_f, NPOP * LXYZ * sizeof(double), cudaMemcpyDeviceToHost));

    FILE *f_f = fopen(f_filename, "wb");
    if (f_f == NULL) {
        fprintf(stderr, "Cannot open file: %s\n", f_filename);
        exit(EXIT_FAILURE);
    }

    size_t written_elements;
    written_elements = fwrite(h_f, sizeof(double), NPOP * LXYZ, f_f);
    if (written_elements != NPOP * LXYZ) {
        fprintf(stderr, "Error writing h_f data: wrote %zu of %d elements\n", 
                written_elements, NPOP * LXYZ);
        fclose(f_f);
        exit(EXIT_FAILURE);
    }
    
    if (fflush(f_f) != 0) {
        fprintf(stderr, "Error flushing f file\n");
        fclose(f_f);
        exit(EXIT_FAILURE);
    }
    
    fclose(f_f);
}
