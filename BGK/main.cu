#include "rb3d.h"

double random_0_to_1() {
    return (double)rand() / RAND_MAX;
}

int main() {
    printf( "Non-dimensional time, %.6e\n", LY / u0 /2);
    printf( "Driving force, %.6e\n", force_in_x);


    struct timeval start_time, end_time;
    double cpu_time_used = 0.0, gpu_time_used = 0.0;
    cudaEvent_t gpu_start, gpu_stop;
    CHECK_CUDA_ERROR(cudaEventCreate(&gpu_start));
    CHECK_CUDA_ERROR(cudaEventCreate(&gpu_stop));
    gettimeofday(&start_time, NULL);

    initialize();

    CHECK_CUDA_ERROR(cudaMemcpy(h_rho, d_rho, LXYZ * sizeof(double), cudaMemcpyDeviceToHost));
    CHECK_CUDA_ERROR(cudaMemcpy(h_ux,  d_ux,  LXYZ * sizeof(double), cudaMemcpyDeviceToHost));
    CHECK_CUDA_ERROR(cudaMemcpy(h_uy,  d_uy,  LXYZ * sizeof(double), cudaMemcpyDeviceToHost));
    CHECK_CUDA_ERROR(cudaMemcpy(h_uz,  d_uz,  LXYZ * sizeof(double), cudaMemcpyDeviceToHost));

    dim3 grid((LX + BLOCK_X - 1) / BLOCK_X, (LY + BLOCK_Y - 1) / BLOCK_Y, (LZ + BLOCK_Z - 1) / BLOCK_Z);
    dim3 block(BLOCK_X, BLOCK_Y, BLOCK_Z);

    int npforcing = 250000;
    srand(12345);
    for (int istep = continue_step + 1; istep <= NEND; istep++) {
        // GPU timing start
        CHECK_CUDA_ERROR(cudaEventRecord(gpu_start));

        if (istep < npforcing){
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_istep,   &istep, sizeof(int)));
            
            if (istep % 2000 ==0) {
                double rand_phase = random_0_to_1();
                CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_rand_phase,   &rand_phase, sizeof(double)));
                printf( "Check_rand_phase_host,%.3e\n", rand_phase);
            }
            forcing_GP<<<grid, block>>>(d_force_realx, d_force_realy, d_force_realz);
        } else {
            forcing<<<grid, block>>>(d_force_realx, d_force_realy, d_force_realz);
        }

        collision_BGK<<<grid, block>>>(d_f, d_force_realx, d_force_realy, d_force_realz, 
                                     d_rho, d_ux, d_uy, d_uz);   // Collision

        CHECK_CUDA_ERROR(cudaGetLastError());
        CHECK_CUDA_ERROR(cudaDeviceSynchronize());

        streaming<<<grid, block>>>(d_f, d_f_temp);   // Streaming
        CHECK_CUDA_ERROR(cudaDeviceSynchronize());
        streaming2<<<grid, block>>>(d_f, d_f_temp); 
        CHECK_CUDA_ERROR(cudaMemcpy(d_f, d_f_temp, NPOP * LXYZ * sizeof(double), cudaMemcpyDeviceToDevice));

        CHECK_CUDA_ERROR(cudaGetLastError());
        CHECK_CUDA_ERROR(cudaDeviceSynchronize());

        macrovar<<<grid, block>>>(d_f, d_force_realx, d_force_realy, d_force_realz, 
                                 d_rho, d_ux, d_uy, d_uz);
        CHECK_CUDA_ERROR(cudaGetLastError());
        CHECK_CUDA_ERROR(cudaDeviceSynchronize());

        // GPU timing end
        CHECK_CUDA_ERROR(cudaEventRecord(gpu_stop));
        CHECK_CUDA_ERROR(cudaEventSynchronize(gpu_stop));

        float gpu_milliseconds = 0.0f;
        CHECK_CUDA_ERROR(cudaEventElapsedTime(&gpu_milliseconds, gpu_start, gpu_stop));
        gpu_time_used += (double)gpu_milliseconds / 1000.0;

        // CPU timing start
        clock_t cpu_start = clock();

        char dirname[256];
        snprintf(dirname, sizeof(dirname), "data");

        if (istep % NDIAG == 0) {
            diag_flow(istep);
        }
        if (istep % NFLOWOUT == 0 || istep == NEND) {
            output_flow(istep);
        }
        if (istep % NNUOUT == 0) {
            // output_nu(istep);
            output_profile(istep);
        }
        if (istep == NEND) {
            output_fg(istep);
        }

        // CPU timing end
        cpu_time_used += (double)(clock() - cpu_start) / CLOCKS_PER_SEC;
    }

    // calculate total time
    gettimeofday(&end_time, NULL);
    double total_time = (end_time.tv_sec - start_time.tv_sec) + 
                       (end_time.tv_usec - start_time.tv_usec) / 1000000.0;

    printf("Simulation completed.\n");
    printf("Total iterations: %d\n", NEND);
    printf("CPU computation time: %.2f seconds\n", cpu_time_used);
    printf("GPU computation time: %.2f seconds\n", gpu_time_used);
    printf("Total computation time: %.2f seconds\n", total_time);

    CHECK_CUDA_ERROR(cudaEventDestroy(gpu_start));
    CHECK_CUDA_ERROR(cudaEventDestroy(gpu_stop));

    CHECK_CUDA_ERROR(cudaFree(d_f));
    CHECK_CUDA_ERROR(cudaFree(d_f_temp));
    CHECK_CUDA_ERROR(cudaFree(d_rho));
    CHECK_CUDA_ERROR(cudaFree(d_ux));
    CHECK_CUDA_ERROR(cudaFree(d_uy));
    CHECK_CUDA_ERROR(cudaFree(d_uz));
    CHECK_CUDA_ERROR(cudaFree(d_force_realx));
    CHECK_CUDA_ERROR(cudaFree(d_force_realy));
    CHECK_CUDA_ERROR(cudaFree(d_force_realz));

    free(h_f);
    free(h_rho);
    free(h_ux);
    free(h_uy);
    free(h_uz);

    CHECK_CUDA_ERROR(cudaDeviceReset());

    return 0;
}