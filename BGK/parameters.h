#ifndef PARAMETERS_H
#define PARAMETERS_H

// Continue computation
#define CONTINUE_STEP 0

// Grid dimensions
#define LX 400
#define LY 199
#define LZ 200
#define LXYZ (LX * LY * LZ)
#define NPOP 27

// Time steps
#define NEND 10
#define NDIAG 10
#define NFLOWOUT 5000
#define NNUOUT 10

// Physical parameters
// #define RAYLEIGH 1.0e6
// #define PRANDTL 0.71
// #define T_HOT 1.0
// #define T_COLD 0.0
#define GRAVITY 1.0
#define BETA (0.1 / (LY * 5.0))

// Block dimensions for CUDA kernels
#define BLOCK_X 8
#define BLOCK_Y 4
#define BLOCK_Z 4

#endif 