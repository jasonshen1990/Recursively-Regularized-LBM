#ifndef PARAMETERS_H
#define PARAMETERS_H

// Continue computation
#define CONTINUE_STEP 0
#define USE_LES 0   // 1: use LES sub-grid shear stress model    or    0:  DNS

// Grid dimensions
#define LX 64
#define LY 64
#define LZ 64
#define LXYZ (LX * LY * LZ)
#define NPOP 27

// Time steps
#define NEND 1000000
#define NDIAG 1000
#define NFLOWOUT 100000
#define NOUT 1000

// Block dimensions for CUDA kernels
#define BLOCK_X 8
#define BLOCK_Y 8
#define BLOCK_Z 4

#endif 