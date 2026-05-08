# Recursively-Regularized-LBM
Implementations of recursively regularized LBM for simulating various flows, will be fully availabe after the acceptace of our manuscript.

Currently, methods using 3rd-order equilibrium and 3rd non-equilibrium are available.

- Updates will be made accordingly, including the CUDA-C++ version.

## Taylor-Green vortex flow
- Quasi two-dimensional simulation using 3D code, periodic in all the directions

- $2\times2$ vortices

<img src="Taylor_Green/Vorticity.png" alt="TGvortex" width="299"/>

## Lid-driven cavity flow

## Laminar channel flow

## Turbulent channel flow （DNS or LES）
BKG is the DNS code for turbulent channel flow using SRT.

In "RR27V33_code_DNS_LES/parameters.h", use 0 for DNS and 1 for LES-WALE.

If you have any issues, please feel free to contact me at shenj@sustech.edu.cn.
