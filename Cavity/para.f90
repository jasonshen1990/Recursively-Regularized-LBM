      subroutine para
      use var_inc
      implicit none

      integer i,j,k  

      pi  = 4.0*atan(1.0)   
      pi2 = 2.0*pi
      RT = 1./3.; cs = sqrt(RT)
      cs3 = RT*cs; cs4 = RT*RT
      cs5 = RT*RT*cs; cs6 = RT*RT*RT


      istep0 = 0
      istep00 = 1 
      imovie = 0

!      istpload = 1001       !!!!THIS NEEDS TO BE CHANGED WHEN STARTING NEW RUNS
      nsteps = 50

      newrun = .true.
!       newrun = .false.
!      newinitflow = .true.
!       newinitflow = .false.
      u0 = 0.1d0
      visc = u0*real(nx)/Re
      tau = 3.d0*visc + 0.5d0
      omega = 1.d0/(3.d0*visc + 0.5d0)
      force_in_y = 0.0 !8.d0*visc*u0/real(nx)/real(nx)

      ww0 = 8.d0/27.d0
      ww1 = 2.d0/27.d0
      ww2 = 1.d0/54.d0
      ww3 = 1.d0/216.d0

      !Dir:   0  1  2  3   4  5   6  7   8   9  10  11  12  13  14 15  16  17  18  19  20  21  22  23  24  25  26
      cix = (/0, 1,-1, 0,  0, 0,  0, 1,  1, -1, -1,  1,  1, -1, -1, 0,  0,  0,  0,  1,  1,  1,  1, -1, -1, -1, -1/)
      ciy = (/0, 0, 0, 1, -1, 0,  0, 1, -1,  1, -1,  0,  0,  0,  0, 1,  1, -1, -1,  1,  1, -1, -1,  1,  1, -1, -1/)
      ciz = (/0, 0, 0, 0,  0, 1, -1, 0,  0,  0,  0,  1, -1,  1, -1, 1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1/)
      tp = (/ww0,ww1,ww1,ww1,ww1,ww1,ww1,ww2,ww2,ww2,ww2,ww2,ww2,ww2,ww2,ww2,ww2,ww2,ww2,ww3,ww3,ww3,ww3,ww3,ww3,ww3,ww3/)
      !Arrray used for finding the opposite velocity
      ipopp=(/0, 2, 1, 4,  3, 6,  5,10,  9,  8,  7, 14, 13, 12, 11,18, 17, 16, 15, 26, 25, 24, 23, 22, 21, 20, 19/)
      !Arrray used for swap collision
      ! ipswap=(/2,4,6,8,10,12,14,16,18,20,22,24,26/)
      ! ipstay=(/0,1,3,5,7,9,11,13,15,17,19,21,23,25/) 

      ly = ny         !local division of dist for procs in y dir
      lz = nz         !local division of dist for procs in z dir


      end subroutine para
!==================================================================

      subroutine allocarray
      use var_inc
      implicit none

      allocate (f(0:npop-1,lx,ly,lz))
      allocate (rho(lx,ly,lz))
      allocate (ux(lx,ly,lz))
      allocate (uy(lx,ly,lz))
      allocate (uz(lx,ly,lz))

      allocate(force_realx(lx,ly,lz))
      allocate(force_realy(lx,ly,lz))
      allocate(force_realz(lx,ly,lz))

      end subroutine allocarray
!==================================================================



