!=============================================================================
!@subroutine para
!@desc Initializes parameters, variables, and constants
!=============================================================================      
      subroutine para
      use var_inc
      ! use mpi
      implicit none

!      logical dirExist
      integer i,j,k
      real vmax

      nsteps = 60000

      s = 1 !Global switch used for 2-array algorithm
!=======================================================
! Physical constants and parameters
!=======================================================
      pi = 4.0*atan(1.0) 
      pi2 = 2.0*pi
      RT = 1/3.
      cs = sqrt(RT)
      cs3 = RT*cs; cs4 = RT*RT
      cs5 = RT*RT*cs; cs6 = RT*RT*RT
      
!====================================      
      ! Turbulent channel flow
!==================================== 
      ! specify the force magnitude
!     visc = 0.0032
!     Rstar = 200.0
      
      Rstar = 50.0
      ustar = 0.05
      visc = ustar*real(nx)/Rstar
      tau = 3.0*visc + 0.5
      force_in_y = 8.d0*visc*ustar/real(nx)/real(nx)

      ! force_in_y = 2.*rho0*ustar*ustar/real(nx)
      ! print*, force_in_y

      ystar = visc/ustar
      force_mag = 1.0


 
      
      ! print*, tau
      ! read(*,*)
!=======================================================
! Initialize MRT related constants
!=======================================================
     

      ww0 = 8.d0/27.d0
      ww1 = 2.d0/27.d0
      ww2 = 1.d0/54.d0
      ww3 = 1.d0/216.d0


      !Dir:   0 1  2 3  4 5  6 7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
      cix = (/0,1,-1,0, 0,0, 0,1, 1,-1,-1, 1, 1,-1,-1, 0, 0, 0, 0, 1, 1, 1, 1,-1,-1,-1,-1/)
      ciy = (/0,0, 0,1,-1,0, 0,1,-1, 1,-1, 0, 0, 0, 0, 1, 1,-1,-1, 1, 1,-1,-1, 1, 1,-1,-1/)
      ciz = (/0,0, 0,0, 0,1,-1,0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1/)
      tp = (/ww0,ww1,ww1,ww1,ww1,ww1,ww1,ww2,ww2,ww2,ww2,ww2,ww2,ww2,ww2,ww2,ww2,ww2,ww2,ww3,ww3,ww3,ww3,ww3,ww3,ww3,ww3/)

      !Arrray used for finding the opposite velocity
      ipopp=(/0,2,1,4,3,6,5,10,9,8,7,14,13,12,11,18,17,16,15,26,25,24,23,22,21,20,19/)
      !       0,1,2,3,4,5,6,7,8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26
      ipref=(/0,2,1,4,3,6,5,9,10,7, 8,13,14,11,12,15,16,17,18,23,24,25,26,19,20,21,22/)
      !Arrray used for swap collision
      ipswap=(/2,4,6,8,10,12,14,16,18,20,22,24,26/)
      ipstay=(/0,1,3,5,7,9,11,13,15,17,19,21,23,25/) 

!=======================================================
! Create MPI topology
!=======================================================
      ly = ny; lz = nz

      
      end subroutine para
!=============================================================================
!@subroutine allocarray
!@desc Allocates all globally used arrays in the var_inc module
!=============================================================================  
      subroutine allocarray
      use var_inc
      implicit none

      allocate (f(0:npop-1,lx,ly,lz))
      ! allocate (f_post_coll(0:npop-1,lx,ly,lz,2))
      allocate (rho(lx,ly,lz))
      allocate (rhop(lx,ly,lz))
      allocate (ux(lx,ly,lz))
      allocate (uy(lx,ly,lz))
      allocate (uz(lx,ly,lz))
      allocate (ox(lx,ly,lz))
      allocate (oy(lx,ly,lz))
      allocate (oz(lx,ly,lz))
      allocate (kx(lx+2,ly+lyext,lz))
      allocate (ky(lx+2,ly+lyext,lz))
      allocate (kz(lx+2,ly+lyext,lz))
      allocate (k2(lx+2,ly+lyext,lz))
      allocate (ik2(lx+2,ly+lyext,lz))
      allocate (vx(lx+2,ly+lyext,lz))
      allocate (vy(lx+2,ly+lyext,lz))
      allocate (vz(lx+2,ly+lyext,lz))
      allocate (wx(lx+2,ly+lyext,lz))
      allocate (wy(lx+2,ly+lyext,lz))
      allocate (wz(lx+2,ly+lyext,lz))
      allocate (ibnodes(0:lx+1,0:ly+1,0:lz+1))
      allocate(force_realx(lx,ly,lz))
      allocate(force_realy(lx,ly,lz))
      allocate(force_realz(lx,ly,lz))
      allocate(hydrofx(lx,ly,lz))
      allocate(hydrofy(lx,ly,lz))
      allocate(hydrofz(lx,ly,lz))
      allocate(hydrofx2(lx,ly,lz))
      allocate(hydrofy2(lx,ly,lz))
      allocate(hydrofz2(lx,ly,lz))
      allocate(pworkx(lx,ly,lz))
      allocate(pworky(lx,ly,lz))
      allocate(pworkz(lx,ly,lz))
      ibnodes = -1


      end subroutine allocarray
!=============================================================================
!@subroutine constructMPItypes
!@desc Creates any custom MPI datatypes used
!=============================================================================  
      ! subroutine constructMPItypes
      ! use mpi
      ! use var_inc
      ! implicit none

      ! call MPI_TYPE_CONTIGUOUS(4, MPI_INTEGER, MPI_IPF_NODE, ierr)
      ! call MPI_TYPE_COMMIT(MPI_IPF_NODE, ierr)
      
      ! call MPI_TYPE_CONTIGUOUS(3, MPI_INTEGER, MPI_FILL_NODE, ierr)
      ! call MPI_TYPE_COMMIT(MPI_FILL_NODE, ierr)

      ! end subroutine constructMPItypes
!=============================================================================
!@subroutine makedir
!@desc Makes sure that the given directory exists, if it does not it attempts
!      to create it.
!@param dirPath = character array which holds the directory path that needs
!                 to be checked
!=============================================================================  
      subroutine makedir (dirPath)
      use var_inc

      implicit none
      logical dirExist
      character(len=120):: dirPath !Needs to be the same as declared length!

      inquire( file=trim(DirPath)//'/.', exist=dirExist )  ! Works with gfortran
      ! inquire(directory = trim(dirPath), exist = dirExist)  ! Works with ifort (yellowstone)
      if(.NOT.dirExist)then
        write(*,*) trim(dirPath)//' not found. Creating...'
        call system('mkdir -p '// trim(dirPath)) !Execute system command to create directory
      endif

      end subroutine makedir
