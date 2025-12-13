      module var_inc
      implicit none

      integer,parameter:: nx = 128, ny = 128, nz = 2
      integer,parameter:: lx = nx
      integer,parameter:: lxh = lx/2, lyh = ny/2
      integer,parameter:: nxh = nx/2, nyh = ny/2, nzh = nz/2
      integer,parameter:: npop = 27
      integer,parameter:: iflowseed = 54321
      integer,parameter:: NTAB = 32 
      integer,parameter:: ndiag = 500, nstat = 500, nspec=100
      integer,parameter:: nflowout = 1000, npartout = 1000   
      integer,parameter:: nmovieout = 20000000, nsij = 1000     
      integer,parameter:: irelease = 200000
      integer,parameter:: iprocrate = 2  

      real,parameter:: rho0 = 1.0
      real,parameter:: Re = 1000.d0

      integer istep, istep0, istep00, nsteps, istpload, imovie   
      integer lz, ly

      logical newrun

      real alpha17,beta17,beta18,beta19,gamma20,gamma21,gamma22
      real alpha26,beta26,esi23,esi24,esi25
      real pi, pi2, anu, visc, tau, u0in, u0, omega
      real force_in_y
      real ww0, ww1, ww2, ww3
 
      integer,dimension(0:npop-1) :: cix, ciy, ciz, ipopp
      real, dimension(0:npop-1) :: tp 
      real, dimension(0:npop-1) :: H000, H100, H010, H001  ! zeroth and first order
      real, dimension(0:npop-1) :: H110, H101, H011, H200, H020, H002  ! Second
      real, dimension(0:npop-1) :: H210, H201, H120, H021, H102, H012, H111  ! Third
      real, dimension(0:npop-1) :: H220, H202, H022, H211, H121, H112 ! Fourth
      real, dimension(0:npop-1) :: H221, H212, H122, H222  ! Fifth and sixth
      real RT, cs, cs3, cs4, cs5, cs6

      real,allocatable,dimension(:,:,:):: force_realx,force_realy, &
                                            force_realz

      real,allocatable,dimension(:,:,:,:):: f
      real,allocatable,dimension(:,:,:):: rho
      real,allocatable,dimension(:,:,:):: ux, uy, uz


      character(len=80):: dirgenr, dirinitflow
      character(len=80):: dirdiag, dirstat   
      character(len=80):: dircntdflow, dircntdpart
      character(len=80):: dirflowout, dirpartout    
      character(len=80):: dirmoviedata


      end module var_inc
