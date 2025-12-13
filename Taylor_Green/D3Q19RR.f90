program Taylor
    implicit none
    ! integer, parameter :: num_cases = 6
    ! integer :: L_values(num_cases) = [50, 80, 100, 150, 200, 300]
    ! integer :: iit_values(num_cases) = [1000, 1600, 2000, 3000, 4000, 6000]
    integer, parameter :: num_cases = 4
    integer :: L_values(num_cases) = [50, 70,80,100]
    integer :: iit_values(num_cases) = [10000, 14000,16000,20000]
    integer :: case_idx
    
    do case_idx = 1, num_cases
        call run_simulation(L_values(case_idx), iit_values(case_idx))
    end do

contains
subroutine run_simulation(L, iit)
    implicit none
    integer, intent(in) :: L, iit
    character(len=100) :: filename_v, filename_va, filename_main
    integer, parameter :: LZ=2
    integer, parameter :: ippQ = 19
    integer :: i, j, k, ip, jj1, jj2
    integer :: it, ipre, jpre, kpre
    real(kind=8) :: start_time, end_time, time
    real(kind=8) :: u9, v9, w9, rho9, f9, wwp9, ix9, iy9, iz9
    real(kind=8), parameter :: pi = 3.14159265358979323846d0
    real(kind=8) :: cs = 1.0d0 / SQRT(3.0d0)
    real(kind=8) :: RT = 1.0d0 / 3.0d0
    real(kind=8) :: cs2, cs3, cs4, cs5, cs6, tau, k_val=2.0d0
    real(kind=8) :: Re0 , visc , U0 , rho0 , dt, D, muV 
    real(kind=8), allocatable :: f(:,:,:,:), fpost(:,:,:,:), feq(:,:,:,:)
    real(kind=8), allocatable :: feq9(:), fneq9(:)

    real(kind=8), allocatable :: rho(:,:,:), u(:,:,:), v(:,:,:), w(:,:,:)
    real(kind=8), allocatable :: p(:,:,:), utheo(:,:,:), vtheo(:,:,:),wtheo(:,:,:), ptheo(:,:,:)
    real(kind=8), allocatable :: dudxi(:,:,:), dudyi(:,:,:), dvdxi(:,:,:), dvdyi(:,:,:)
    real(kind=8), allocatable :: rhoui(:,:,:), rhovi(:,:,:), drhoudxi(:,:,:)
    real(kind=8), allocatable :: drhoudyi(:,:,:), drhovdxi(:,:,:), drhovdyi(:,:,:)
    real(kind=8), allocatable :: dudx(:,:,:), dudy(:,:,:), dvdx(:,:,:), dvdy(:,:,:)
    real(kind=8), allocatable :: d2udxy(:,:,:), d2vdxy(:,:,:), d2vdy2(:,:,:)
    real(kind=8), allocatable :: rhou(:,:,:), rhov(:,:,:),rhow(:,:,:)
    real(kind=8), allocatable :: drhoudx(:,:,:), drhoudy(:,:,:), drhovdx(:,:,:), drhovdy(:,:,:)
    real(kind=8), allocatable :: drhodx(:,:,:), drhody(:,:,:), udrhodx(:,:,:), vdrhody(:,:,:)
    real(kind=8), allocatable :: d2udxytheo(:,:,:), d2vdxytheo(:,:,:)
    real uu, vv, ww
    real uv, eu, sec, thir, four, fifsix
    !equilibrium recursive relation coefficients
    real a110, a101, a011, a200, a020, a002  ! Second
    real a210, a201, a120, a021, a102, a012, a111  ! Third
    real a210012p, a210012m, a102120p, a102120m, a021201p, a021201m

    real a220, a202, a022, a211, a121, a112 ! Fourth
    real a221, a212, a122, a222  ! Fifth and sixth
    !non-equilibrium recursive relation coefficients
    real ane110, ane101, ane011, ane200, ane020, ane002  ! Second
    real ane210, ane201, ane120, ane021, ane102, ane012, ane111  ! Third
    real ane210012p, ane210012m, ane102120p, ane102120m, ane021201p, ane021201m

    real ane220, ane202, ane022, ane211, ane121, ane112 ! Fourth
    real ane221, ane212, ane122, ane222  ! Fifth and sixth
    real, dimension(ippQ) :: H000, H100, H010, H001  ! zeroth and first order
    real, dimension(ippQ) :: H110, H101, H011, H200, H020, H002  ! Second
    real, dimension(ippQ) :: H210, H201, H120, H021, H102, H012, H111  ! Third
    real, dimension(ippQ) :: H220, H202, H022, H211, H121, H112 ! Fourth
    real, dimension(ippQ) :: H221, H212, H122, H222  ! Fifth and sixth
    real, dimension(ippQ) :: H210012p, H210012m, H102120p, H102120m, H021201p, H021201m

    real(kind=8) :: termi1, termi2, termi3, termi31, termi32, termi33, termi34
    real(kind=8) :: termi4, termi41, termi42, termi43, termi44, termi45, termi46, fneq
    real(kind=8) :: term1, term2, term3, term31, term32, term33
    real(kind=8) :: term4, term41, term42, term43, term5, S_alpha
    real(kind=8) :: xx, yy, zz, x, y, z
    real(kind=8), dimension(L) :: vnum, vtheoline, pnum, udrhodxnum, vdrhodynum, d2vdy2num, d2vdxynum, dvdxnum, dvdynum 
    integer :: ip1, ip2, im1, im2, jp1, jp2, jm1, jm2, kp1, km1
    integer :: ix(ippQ) = [0,1,-1,0, 0,0, 0,1, 1,-1,-1, 1, 1,-1,-1, 0, 0, 0, 0]
    integer :: iy(ippQ) = [0,0, 0,1,-1,0, 0,1,-1, 1,-1, 0, 0, 0, 0, 1, 1,-1,-1]
    integer :: iz(ippQ) = [0,0, 0,0, 0,1,-1,0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1]

    real(kind=8) :: wwp(ippQ) = [1.0d0/3.0d0, 1.0d0/18.0d0, 1.0d0/18.0d0, 1.0d0/18.0d0, 1.0d0/18.0d0, 1.0d0/18.0d0, 1.0d0/18.0d0, &
                              1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0, &
                              1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0]
    
    Re0 = 200.0d0/cs
    U0 = 1.0d0*cs/ 100.0d0
    visc = U0 * L / Re0
    rho0 = 1.0d0
    cs2 = RT
    cs4 = RT**2
    cs3 = RT*cs
    cs5 = RT*RT*cs
    cs6 = RT*RT*RT
    tau = visc/cs2 + 0.5d0
    dt = 1.0d0
    D = 3.0d0 
    muV = 0.05d0

    write(filename_v, '(A,I0,A)') 'D3Q19RR33BMa0.01L', L_values(case_idx), 'v.txt'
    write(filename_va, '(A,I0,A)') 'D3Q19RR33BMa0.01L', L_values(case_idx), 'vA.txt'
    write(filename_main, '(A,I0,A)') 'D3Q19RR33BMa0.01L', L_values(case_idx), '.txt'

    ! 分配数组空间
    allocate(f(ippQ, L, L, LZ))
    allocate(feq9(ippQ), fneq9(ippQ))
    allocate(fpost(ippQ, L, L, LZ))
    allocate(feq(ippQ, L, L, LZ))
    allocate(rho(L, L, LZ), u(L, L, LZ), v(L, L, LZ),w(L, L, LZ), p(L, L, LZ))
    allocate(utheo(L, L, LZ), vtheo(L, L, LZ), ptheo(L, L, LZ))
    allocate(dudxi(L, L, LZ), dudyi(L, L, LZ), dvdxi(L, L, LZ), dvdyi(L, L, LZ))
    allocate(rhoui(L, L, LZ), rhovi(L, L, LZ), drhoudxi(L, L, LZ), drhoudyi(L, L, LZ), drhovdxi(L, L, LZ), drhovdyi(L, L, LZ))
    allocate(dudx(L, L, LZ), dudy(L, L, LZ), dvdx(L, L, LZ), dvdy(L, L, LZ),d2udxy(L, L, LZ),d2vdxy(L, L, LZ),d2vdy2(L, L, LZ))
    allocate(rhou(L, L, LZ), rhov(L, L, LZ),rhow(L, L, LZ), drhoudx(L, L, LZ))
    allocate(drhoudy(L, L, LZ), drhovdx(L, L, LZ), drhovdy(L, L, LZ))
    allocate(drhodx(L, L, LZ), drhody(L, L, LZ), udrhodx(L, L, LZ), vdrhody(L, L, LZ))
    allocate(d2udxytheo(L, L, LZ),d2vdxytheo(L, L, LZ))
    ! Initialization
    call cpu_time(start_time)

    !Hermite_polynomial
    H000 = 1.
    do ip = 1, ippQ 
        H100(ip) = ix(ip) / cs
        H010(ip) = iy(ip) / cs
        H001(ip) = iz(ip) / cs
        H110(ip) = ix(ip) * iy(ip) / RT
        H101(ip) = ix(ip) * iz(ip) / RT
        H011(ip) = iy(ip) * iz(ip) / RT
        H200(ip) = ix(ip) * ix(ip) / RT - 1
        H020(ip) = iy(ip) * iy(ip) / RT - 1
        H002(ip) = iz(ip) * iz(ip) / RT - 1
        !third
        H210012p(ip)=(ix(ip) * ix(ip) * iy(ip) / cs3 - iy(ip) / cs) + (iz(ip) * iz(ip) * iy(ip) / cs3 - iy(ip) / cs)
        H210012m(ip)=(ix(ip) * ix(ip) * iy(ip) / cs3 - iy(ip) / cs) - (iz(ip) * iz(ip) * iy(ip) / cs3 - iy(ip) / cs)
        H102120p(ip)=(iz(ip) * iz(ip) * ix(ip) / cs3 - ix(ip) / cs) + (iy(ip) * iy(ip) * ix(ip) / cs3 - ix(ip) / cs)
        H102120m(ip)=(iz(ip) * iz(ip) * ix(ip) / cs3 - ix(ip) / cs) - (iy(ip) * iy(ip) * ix(ip) / cs3 - ix(ip) / cs)
        H021201p(ip)=(iy(ip) * iy(ip) * iz(ip) / cs3 - iz(ip) / cs) + (ix(ip) * ix(ip) * iz(ip) / cs3 - iz(ip) / cs)
        H021201m(ip)=(iy(ip) * iy(ip) * iz(ip) / cs3 - iz(ip) / cs) - (ix(ip) * ix(ip) * iz(ip) / cs3 - iz(ip) / cs)
    end do

    do k = 1, LZ
       do j = 1, L
           do i = 1, L
               u(i,j,k) = -U0*cos(k_val*pi*real(i-0.5, kind=8)/L)*sin(k_val*pi*real(j-0.5, kind=8)/L)
               v(i,j,k) = U0*cos(k_val*pi*real(j-0.5, kind=8)/L)*sin(k_val*pi*real(i-0.5, kind=8)/L)
               w(i,j,k) = 0.d0
               rho(i,j,k) = 1.0d0 - U0*U0/2*cos(k_val*pi/L*(real(i-0.5, kind=8) -&
                           real(j-0.5, kind=8)))*cos(k_val*pi/L*(real(i-0.5, kind=8) +&
                           real(j-0.5, kind=8)))/cs2
    
           end do
       end do
    end do

    do k = 1, LZ
        do j = 1, L
            do i = 1, L
                x = real(i-0.5, kind=8)
                y = real(j-0.5, kind=8)
    
                dudxi(i,j,k) = (k_val*pi/L)*U0*sin(k_val*pi*x/L)*sin(k_val*pi*y/L)
                dudyi(i,j,k) = - (k_val*pi/L)*U0*cos(k_val*pi*x/L)*cos(k_val*pi*y/L)
                dvdxi(i,j,k) = (k_val*pi/L)*U0*cos(k_val*pi*y/L)*cos(k_val*pi*x/L)
                dvdyi(i,j,k) = - (k_val*pi/L)*U0*sin(k_val*pi*y/L)*sin(k_val*pi*x/L)
    
                drhoudxi(i,j,k) = (-U0**2 * k_val * pi / (cs**2 * L) * sin(2*k_val*pi*x/L)) * &
                        (-U0 * cos(k_val*pi*x/L) * sin(k_val*pi*y/L)) + &
                        (1.0d0 - U0**2 / (2 * cs**2) * cos(k_val*pi*(x-y)/L) * cos(k_val*pi*(x+y)/L)) *&
                         (U0 * k_val * pi / L * sin(k_val*pi*x/L) * sin(k_val*pi*y/L))
                drhoudyi(i,j,k) = (-U0**2 * k_val * pi / (cs**2 * L) * (-sin(2*k_val*pi*y/L))) * &
                        (-U0 * cos(k_val*pi*x/L) * sin(k_val*pi*y/L)) + &
                        (1.0d0 - U0**2 / (2 * cs**2) * cos(k_val*pi*(x-y)/L) * cos(k_val*pi*(x+y)/L)) *&
                         (-U0 * cos(k_val*pi*y/L) * k_val * pi / L * cos(k_val*pi*x/L))
                drhovdxi(i,j,k) = (-U0**2 * k_val * pi / (cs**2 * L) * sin(2*k_val*pi*x/L)) * &
                        (U0 * sin(k_val*pi*x/L) * cos(k_val*pi*y/L)) + &
                        (1.0d0 - U0**2 / (2 * cs**2) * cos(k_val*pi*(x-y)/L) * cos(k_val*pi*(x+y)/L)) *&
                         (U0 * cos(k_val*pi*x/L) * k_val * pi / L * cos(k_val*pi*y/L))
                drhovdyi(i,j,k)  = (-U0**2 * k_val * pi / (cs**2 * L) * (-sin(2*k_val*pi*y/L))) *&
                        (U0 * sin(k_val*pi*x/L) * cos(k_val*pi*y/L)) + &
                        (1.0d0 - U0**2 / (2 * cs**2) * cos(k_val*pi*(x-y)/L) * cos(k_val*pi*(x+y)/L)) *&
                         (-U0 * sin(k_val*pi*y/L) * k_val * pi / L * sin(k_val*pi*x/L))
            end do
        end do
    end do

    do k = 1, LZ
        do j = 1, L
            do i = 1, L           
                do ip = 1, ippQ
                        u9 = u(i,j,k)
                        v9 = v(i,j,k)
                        w9 = w(i,j,k)
                        rho9 = rho(i,j,k)

                        ip1 = mod(i, L) + 1
                        ip2 = mod(i+1, L) + 1
                        im1 = mod(i-2+L, L) + 1 
                        im2 = mod(i-3+L, L) + 1
    
                        jp1 = mod(j, L) + 1
                        jp2 = mod(j+1, L) + 1
                        jm1 = mod(j-2+L, L) + 1
                        jm2 = mod(j-3+L, L) + 1
    
                        dudx(i,j,k) = - (u(ip2,j,k) - 8.0*u(ip1,j,k) + 8.0*u(im1,j,k) - u(im2,j,k)) / 12.0
                        dudy(i,j,k) = - (u(i,jp2,k) - 8.0*u(i,jp1,k) + 8.0*u(i,jm1,k) - u(i,jm2,k)) / 12.0
                        dvdx(i,j,k) = - (v(ip2,j,k) - 8.0*v(ip1,j,k) + 8.0*v(im1,j,k) - v(im2,j,k)) / 12.0
                        dvdy(i,j,k) = - (v(i,jp2,k) - 8.0*v(i,jp1,k) + 8.0*v(i,jm1,k) - v(i,jm2,k)) / 12.0

                    !equilibrium distribution function & a^eq
                        uu = u9 / cs
                        vv = v9 / cs
                        ww = w9 / cs
                        a110 = uu * vv; a101 = uu * ww; a011 = vv * ww
                        a200 = uu * uu; a020 = vv * vv; a002 = ww * ww  ! Second
                        a210 = a200 * vv
                        a201 = a200 * ww
                        a120 = a020 * uu
                        a021 = a020 * ww
                        a102 = a002 * uu
                        a012 = a002 * vv
                        a111 = a110 * ww ! Third
                        a210012p = a210 + a012
                        a102120p = a102 + a120
                        a021201p = a021 + a201
                        a210012m = a210 - a012
                        a102120m = a102 - a120
                        a021201m = a021 - a201
                        eu = H100(ip)*uu + H010(ip)*vv + H001(ip)*ww
                        sec = 2. * ( H110(ip) * a110 + H101(ip) * a101 + H011(ip) * a011 ) &
                            + H200(ip) * a200 + H020(ip) * a020 + H002(ip) * a002
                        thir = ( H210012p(ip) ) * ( a210012p ) + &
                               ( H102120p(ip) ) * ( a102120p ) + &
                               ( H021201p(ip) ) * ( a021201p ) + &
                               ( H210012m(ip) ) * ( a210012m ) / 3.d0 + &
                               ( H102120m(ip) ) * ( a102120m ) / 3.d0 + &
                               ( H021201m(ip) ) * ( a021201m ) / 3.d0
                        feq9(ip) = 1. + eu + 0.5 * ( sec + thir )
                        feq9(ip) = wwp(ip) * rho9 * feq9(ip)               
                !non-equilibrium distribution function & a^neq
                    a110 = rho9 * uu * vv; a101 = rho9 * uu * ww; a011 = rho9 * vv * ww
                    a200 = rho9 * uu * uu; a020 = rho9 * vv * vv; a002 = rho9 * ww * ww  ! Second   
                    ane110 = 0.0; ane101 = 0.0; ane011 = 0.0
                    ane200 = 0.0; ane020 = 0.0; ane002 = 0.0
                    ! f9 = f9(ip)
                    ! ane110 = ane110 + H110(ip) * f9
                    ! ane101 = ane101 + H101(ip) * f9
                    ! ane011 = ane011 + H011(ip) * f9
                    ! ane200 = ane200 + H200(ip) * f9
                    ! ane020 = ane020 + H020(ip) * f9
                    ! ane002 = ane002 + H002(ip) * f9

                    ane110 = -tau * rho9 * (dudy(i,j,k)+dvdx(i,j,k))
                    ane101 = 0.0 !-tau * rho9 * 0.5 * (0.0+0.0)   (dudz=0; dwdx=0)
                    ane011 = 0.0 !-tau * rho9 * 0.5 * (0.0+0.0)   (dvdz=0; dwdy=0)
                    ane200 = -tau * rho9 * 2.0 * dudx(i,j,k)
                    ane020 = -tau * rho9 * 2.0 * dvdy(i,j,k)
                    ane002 = 0.0
                    ane210 = ane200 * vv + 2. * uu * ane110
                    ane201 = ane200 * ww + 2. * uu * ane101
                    ane120 = ane020 * uu + 2. * vv * ane110
                    ane021 = ane020 * ww + 2. * vv * ane011
                    ane102 = ane002 * uu + 2. * ww * ane101
                    ane012 = ane002 * vv + 2. * ww * ane011
                    ane210012p = ane210 + ane012
                    ane102120p = ane102 + ane120
                    ane021201p = ane021 + ane201
                    ane210012m = ane210 - ane012
                    ane102120m = ane102 - ane120
                    ane021201m = ane021 - ane201
                    
                        wwp9 = wwp(ip)
                          
                          sec = 2. * ( H110(ip) * ane110 + H101(ip) * ane101 + H011(ip) * ane011 ) &
                              + H200(ip) * ane200 + H020(ip) * ane020 + H002(ip) * ane002
                          thir = ( H210012p(ip) ) * ( ane210012p ) + &
                               ( H102120p(ip) ) * ( ane102120p ) + &
                               ( H021201p(ip) ) * ( ane021201p ) + &
                               ( H210012m(ip) ) * ( ane210012m ) / 3.d0 + &
                               ( H102120m(ip) ) * ( ane102120m ) / 3.d0 + &
                               ( H021201m(ip) ) * ( ane021201m ) / 3.d0
                        
                          fneq9(ip) = wwp9 * 0.5 * ( sec + thir )
                          f(ip,i,j,k)= feq9(ip)+fneq9(ip)
                end do
            end do
        end do
    end do

    do it = 1, iit

        dudx = 0.0
        dudy = 0.0
        dvdx = 0.0
        dvdy = 0.0
        rhou = 0.0
        rhov = 0.0
        drhoudx = 0.0
        drhoudy = 0.0
        drhovdx = 0.0
        drhovdy = 0.0

        ! Compute rho*u and rho*v
        do k = 1, LZ
            do j = 1, L
                do i = 1, L
                    do ip = 1, ippQ
                        rhou(i,j,k) = rhou(i,j,k) + ix(ip) * f(ip,i,j,k)
                        rhov(i,j,k) = rhov(i,j,k) + iy(ip) * f(ip,i,j,k)
                        rhow(i,j,k) = rhow(i,j,k) + iz(ip) * f(ip,i,j,k)
                    end do
                end do
            end do
        end do

        ! Compute spatial derivatives
        do k = 1, LZ
            do j = 1, L
                do i = 1, L
                    ip1 = mod(i, L) + 1
                    ip2 = mod(i+1, L) + 1
                    im1 = mod(i-2+L, L) + 1 
                    im2 = mod(i-3+L, L) + 1
    
                    jp1 = mod(j, L) + 1
                    jp2 = mod(j+1, L) + 1
                    jm1 = mod(j-2+L, L) + 1
                    jm2 = mod(j-3+L, L) + 1
    
                    dudx(i,j,k) = - (u(ip2,j,k) - 8.0*u(ip1,j,k) + 8.0*u(im1,j,k) - u(im2,j,k)) / 12.0
                    dudy(i,j,k) = - (u(i,jp2,k) - 8.0*u(i,jp1,k) + 8.0*u(i,jm1,k) - u(i,jm2,k)) / 12.0
                    dvdx(i,j,k) = - (v(ip2,j,k) - 8.0*v(ip1,j,k) + 8.0*v(im1,j,k) - v(im2,j,k)) / 12.0
                    dvdy(i,j,k) = - (v(i,jp2,k) - 8.0*v(i,jp1,k) + 8.0*v(i,jm1,k) - v(i,jm2,k)) / 12.0
                    
                    drhoudx(i,j,k) = - (rhou(ip2,j,k) - 8.0*rhou(ip1,j,k) + 8.0*rhou(im1,j,k) - rhou(im2,j,k)) / 12.0
                    drhoudy(i,j,k) = - (rhou(i,jp2,k) - 8.0*rhou(i,jp1,k) + 8.0*rhou(i,jm1,k) - rhou(i,jm2,k)) / 12.0
                    drhovdx(i,j,k) = - (rhov(ip2,j,k) - 8.0*rhov(ip1,j,k) + 8.0*rhov(im1,j,k) - rhov(im2,j,k)) / 12.0
                    drhovdy(i,j,k) = - (rhov(i,jp2,k) - 8.0*rhov(i,jp1,k) + 8.0*rhov(i,jm1,k) - rhov(i,jm2,k)) / 12.0
                end do
            end do
        end do

        ! collision
        do k = 1, LZ
            do j = 1, L
                do i = 1, L
                    u9 = u(i,j,k)
                    v9 = v(i,j,k)
                    w9 = w(i,j,k)
                    rho9 = rho(i,j,k)
                    !equilibrium distribution function & a^eq
                    uu = u9 / cs
                    vv = v9 / cs
                    ww = w9 / cs
                    a110 = uu * vv; a101 = uu * ww; a011 = vv * ww
                    a200 = uu * uu; a020 = vv * vv; a002 = ww * ww  ! Second
                    a210 = a200 * vv
                    a201 = a200 * ww
                    a120 = a020 * uu
                    a021 = a020 * ww
                    a102 = a002 * uu
                    a012 = a002 * vv
                    a111 = a110 * ww ! Third
                    a210012p = a210 + a012
                    a102120p = a102 + a120
                    a021201p = a021 + a201
                    a210012m = a210 - a012
                    a102120m = a102 - a120
                    a021201m = a021 - a201
                    do ip = 1, ippQ
                        wwp9 = wwp(ip)
                        eu = H100(ip)*uu + H010(ip)*vv + H001(ip)*ww
                        sec = 2. * ( H110(ip) * a110 + H101(ip) * a101 + H011(ip) * a011 ) &
                            + H200(ip) * a200 + H020(ip) * a020 + H002(ip) * a002
                        thir = ( H210012p(ip) ) * ( a210012p ) + &
                               ( H102120p(ip) ) * ( a102120p ) + &
                               ( H021201p(ip) ) * ( a021201p ) + &
                               ( H210012m(ip) ) * ( a210012m ) / 3.d0 + &
                               ( H102120m(ip) ) * ( a102120m ) / 3.d0 + &
                               ( H021201m(ip) ) * ( a021201m ) / 3.d0
                        feq9(ip) = 1. + eu + 0.5 * ( sec + thir )             
                        feq9(ip) = wwp9 * rho9 * feq9(ip)
                    end do
                ! non-equilibrium distribution function & a^neq
                    uu = u9 / cs
                    vv = v9 / cs
                    ww = w9 / cs
                    a110 = rho9 * uu * vv; a101 = rho9 * uu * ww; a011 = rho9 * vv * ww
                    a200 = rho9 * uu * uu; a020 = rho9 * vv * vv; a002 = rho9 * ww * ww  ! Second   
                    ane110 = 0.0; ane101 = 0.0; ane011 = 0.0
                    ane200 = 0.0; ane020 = 0.0; ane002 = 0.0
                    do ip = 1, ippQ
                        f9 = f(ip,i,j,k)
                        ane110 = ane110 + H110(ip) * f9
                        ane101 = ane101 + H101(ip) * f9
                        ane011 = ane011 + H011(ip) * f9
                        ane200 = ane200 + H200(ip) * f9
                        ane020 = ane020 + H020(ip) * f9
                        ane002 = ane002 + H002(ip) * f9
                    end do
                    ane110 = (ane110 - a110)
                    ane101 = (ane101 - a101)
                    ane011 = (ane011 - a011)
                    ane200 = (ane200 - a200)
                    ane020 = (ane020 - a020)
                    ane002 = (ane002 - a002)  

                    ane210 = ane200 * vv + 2. * uu * ane110
                    ane201 = ane200 * ww + 2. * uu * ane101
                    ane120 = ane020 * uu + 2. * vv * ane110
                    ane021 = ane020 * ww + 2. * vv * ane011
                    ane102 = ane002 * uu + 2. * ww * ane101
                    ane012 = ane002 * vv + 2. * ww * ane011
                    ane210012p = ane210 + ane012
                    ane102120p = ane102 + ane120
                    ane021201p = ane021 + ane201
                    ane210012m = ane210 - ane012
                    ane102120m = ane102 - ane120
                    ane021201m = ane021 - ane201

                    do ip = 1, ippQ
                        wwp9 = wwp(ip)
                          sec = 2. * ( H110(ip) * ane110 + H101(ip) * ane101 + H011(ip) * ane011 ) &
                              + H200(ip) * ane200 + H020(ip) * ane020 + H002(ip) * ane002
                          thir = ( H210012p(ip) ) * ( ane210012p ) + &
                               ( H102120p(ip) ) * ( ane102120p ) + &
                               ( H021201p(ip) ) * ( ane021201p ) + &
                               ( H210012m(ip) ) * ( ane210012m ) / 3.d0 + &
                               ( H102120m(ip) ) * ( ane102120m ) / 3.d0 + &
                               ( H021201m(ip) ) * ( ane021201m ) / 3.d0
                          fneq9(ip) = wwp9 * 0.5 * ( sec + thir ) 
                    end do
                ! collision
                    do ip = 1, ippQ
                        fpost(ip,i,j,k) = feq9(ip) + (1.- 1./tau) * fneq9(ip)
                    end do
                end do
            end do
        end do
         ! Streaming
        do k = 1, LZ
            do j = 1, L
                do i = 1, L
                    do ip = 1, ippQ
                        ipre = i - ix(ip)
                        jpre = j - iy(ip)
                        kpre = k - iz(ip)
                        
                        if(ipre < 1) ipre = ipre + L
                        if(ipre > L) ipre = ipre - L
                        if(jpre < 1) jpre = jpre + L
                        if(jpre > L) jpre = jpre - L
                        if(kpre < 1) kpre = kpre + LZ
                        if(kpre > 2) kpre = kpre - LZ
                        f(ip,i,j,k) = fpost(ip,ipre,jpre,kpre)
                    end do
                    f(1,i,j,k) = fpost(1,i,j,k)  
                end do
            end do
        end do
        ! macrovarible
        do k = 1, LZ
            do j = 1, L
                do i = 1, L
                    rho9 = 0.0d0
                    u9 = 0.0d0
                    v9 = 0.0d0
                    w9 = 0.0d0
                    do ip = 1, ippQ
                        rho9 = rho9 + f(ip,i,j,k)
                        u9 = u9 + f(ip,i,j,k)*ix(ip)
                        v9 = v9 + f(ip,i,j,k)*iy(ip)
                        w9 = w9 + f(ip,i,j,k)*iz(ip)
                    end do
                    rho(i,j,k) = rho9
                    u(i,j,k) = u9 / rho9
                    v(i,j,k) = v9 / rho9
                    w(i,j,k) = w9 / rho9
                    p(i,j,k) = (rho9-1.0d0) * cs2
                    
                end do
            end do
        end do
        if (it == iit) then            
            ! compute error term with NS
            do k = 1, LZ
                do j = 1, L
                    do i = 1, L
                        ip1 = mod(i, L) + 1
                        ip2 = mod(i+1, L) + 1
                        im1 = mod(i-2+L, L) + 1 
                        im2 = mod(i-3+L, L) + 1
    
                        jp1 = mod(j, L) + 1
                        jp2 = mod(j+1, L) + 1
                        jm1 = mod(j-2+L, L) + 1
                        jm2 = mod(j-3+L, L) + 1
    
                        dudx(i,j,k) = - (u(ip2,j,k) - 8.0*u(ip1,j,k) + 8.0*u(im1,j,k) - u(im2,j,k)) / 12.0
                        dudy(i,j,k) = - (u(i,jp2,k) - 8.0*u(i,jp1,k) + 8.0*u(i,jm1,k) - u(i,jm2,k)) / 12.0
                        dvdx(i,j,k) = - (v(ip2,j,k) - 8.0*v(ip1,j,k) + 8.0*v(im1,j,k) - v(im2,j,k)) / 12.0
                        dvdy(i,j,k) = - (v(i,jp2,k) - 8.0*v(i,jp1,k) + 8.0*v(i,jm1,k) - v(i,jm2,k)) / 12.0
    
                        drhodx(i, j,k) = (-rho(ip2, j,k) + 8.0*rho(ip1, j,k) - 8.0*rho(im1, j,k) + rho(im2, j,k)) / 12.0 
                        drhody(i, j,k) = (-rho(i, jp2,k) + 8.0*rho(i, jp1,k) - 8.0*rho(i, jm1,k) + rho(i, jm2,k)) / 12.0 
                        udrhodx(i, j,k) = drhodx(i, j,k) * u(i,j,k)
                        vdrhody(i, j,k) = drhody(i, j,k) * v(i,j,k)
    
                        d2udxy(i,j,k) = (u(ip2, jp2,k) - 8.0 * u(ip1, jp2,k) + 8.0 * u(im1, jp2,k) - u(im2, jp2,k) &
                        - 8.0 * u(ip2, jp1,k) + 64.0 * u(ip1, jp1,k) - 64.0 * u(im1, jp1,k) + 8.0 * u(im2, jp1,k) &
                        + 8.0 * u(ip2, jm1,k) - 64.0 * u(ip1, jm1,k) + 64.0 * u(im1, jm1,k) - 8.0 * u(im2, jm1,k) &
                        - u(ip2, jm2,k) + 8.0 * u(ip1, jm2,k) - 8.0 * u(im1, jm2,k) + u(im2, jm2,k)) / 144.0
    
                        d2vdxy(i,j,k) = (v(ip2, jp2,k) - 8.0 * v(ip1, jp2,k) + 8.0 * v(im1, jp2,k) - v(im2, jp2,k) &
                        - 8.0 * v(ip2, jp1,k) + 64.0 * v(ip1, jp1,k) - 64.0 * v(im1, jp1,k) + 8.0 * v(im2, jp1,k) &
                        + 8.0 * v(ip2, jm1,k) - 64.0 * v(ip1, jm1,k) + 64.0 * v(im1, jm1,k) - 8.0 * v(im2, jm1,k) &
                        - v(ip2, jm2,k) + 8.0 * v(ip1, jm2,k) - 8.0 * v(im1, jm2,k) + v(im2, jm2,k)) /144.0
    
                        d2vdy2(i, j,k) = (-v(i, jp2,k) + 16.0*v(i, jp1,k) - 30.0*v(i, j,k) + 16.0*v(i, jm1,k) - v(i, jm2,k)) / 12.0 
                    end do
                end do
            end do

            open(unit=20, file=filename_v, status='replace', action='write')
            write(20, '(A)') '# X, Y, U, V, P,dUdx, dUdy, dVdx, dVdy, d2Udxy, d2Vdxy, Udrhodx, Vdrhody'
                do k = 1, LZ
                    do j = 1, L
                        do i = 1, L
                            write(20, '(13F12.6)') real(i - 0.5, kind=8), real(j - 0.5, kind=8), &
                                 u(i,j,1), v(i,j,1), p(i,j,1), &
                                 dudx(i,j,1), dudy(i,j,1), dvdx(i,j,1), dvdy(i,j,1), &
                                 d2udxy(i,j,1), d2vdxy(i,j,1), &
                                 udrhodx(i,j,1), vdrhody(i,j,1)
                        end do
                    end do
                end do
            close(unit=20)
        endif
    end do
   
    ! theory solution
    utheo = 0.0d0
    vtheo = 0.0d0
    ptheo = 0.0d0

    open(unit=30, file=filename_va, status='replace', action='write')
    write(30, '(A)') '# X, Y, U, V, P, d2Udxy, d2Vdxy'
    do k = 1,LZ
        do j = 1, L
            yy = real(j - 0.5, kind=8)
            do i = 1, L           
                xx = real(i - 0.5, kind=8)         
    
                utheo(i,j,k) = -U0*cos(k_val*pi*xx/L)*sin(k_val*pi*yy/L)*exp(-2.0d0*k_val**2*pi**2*visc*iit/L**2)
                vtheo(i,j,k) = U0*cos(k_val*pi*yy/L)*sin(k_val*pi*xx/L)*exp(-2.0d0*k_val**2*pi**2*visc*iit/L**2)
                ptheo(i,j,k) = (-U0*U0/2.0d0*cos(k_val*pi/dble(L)*(xx - yy)) &
                *cos(k_val*pi/dble(L)*(xx + yy)))*exp(-4.0d0*k_val**2*pi**2*visc*dble(iit)/dble(L)**2)
    
                d2udxytheo(i,j,k) = (k_val*pi/L)**2*U0*sin(k_val*pi*xx/L)*&
                                    cos(k_val*pi*yy/L)*exp(-2.0d0*k_val**2*pi**2*visc*iit/L**2)
                d2vdxytheo(i,j,k) = -(k_val*pi/L)**2*U0*sin(k_val*pi*yy/L)*&
                                    cos(k_val*pi*xx/L)*exp(-2.0d0*k_val**2*pi**2*visc*iit/L**2)
    
                write(30, '(7F12.6)') real(i - 0.5, kind=8), real(j - 0.5, kind=8), utheo(i,j,1), vtheo(i,j,1),&
                                        ptheo(i,j,1), d2udxytheo(i,j,1), d2vdxytheo(i,j,1)                                  
            end do
        end do
    end do
    call cpu_time(end_time)
    time = end_time - start_time

    ! midline solution comparison
    open(unit=10, file=filename_main, status='replace', action='write')
        write(10, *) 'Execution time: ', time, ' seconds'
        write(10, '(A)') '# Position, Numerical V, Theoretical V,P, Udrhodx, Vdrhody, mud2Vdy2, mud2Vdxy, dvdx, dvdy'
        
        jj1 = L/2; jj2 = jj1 + 1;  
        
        do i = 1, L
            xx = real(i - 0.5, kind=8)
            
            vnum(i) = 0.5 * (v(i,jj1,1) + v(i,jj2,1)) 
            vtheoline(i) = -U0*sin(k_val*pi*xx/L)*exp(-2.0d0*k_val**2*pi**2*visc*iit/L**2)
            pnum(i) = 0.5 * (p(i,jj1,1) + p(i,jj2,1))
            udrhodxnum(i) = 0.5 * (udrhodx(i,jj1,1) + udrhodx(i,jj2,1))
            vdrhodynum(i) = 0.5 * (vdrhody(i,jj1,1) + vdrhody(i,jj2,1))
            d2vdy2num(i) = visc * 0.5 * (d2vdy2(i,jj1,1) + d2vdy2(i,jj2,1)) 
            d2vdxynum(i) = visc * 0.5 * (d2vdxy(i,jj1,1) + d2vdxy(i,jj2,1)) 
            dvdxnum(i) =  0.5 * (dvdx(i,jj1,1) + dvdx(i,jj2,1))
            dvdynum(i) =  0.5 * (dvdy(i,jj1,1) + dvdy(i,jj2,1))

            write(10, '(I5, 10F20.10)') i, vnum(i), vtheoline(i), pnum(i), udrhodxnum(i), vdrhodynum(i),&
                                       d2vdy2num(i), d2vdxynum(i), dvdxnum(i), dvdynum(i)
        end do
        
    close(unit=10)
    ! write(*,*) u(50,50,1)

    deallocate(f, fpost, feq, rho, u, v, p, utheo, vtheo, ptheo)
    deallocate(feq9, fneq9)

    deallocate(dudxi, dudyi, dvdxi, dvdyi, rhoui, rhovi, drhoudxi, drhoudyi, drhovdxi, drhovdyi)
    deallocate(dudx, dudy, dvdx, dvdy, d2udxy, d2vdxy, d2vdy2, d2udxytheo, d2vdxytheo)
    deallocate(rhou, rhov, drhoudx, drhoudy, drhovdx, drhovdy)
    deallocate(drhodx, drhody, udrhodx, vdrhody)
    end subroutine run_simulation

end program Taylor