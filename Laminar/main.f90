!=============================================================================
!       MPI F90 Code for Simulating turbulent flows over canopy
!
!       Using LBM D3Q27 MRT Model.
!
!       MPI version is designed to operate on NK's clusters
!
!       The DHIT code was originally written by Jason, Dec 2024.
!
!=============================================================================
    PROGRAM main
    ! use mpi
    use var_inc
    implicit none
!      real xx,yy
    integer :: i, j, k 
    integer :: ix,iy,iz
	real :: start_time, end_time, time

    call cpu_time(start_time)

    !Initialize all variables used
    !@file para.f90
    call para
    !Allocate all global arrays used in var_inc module
    !@file para.f90
    call allocarray
    !Construct custom MPI data types
    !@file para.f90
    call Hermite_polynomial

!=======================================================
! FLOW INITIALIZATION/ LOADING
!=======================================================
         ux = 0.0
         uy = 0.0
         uz = 0.0
    ! Calculate forcing values used in collision_MRT
    !@file collision.f90
    call FORCING

    !Initialize particle populations on each node
    !@file partlib.f90
    call initpop

    istep = 0

    ! save initial pre-relaxed flow      
    call macrovar !@file collison.f90
    istep0 = 0
    istep = istep0

    if(myid.eq.0)write(*,*)'Loading successful'

!=======================================================
! END OF FLOW INITIALIZATION/ LOADING
!=======================================================
    ! Calculate forcing values used in collision_MRT
    ! Forcing is independent of time so only call once here
    !@file collision.f90
    ! call FORCING

    !Create building boundary
    !@file building.f90
    ! call building

    !Call macrovar to initialize velocity variables
    call macrovar !@file collision.f90

! MAIN LOOP
!=======================================================
    do istep = istep0+1,istep0+nsteps 

        if(myid.eq.0 .and. mod(istep,1000).eq.0)then
            write(*,*) 'istep=',istep
        endif

        call collision_RR
            
        ! call streaming
        call macrovar
        if (mod(istep,100).eq.0) write(*,*) 'uy',istep, uy(nx/2,ny/2,nz/2)/ustar
        if (mod(istep,100).eq.0) write(*,*) 'uybotto',istep, uy(1,ny/2,nz/2)/ustar
        if (mod(istep,100).eq.0) write(*,*) 'uytoppp',istep, uy(nx,ny/2,nz/2)/ustar
        if (mod(istep,100).eq.0) write(*,*) 'checkrho',istep, rho(nx/2,ny/2,nz/2)

        if (mod(istep, 5000).eq.0) call outputflow

        ! if (istep == nsteps) then
        !     close(output_unit)
        !     file_opened = .false.
        ! endif

    end do
!=======================================================
! END OF MAIN LOOP
!=======================================================
    call cpu_time(end_time)
    time = end_time - start_time
    write(*,*) time

    END PROGRAM main

    subroutine outputflow
    use var_inc
    implicit none

    character(len=160) :: fnm, chartemp
    character(len=256) :: command
    logical :: dir_exists

    ! Check if output directory exists, create it if not (using SYSTEM)
    dir_exists = .false.
    call system('if exist output (set dir_exists=true)')

    if (.not. dir_exists) then
        print *, 'Creating output directory...'
        call system('mkdir output')
    end if

    ! Write ux data to file
    write(chartemp, '(i8)') istep
    fnm = 'output/ux'//trim(adjustl(chartemp))//'.dat'
    open(1, file=trim(fnm), status = 'unknown', form='unformatted')
    write(1) ux
    close(1)

    ! Write uy data to file
    write(chartemp, '(i8)') istep
    fnm = 'output/uy'//trim(adjustl(chartemp))//'.dat'
    open(1, file=trim(fnm), status = 'unknown', form='unformatted')
    write(1) uy
    close(1)
  
    ! Write uz data to file
    write(chartemp, '(i8)') istep
    fnm = 'output/uz'//trim(adjustl(chartemp))//'.dat'
    open(1, file=trim(fnm), status = 'unknown', form='unformatted')
    write(1) uz
    close(1)

    return
end subroutine outputflow

