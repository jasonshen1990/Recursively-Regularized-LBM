! -------------------------------------------------------------------------
!       F90 Code for Simulating lid-driven cavity flow
!       Using LBM D3Q27 RR Model.
! -------------------------------------------------------------------------
! this code is identical with test23 code, except that new subroutines
! for statistics of particle and fluid rms velocities are added 
      PROGRAM main
      use var_inc
      implicit none
      integer:: i,j,k

      call para
      call allocarray
      call Hermite_polynomial

        ux = 0.0
        uy = 0.0
        uz = 0.0
        rho = 1.0

        call initpop

      call FORCING
! main loop

      DO while(istep .lt. 100000)!istep0+nsteps
          if(mod(istep,100).eq.0) then
             write(*,*) 'istep= ',istep+1, uy(128, 64,1), uy(128, 64, 2), uz(128, 64, 2)
          end if
    
          call collision_RR
          call streaming
          call macrovar
    
          istep = istep+1
    
          if(mod(istep,5000).eq.0)  call outputflow
      END DO

      call outputflow
! main loop ends here  

101   continue

      END PROGRAM main

!*******************outputflow*************
      subroutine outputflow
      use var_inc
      implicit none

      character(len=200) :: fnm,chartemp
      character(len=256) :: command
      logical :: dir_exists

      ! Check if output directory exists, create it if not (using SYSTEM)
      dir_exists = .false.
      call system('if exist 32re1000 (set dir_exists=true)')

      if (.not. dir_exists) then
          print *, 'Creating output directory...'
          call system('mkdir 32re1000')
      end if

       write(chartemp,'(i6)') istep
       fnm = '32re1000/rho'//trim(adjustl(chartemp))//'.dat'
       open(1,file=trim(fnm),status = 'unknown',form = 'unformatted')
       write(1) rho(:,:,:)
       close(1)

       fnm = '32re1000/ux'//trim(adjustl(chartemp))//'.dat'
       open(1,file=trim(fnm),status = 'unknown',form = 'unformatted')
       write(1) ux(:,:,:)
       close(1)
       
       fnm = '32re1000/uy'//trim(adjustl(chartemp))//'.dat'
       open(1,file=trim(fnm),status = 'unknown',form = 'unformatted')
       write(1) uy(:,:,:)
       close(1)

       fnm = '32re1000/uz'//trim(adjustl(chartemp))//'.dat'
       open(1,file=trim(fnm),status = 'unknown',form = 'unformatted')
       write(1) uz(:,:,:)
       close(1)

160  format(2x,8(1pe16.6))

      return
      end
