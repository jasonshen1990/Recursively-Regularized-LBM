!=============================================================================
!@subroutine collision_RR
!@desc Executes collision and propagation of the fluid. Current follows the
!      two-array algorythm which fuses collision and propagation together
!      into a single step. Calculations have be movef to element-wise
!      operations to help increase performance. May change in the future.
!=============================================================================
      subroutine collision_RR  
      ! use mpi
      use var_inc
      implicit none

      real, dimension(0:npop-1) :: f9 , fneqRR
      real, dimension(0:npop-1) :: Fbar
      real, dimension(0:npop-1,lx,ly)     :: tmpzpS, tmpzmS
      real, dimension(0:npop-1,lx,0:lz+1) :: tmpypS,tmpymS
      real rho9, ux9, uy9, uz9, ux9s, uy9s, uz9s, u9,v9,w9
      integer ip, jp, ix, iy, iz, imove, jmove, kmove, ipi, i ,j ,k
      real stress, stress1, stress2, stress3, stress_xy, stress_xy1, stress_xy2, stress_xy3, stresswalla1, stresswallb
      real ss, ss1, stheory, stheta, stresswalla2, theta0, ycc, theory, error,theta
      real fx9,fy9,fz9,G1,G2,G3
      real*8 feq, feqRR, feqRR_4
      real fchange
      integer :: output_unit, iconf, iout
      real fneq7,fneq8,fneq9,fneq10,fneq19,fneq20,fneq21,fneq22,fneq23,fneq24,fneq25,fneq26

      !Index through fluid nodes
      do iz = 1,nz
      do iy = 1,ny
      do ix = 1,nx
      !If we are not in a solid particle execute collision
      ! if(ibnodes(ix,iy,iz) > 0)goto 111

        f9(:) = f(:,ix,iy,iz)
        rho9 = rho(ix,iy,iz)
        ux9 = ux(ix,iy,iz)
        uy9 = uy(ix,iy,iz)
        uz9 = uz(ix,iy,iz)
        ux9s = ux9*ux9
        uy9s = uy9*uy9
        uz9s = uz9*uz9

        !Merging the forcing term with the collision operator
        fx9 = force_realx(ix,iy,iz)
        fy9 = force_realy(ix,iy,iz)
        fz9 = force_realz(ix,iy,iz)
        G3 = ux9*fx9 + uy9*fy9 + uz9*fz9

        Fbar(0) = -0.5*G3/RT

        do ip=1,npop-1
          G1 = cix(ip)*fx9 + ciy(ip)*fy9 + ciz(ip)*fz9
          Fbar(ip) = 0.5*(G1 - G3)/RT
        enddo

        fneqRR = 0.0
        call comptfneqRR(rho9, ux9, uy9, uz9, f9, Fbar, fneqRR)
        do ip = 0, npop-1
          fchange = (1.- 1./tau) * fneqRR(ip) &
          + Fbar(ip)*feqRR_4(ip, rho9, ux9, uy9, uz9)
          f9(ip) = feqRR(ip,rho9,ux9,uy9,uz9) + fchange
        enddo
        
        do ipi = 1,14
          ip = ipstay(ipi)  
          imove = ix + cix(ip) 
          jmove = iy + ciy(ip)
          kmove = iz + ciz(ip)
          if(imove < 1 .or. imove > lx)then
            f(ipopp(ip),ix,iy,iz) = f9(ip) !+ 0.5*Fbar(ip)
          elseif(jmove < 1)then
             tmpymS(ip,imove,kmove) =  f9(ip) !+ 0.5*Fbar(ip)
          elseif(jmove > ly)then
            tmpypS(ip,imove,kmove) = f9(ip) !+ 0.5*Fbar(ip)
          elseif(kmove < 1)then
            tmpzmS(ip,imove,jmove) = f9(ip) !+ 0.5*Fbar(ip)
          elseif(kmove > lz)then
            tmpzpS(ip,imove,jmove) = f9(ip) !+ 0.5*Fbar(ip)
          else
            f(ipopp(ip),ix,iy,iz) = f9(ip) !+ 0.5*Fbar(ip)
          endif
        enddo

! 111     continue        
        do ipi = 1,13
          ip = ipswap(ipi)
          imove = ix + cix(ip) 
          jmove = iy + ciy(ip)
          kmove = iz + ciz(ip)
          if(imove < 1 .or. imove > lx)then
            f(ipopp(ip),ix,iy,iz) = f9(ip) !+ 0.5*Fbar(ip)
          elseif(jmove < 1)then
             tmpymS(ip,imove,kmove) =  f9(ip) !+ 0.5*Fbar(ip)
          elseif(jmove > ly)then
            tmpypS(ip,imove,kmove) = f9(ip) !+ 0.5*Fbar(ip)
          elseif(kmove < 1)then
            tmpzmS(ip,imove,jmove) = f9(ip) !+ 0.5*Fbar(ip)
          elseif(kmove > lz)then
            tmpzpS(ip,imove,jmove) = f9(ip) !+ 0.5*Fbar(ip)
          else
            f(ipopp(ip),ix,iy,iz) = f(ip,imove,jmove,kmove)
            f(ip,imove,jmove,kmove) = f9(ip) !+ 0.5*Fbar(ip)
          endif
        enddo

      end do !x
      end do !y
      end do !z

      call collisionExchnge(tmpymS,tmpypS,tmpzmS,tmpzpS)

      DO iz = 1,lz
      DO iy = 1,ly
      DO ix = 1,1
          f9 = f(:, ix, iy, iz)
          rho9 = f9(0)+f9(1)+f9(2)+f9(3)+f9(4)+f9(5)+f9(6)+f9(7)+f9(8)  &
                +f9(9)+f9(10)+f9(11)+f9(12)+f9(13)+f9(14)+f9(15)+f9(16) &
                +f9(17)+f9(18)+f9(19)+f9(20)+f9(21)+f9(22)+f9(23)+f9(24) &
                +f9(25)+f9(26)
          u9 = f9(1)-f9(2)+f9(7)+f9(8)-f9(9)-f9(10)+f9(11)+f9(12)-f9(13) &
               -f9(14)+f9(19)+f9(20)+f9(21)+f9(22)-f9(23)-f9(24)-f9(25)-f9(26)
          v9 = f9(3)-f9(4)+f9(7)-f9(8)+f9(9)-f9(10)+f9(15)+f9(16)-f9(17) &
               -f9(18)+f9(19)+f9(20)-f9(21)-f9(22)+f9(23)+f9(24)-f9(25)-f9(26)
          w9 = f9(5)-f9(6)+f9(11)-f9(12)+f9(13)-f9(14)+f9(15)-f9(16)+   &
               f9(17)-f9(18)+f9(19)-f9(20)+f9(21)-f9(22)+f9(23)-f9(24) +  &
               f9(25)-f9(26)
          u9 = (u9+ force_realx(ix,iy,iz)/2.) / rho9
          v9 = (v9+ force_realy(ix,iy,iz)/2.) / rho9
          w9 = (w9+ force_realz(ix,iy,iz)/2.) / rho9
          fneqRR = 0.0
          call comptfneqRR(rho9, u9, v9, w9, f9,Fbar, fneqRR)
!--------Perform collision--------
          do ip = 0, npop-1
            if (cix(ip).gt.0) then
            f(ip,ix,iy,iz) = feqRR(ip,rho9,u9,v9,w9) + fneqRR(ip) - Fbar(ip)*feqRR_4(ip,rho9,u9,v9,w9)
            end if
          enddo
      END DO
      END DO
      END DO

      DO iz = 1,lz
      DO iy = 1,ly
      DO ix = lx,lx
          f9 = f(:, ix, iy, iz)
          rho9 = f9(0)+f9(1)+f9(2)+f9(3)+f9(4)+f9(5)+f9(6)+f9(7)+f9(8)  &
                +f9(9)+f9(10)+f9(11)+f9(12)+f9(13)+f9(14)+f9(15)+f9(16) &
                +f9(17)+f9(18)+f9(19)+f9(20)+f9(21)+f9(22)+f9(23)+f9(24) &
                +f9(25)+f9(26)
          u9 = f9(1)-f9(2)+f9(7)+f9(8)-f9(9)-f9(10)+f9(11)+f9(12)-f9(13) &
               -f9(14)+f9(19)+f9(20)+f9(21)+f9(22)-f9(23)-f9(24)-f9(25)-f9(26)
          v9 = f9(3)-f9(4)+f9(7)-f9(8)+f9(9)-f9(10)+f9(15)+f9(16)-f9(17) &
               -f9(18)+f9(19)+f9(20)-f9(21)-f9(22)+f9(23)+f9(24)-f9(25)-f9(26)
          w9 = f9(5)-f9(6)+f9(11)-f9(12)+f9(13)-f9(14)+f9(15)-f9(16)+   &
               f9(17)-f9(18)+f9(19)-f9(20)+f9(21)-f9(22)+f9(23)-f9(24) +  &
               f9(25)-f9(26)
          u9 = (u9+ force_realx(ix,iy,iz)/2.) / rho9
          v9 = (v9+ force_realy(ix,iy,iz)/2.) / rho9
          w9 = (w9+ force_realz(ix,iy,iz)/2.) / rho9
          fneqRR = 0.0
          call comptfneqRR(rho9, u9, v9, w9, f9,Fbar, fneqRR)
!--------Perform collision--------
          do ip = 0, npop-1
            if (cix(ip).lt.0) then
            f(ip,ix,iy,iz) = feqRR(ip,rho9,u9,v9,w9) + fneqRR(ip) - Fbar(ip)*feqRR_4(ip,rho9,u9,v9,w9)
            end if
          enddo
      END DO
      END DO
      END DO

      ! call collisionExchnge(tmpymS,tmpypS,tmpzmS,tmpzpS)
      !Exchange information with MPI neighbors to update edges

      end subroutine collision_RR
!=============================================================================
!@subroutine collisionExchnge 
!@desc Exchanges updated fluid distributions between neighboring MPI tasks and
!      updates the local domains edges with distributions recieved
!@param tmpymSi,tmpypSi,tmpzmSi,tmpzpSi = real array; contains fluid
!       distributions that need to be passed to nieghboring MPI tasks
!=============================================================================
      subroutine collisionExchnge(tmpymSi,tmpypSi,tmpzmSi,tmpzpSi)
      use var_inc
      ! use mpi
      implicit none
! 
      integer ilen
      ! integer status_array(MPI_STATUS_SIZE,4), req(4)
      real, dimension(0:npop-1,lx,0:lz+1) ::  tmpypSi, tmpymSi
      real, dimension(0:npop-1,lx,ly)     :: tmpzpSi, tmpzmSi
      real, dimension(9,lx,0:lz+1):: tmpypS, tmpymS, tmpypR, tmpymR
      real, dimension(9,lx,ly):: tmpzpS, tmpzmS, tmpzmR, tmpzpR
      real, dimension(0:npop-1,ly,lz)  :: tmpxB, tmpxT
      real*8 rho9, u9, v9, w9, feq, feqRR
      integer ix, iy, iz, ip
      real, dimension(0:npop-1) :: f9, fneqRR
! 
      !Create send buffers
      !Note that we only pass the distributions needed to decrease message size
      !from y out to y in
      tmpypS(1,:,:) = tmpypSi(3,:,:)
      tmpypS(2,:,:) = tmpypSi(7,:,:)
      tmpypS(3,:,:) = tmpypSi(9,:,:)
      tmpypS(4,:,:) = tmpypSi(15,:,:)
      tmpypS(5,:,:) = tmpypSi(16,:,:)
      tmpypS(6,:,:) = tmpypSi(19,:,:)
      tmpypS(7,:,:) = tmpypSi(20,:,:)
      tmpypS(8,:,:) = tmpypSi(23,:,:)
      tmpypS(9,:,:) = tmpypSi(24,:,:)
! 
      !from y in to y out
      tmpymS(1,:,:) = tmpymSi(4,:,:)
      tmpymS(2,:,:) = tmpymSi(8,:,:)
      tmpymS(3,:,:) = tmpymSi(10,:,:)
      tmpymS(4,:,:) = tmpymSi(17,:,:)
      tmpymS(5,:,:) = tmpymSi(18,:,:)
      tmpymS(6,:,:) = tmpymSi(21,:,:)
      tmpymS(7,:,:) = tmpymSi(22,:,:)
      tmpymS(8,:,:) = tmpymSi(25,:,:)
      tmpymS(9,:,:) = tmpymSi(26,:,:)
! 
      !Send/Recieve updated distributions with Y MPI neighbors
      tmpymR = tmpypS
      tmpypR = tmpymS
! 
      !Update local domain Y edge nodes
      !Note we must account for wall bounce back here!
      f(3,:,1,:) = tmpymR(1,:,1:lz)
      f(7,2:lx,1,:) = tmpymR(2,2:lx,1:lz)
      f(9,1:lx-1,1,:) = tmpymR(3,1:lx-1,1:lz)
      f(15,:,1,:) = tmpymR(4,:,1:lz)
      f(16,:,1,:) = tmpymR(5,:,1:lz)
      f(19,2:lx,1,2:lz) = tmpymR(6,2:lx,2:lz)
      f(20,2:lx,1,1:lz-1) = tmpymR(7,2:lx,1:lz-1)
      f(23,1:lx-1,1,2:lz) = tmpymR(8,1:lx-1,2:lz)
      f(24,1:lx-1,1,1:lz-1) = tmpymR(9,1:lx-1,1:lz-1)
! 
      ! 9, 24, 23 @ 1 in y direction, mirror reflection
! 
      f(4,:,ly,:) = tmpypR(1,:,1:lz)
      f(8,2:lx,ly,:) = tmpypR(2,2:lx,1:lz)
      f(10,1:lx-1,ly,:) = tmpypR(3,1:lx-1,1:lz)
      f(17,:,ly,:) = tmpypR(4,:,1:lz)
      f(18,:,ly,:) = tmpypR(5,:,1:lz)
      f(21,2:lx,ly,2:lz) = tmpypR(6,2:lx,2:lz)
      f(22,2:lx,ly,1:lz-1) = tmpypR(7,2:lx,1:lz-1)
      f(25,1:lx-1,ly,2:lz) = tmpypR(8,1:lx-1,2:lz)
      f(26,1:lx-1,ly,1:lz-1) = tmpypR(9,1:lx-1,1:lz-1)
! 
      ! 10, 25, 26 @ ly in y direction, mirror reflection
      ! 
      !Add corner distributions to Z buffer
      tmpzmSi(16,:,1) = tmpymR(5,:,0)
      tmpzmSi(20,:,1) = tmpymR(7,:,0)
      tmpzmSi(24,:,1) = tmpymR(9,:,0)
! 
      tmpzmSi(18,:,ly) = tmpypR(5,:,0)
      tmpzmSi(22,:,ly) = tmpypR(7,:,0)
      tmpzmSi(26,:,ly) = tmpypR(9,:,0)
! 
! 
      tmpzpSi(15,:,1) = tmpymR(4,:,lz+1)
      tmpzpSi(19,:,1) = tmpymR(6,:,lz+1)
      tmpzpSi(23,:,1) = tmpymR(8,:,lz+1)
! 
      tmpzpSi(17,:,ly) = tmpypR(4,:,lz+1)
      tmpzpSi(21,:,ly) = tmpypR(6,:,lz+1)
      tmpzpSi(25,:,ly) = tmpypR(8,:,lz+1)
! 
      !Add local data to Z send buffers
      tmpzpS(1,:,:) = tmpzpSi(5,:,:)
      tmpzpS(2,:,:) = tmpzpSi(11,:,:)
      tmpzpS(3,:,:) = tmpzpSi(13,:,:)
      tmpzpS(4,:,:) = tmpzpSi(15,:,:)
      tmpzpS(5,:,:) = tmpzpSi(17,:,:)
      tmpzpS(6,:,:) = tmpzpSi(19,:,:)
      tmpzpS(7,:,:) = tmpzpSi(21,:,:)
      tmpzpS(8,:,:) = tmpzpSi(23,:,:)
      tmpzpS(9,:,:) = tmpzpSi(25,:,:)
! 
      tmpzmS(1,:,:) = tmpzmSi(6,:,:)
      tmpzmS(2,:,:) = tmpzmSi(12,:,:)
      tmpzmS(3,:,:) = tmpzmSi(14,:,:)
      tmpzmS(4,:,:) = tmpzmSi(16,:,:)
      tmpzmS(5,:,:) = tmpzmSi(18,:,:)
      tmpzmS(6,:,:) = tmpzmSi(20,:,:)
      tmpzmS(7,:,:) = tmpzmSi(22,:,:)
      tmpzmS(8,:,:) = tmpzmSi(24,:,:)
      tmpzmS(9,:,:) = tmpzmSi(26,:,:)
! 
      tmpzmR = tmpzpS
      tmpzpR = tmpzmS
! 
      !Update local domain Z edge nodes
      !Note we must account for wall bounce back here!
      f(5,:,:,1) = tmpzmR(1,:,:)
      f(11,2:lx,:,1) = tmpzmR(2,2:lx,:)
      f(13,1:lx-1,:,1) = tmpzmR(3,1:lx-1,:)
      f(15,:,:,1) = tmpzmR(4,:,:)
      f(17,:,:,1) = tmpzmR(5,:,:)
      f(19,2:lx,:,1) = tmpzmR(6,2:lx,:)
      f(21,2:lx,:,1) = tmpzmR(7,2:lx,:)
      f(23,1:lx-1,:,1) = tmpzmR(8,1:lx-1,:)
      f(25,1:lx-1,:,1) = tmpzmR(9,1:lx-1,:)
! 
      ! 13, 15, 25 @ 1 in z direction, mirror reflection
! 
      f(6,:,:,lz) = tmpzpR(1,:,:)
      f(12,2:lx,:,lz) = tmpzpR(2,2:lx,:)
      f(14,1:lx-1,:,lz) = tmpzpR(3,1:lx-1,:)
      f(16,:,:,lz) = tmpzpR(4,:,:)
      f(18,:,:,lz) = tmpzpR(5,:,:)
      f(20,2:lx,:,lz) = tmpzpR(6,2:lx,:)
      f(22,2:lx,:,lz) = tmpzpR(7,2:lx,:)
      f(24,1:lx-1,:,lz) = tmpzpR(8,1:lx-1,:)
      f(26,1:lx-1,:,lz) = tmpzpR(9,1:lx-1,:)
! 
      ! 14, 26, 24 @ lz in z direction, mirror reflection

      ! bottom
!       f(1,  1, :, :)  = tmpxB(2,  :, :)
!       f(7,  1, :, :)  = tmpxB(10, :, :)
!       f(8,  1, :, :)  = tmpxB(9,  :, :)
!       f(11, 1, :, :)  = tmpxB(14, :, :)
!       f(12, 1, :, :)  = tmpxB(13, :, :)
!       f(19, 1, :, :)  = tmpxB(26, :, :)
!       f(20, 1, :, :)  = tmpxB(25, :, :)
!       f(21, 1, :, :)  = tmpxB(24, :, :)
!       f(22, 1, :, :)  = tmpxB(23, :, :)
!       DO iz = 1,lz
!       DO iy = 1,ly
!       DO ix = 1,1
!           f9 = f(:, ix, iy, iz)
!           rho9 = f9(0)+f9(1)+f9(2)+f9(3)+f9(4)+f9(5)+f9(6)+f9(7)+f9(8)  &
!                 +f9(9)+f9(10)+f9(11)+f9(12)+f9(13)+f9(14)+f9(15)+f9(16) &
!                 +f9(17)+f9(18)+f9(19)+f9(20)+f9(21)+f9(22)+f9(23)+f9(24) &
!                 +f9(25)+f9(26)
!           u9 = f9(1)-f9(2)+f9(7)+f9(8)-f9(9)-f9(10)+f9(11)+f9(12)-f9(13) &
!                -f9(14)+f9(19)+f9(20)+f9(21)+f9(22)-f9(23)-f9(24)-f9(25)-f9(26)
!           v9 = f9(3)-f9(4)+f9(7)-f9(8)+f9(9)-f9(10)+f9(15)+f9(16)-f9(17) &
!                -f9(18)+f9(19)+f9(20)-f9(21)-f9(22)+f9(23)+f9(24)-f9(25)-f9(26)
!           w9 = f9(5)-f9(6)+f9(11)-f9(12)+f9(13)-f9(14)+f9(15)-f9(16)+   &
!                f9(17)-f9(18)+f9(19)-f9(20)+f9(21)-f9(22)+f9(23)-f9(24) +  &
!                f9(25)-f9(26)
!           u9 = (u9+ force_realx(ix,iy,iz)/2.) / rho9
!           v9 = (v9+ force_realy(ix,iy,iz)/2.) / rho9
!           w9 = (w9+ force_realz(ix,iy,iz)/2.) / rho9
!           fneqRR = 0.0
!           call comptfneqRR(rho9, u9, v9, w9, f9, fneqRR)
! !--------Perform collision--------
!           do ip = 0, npop-1
!             if (cix(ip).gt.0) then
!             f9(ip) = feqRR(ip,rho9,u9,v9,w9) + fneqRR(ip)
!             end if
!           enddo
!       END DO
!       END DO
!       END DO
! !top
!       f(2,  lx, :, :)  = tmpxT(1,  :, :)
!       f(9,  lx, :, :)  = tmpxT(8,  :, :)
!       f(10, lx, :, :)  = tmpxT(7,  :, :)
!       f(13, lx, :, :)  = tmpxT(12, :, :)
!       f(14, lx, :, :)  = tmpxT(11, :, :)
!       f(23, lx, :, :)  = tmpxT(22, :, :)
!       f(24, lx, :, :)  = tmpxT(21, :, :)
!       f(25, lx, :, :)  = tmpxT(20, :, :)
!       f(26, lx, :, :)  = tmpxT(19, :, :)
!       DO iz = 1,lz
!       DO iy = 1,ly
!       DO ix = lx,lx
!           f9 = f(:, ix, iy, iz)
!           rho9 = f9(0)+f9(1)+f9(2)+f9(3)+f9(4)+f9(5)+f9(6)+f9(7)+f9(8)  &
!                 +f9(9)+f9(10)+f9(11)+f9(12)+f9(13)+f9(14)+f9(15)+f9(16) &
!                 +f9(17)+f9(18)+f9(19)+f9(20)+f9(21)+f9(22)+f9(23)+f9(24) &
!                 +f9(25)+f9(26)
!           u9 = f9(1)-f9(2)+f9(7)+f9(8)-f9(9)-f9(10)+f9(11)+f9(12)-f9(13) &
!                -f9(14)+f9(19)+f9(20)+f9(21)+f9(22)-f9(23)-f9(24)-f9(25)-f9(26)
!           v9 = f9(3)-f9(4)+f9(7)-f9(8)+f9(9)-f9(10)+f9(15)+f9(16)-f9(17) &
!                -f9(18)+f9(19)+f9(20)-f9(21)-f9(22)+f9(23)+f9(24)-f9(25)-f9(26)
!           w9 = f9(5)-f9(6)+f9(11)-f9(12)+f9(13)-f9(14)+f9(15)-f9(16)+   &
!                f9(17)-f9(18)+f9(19)-f9(20)+f9(21)-f9(22)+f9(23)-f9(24) +  &
!                f9(25)-f9(26)
!           u9 = (u9+ force_realx(ix,iy,iz)/2.) / rho9
!           v9 = (v9+ force_realy(ix,iy,iz)/2.) / rho9
!           w9 = (w9+ force_realz(ix,iy,iz)/2.) / rho9
!           fneqRR = 0.0
!           call comptfneqRR(rho9, u9, v9, w9, f9, fneqRR)
! !--------Perform collision--------
!           do ip = 0, npop-1
!             if (cix(ip).lt.0) then
!             f9(ip) = feqRR(ip,rho9,u9,v9,w9) + fneqRR(ip)
!             end if
!           enddo
!       END DO
!       END DO
!       END DO

      end subroutine
!=============================================================================
!@subroutine macrovar 
!@desc Calculates macroscopic fluid properties density(rho), x-velocity(ux),
!      y-velocity(uy), and z-velocity(uz)
!=============================================================================
      subroutine macrovar 
      use var_inc
      implicit none 

      integer ip, i, j, k, id

      real  sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum9,sum10,ux9,uy9,uz9
      real  rho9
      real, dimension(0:npop-1) :: f9

      !Calculate denisty and velocities of the fluid
      do k = 1,lz
      do j = 1,ly
      do i = 1,lx
          !If we are not inside a solid particle
          f9 = f(:,i,j,k)

          rho9 = f9(0)+f9(1)+f9(2)+f9(3)+f9(4)+f9(5)+f9(6)+f9(7)+f9(8)  &
                +f9(9)+f9(10)+f9(11)+f9(12)+f9(13)+f9(14)+f9(15)+f9(16) &
                +f9(17)+f9(18)+f9(19)+f9(20)+f9(21)+f9(22)+f9(23)+f9(24) &
                +f9(25)+f9(26)

          sum1 = f9(7) - f9(10)
          sum2 = f9(8) - f9(9)

          sum3 = f9(11) - f9(14)
          sum4 = f9(12) - f9(13)

          sum5 = f9(15) - f9(18)
          sum6 = f9(17) - f9(16)
          
          sum7 = f9(19) - f9(26)
          sum8 = f9(20) - f9(25)
       
          sum9 = f9(22) - f9(23)
          sum10 = f9(21) - f9(24)

          ux9 = f9(1)-f9(2)+sum1+sum2+sum3+sum4+sum7+sum8+sum9+sum10
          uy9 = f9(3)-f9(4)+sum1-sum2+sum5-sum6+sum7+sum8-sum9-sum10
          uz9 = f9(5)-f9(6)+sum3-sum4+sum5+sum6+sum7-sum8-sum9+sum10

          rho(i,j,k) = rho9
          ! print*, istep,rho9
          
          ! ux(i,j,k) = ux9 / rho9
          ! uy(i,j,k) = (uy9 + 0.5*rho9*force_in_y) / rho9
          ! uz(i,j,k) = uz9 /rho9

          ux(i,j,k) = (ux9 + force_realx(i,j,k)/2.) / rho9
          uy(i,j,k) = (uy9 + force_realy(i,j,k)/2.) / rho9
          uz(i,j,k) = (uz9 + force_realz(i,j,k)/2.) / rho9
          ! print*, istep,i,j,k, ux(i,j,k),uy(i,j,k),uz(i,j,k)

          ! read(*,*)
        
          


      enddo !x
      enddo !y
      enddo !z
      end subroutine macrovar
!=============================================================================
!@subroutine rhoupdat
!@desc Updates only density of the local fluid nodes, currently only used in
!      the pre-relaxation process.
!=============================================================================
      subroutine rhoupdat 
      use var_inc
      implicit none 

      integer ip

      rho = f(0,:,:,:)
      do ip = 1,npop-1
        rho = rho + f(ip,:,:,:)
      end do

      end subroutine rhoupdat 
!===================================================================
!@subroutine avedensity
!@desc Calculates the average density for the entire fluid domain and
!removes
!      it from all densities of local nodes. This is to account for loss
!      of mass in the interpolation bounceback.
!=============================================================================
      subroutine avedensity
      use var_inc
      ! use mpi
      implicit none

      integer iz, iy, ix
      integer nfluid0, nfluidtotal
      real rhomean0, rhomean

      rhomean0 = 0.d0
      nfluid0 = count(ibnodes(1:lx,1:ly,1:lz) < 0)
      rhomean0 = sum(rho,MASK = (ibnodes(1:lx,1:ly,1:lz) < 0))
      ! CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      ! CALL MPI_ALLREDUCE(nfluid0,nfluidtotal,1,MPI_INTEGER,MPI_SUM,mpi_comm_world,ierr)
      ! CALL MPI_ALLREDUCE(rhomean0,rhomean,1,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)

      rhomean = rhomean/dfloat(nfluidtotal)

      do iz = 1,lz
        do iy = 1,ly
          do ix = 1,lx
            rho(ix,iy,iz) = rho(ix,iy,iz) - rhomean
          enddo
        enddo
      enddo

      end subroutine avedensity
!==================================================================
      SUBROUTINE FORCING
      ! use mpi
      use var_inc
      implicit none
      integer ixs,ihh,i,j,k,jj,kk
      real x9,y9,z9

      force_realx(:,:,:) = 0.0
      force_realy(:,:,:) = force_in_y*force_mag
      force_realz(:,:,:) = 0.0

      RETURN
      END SUBROUTINE FORCING
!=================================
      SUBROUTINE FORCINGP
      ! use mpi
      use var_inc
      implicit none
      integer ixs,ihh,i,j,k,jj,kk
      real x9,y9,z9,Amp0,beta9,gamma9,Tpd, phase9,Tpdp,ixs0

      Tpd = 2000.
      Tpdp = 1500.
      beta9 = 3.0
      gamma9=4.0
      ixs0 = 2
      Amp0 = 100.00*beta9/real(ny)*sin(pi2*real(istep)/Tpd)
!     phase9 = sin(pi2*real(istep)/Tpdp)
      phase9 = 0.25

      force_realx(:,:,:) = 0.0
      force_realy(:,:,:) = force_in_y
      force_realz(:,:,:) = 0.0
! add some perturbation
      ixs = ixs0
      ihh = lxh/2
      do k=1,lz
      kk = k + indz*lz
      z9 = pi2*(real(kk)-0.5)/real(nz)
      do j=1,ly
      jj = j + indy*ly
      y9 = pi2*(real(jj)-0.5)/real(ny)
      do i=1,ihh
      x9 = pi2*(real(i)-0.5)/real(ihh)
      force_realx(i+ixs,j,k) = force_in_y*0.5*Amp0*real(ihh)*(1.-cos(x9))  &
                                 *cos(beta9*y9)*cos(gamma9*z9)
      force_realy(i+ixs,j,k) = force_in_y*(1.0-Amp0*real(ny)/beta9    &
                       *sin(x9)*sin(beta9*y9)*cos(gamma9*z9))
      force_realz(i+ixs,j,k) = force_in_y*0.5*Amp0*real(nz)/gamma9*   &
                        sin(x9)*cos(beta9*y9)*sin(gamma9*z9)
      end do
      end do
      end do

      ihh = lxh/2
      ixs = nx - ixs0 - ihh
      do k=1,lz
      kk = k + indz*lz
      z9 = pi2*(real(kk)-0.5)/real(nz)
      do j=1,ly
      jj = j + indy*ly
      y9 = pi2*( (real(jj)-0.5)/real(ny) + phase9 )
      do i=1,ihh
      x9 = pi2*(real(i)-0.5)/real(ihh)
      force_realx(i+ixs,j,k) = -force_in_y*0.5*Amp0*real(ihh)*(1.-cos(x9))  &
                                 *cos(beta9*y9)*cos(gamma9*z9)
      force_realy(i+ixs,j,k) = force_in_y*(1.0+Amp0*real(ny)/beta9    &
                       *sin(x9)*sin(beta9*y9)*cos(gamma9*z9))
      force_realz(i+ixs,j,k) = -force_in_y*0.5*Amp0*real(nz)/gamma9*   &
                        sin(x9)*cos(beta9*y9)*sin(gamma9*z9)
      end do
      end do
      end do

      RETURN
      END SUBROUTINE FORCINGP

