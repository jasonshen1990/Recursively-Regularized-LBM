!==================================================================
      subroutine streaming
      use var_inc
      implicit none
 
      real, dimension(0:npop-1)  :: f9, fneqRR
      real, dimension(0:npop-1,lx,lz)  :: tmpxL, tmpxR
      real, dimension(0:npop-1,ly,lz)  :: tmpxB, tmpxT
      real, dimension(0:npop-1,lx,ly)  :: tmpxF, tmpxH
      real rho9, u9, v9, w9
      real*8 feq, feqRR,feqRR_4
      integer ip, ix, iy, iz
      real fx9,fy9,fz9,G1,G2,G3
      real, dimension(0:npop-1) :: Fbar
! save 6 layers for bounce back after streaming(Each two for three directions)
      tmpxB = f(:, 1,  :, :)  
      tmpxT = f(:, lx, :, :)
      tmpxL = f(:, :,  1, :)
      tmpxR = f(:, :, ly, :)

! y+ move Dir 3,7,8,15,17  

      f(3,:,:,:) = cshift(f(3,:,:,:), shift = -1, dim = 2)
      f(7,:,:,:) = cshift(f(7,:,:,:), shift = -1, dim = 2)
      f(9,:,:,:) = cshift(f(9,:,:,:), shift = -1, dim = 2)
      f(15,:,:,:) = cshift(f(15,:,:,:), shift = -1, dim = 2)
      f(16,:,:,:) = cshift(f(16,:,:,:), shift = -1, dim = 2)
      f(19,:,:,:) = cshift(f(19,:,:,:), shift = -1, dim = 2)
      f(20,:,:,:) = cshift(f(20,:,:,:), shift = -1, dim = 2)
      f(23,:,:,:) = cshift(f(23,:,:,:), shift = -1, dim = 2)
      f(24,:,:,:) = cshift(f(24,:,:,:), shift = -1, dim = 2)

! y- move Dir 4,9,10,16,18  
      f(4,:,:,:) = cshift(f(4,:,:,:), shift = 1, dim = 2)
      f(8,:,:,:) = cshift(f(8,:,:,:), shift = 1, dim = 2)
      f(10,:,:,:) = cshift(f(10,:,:,:), shift = 1, dim = 2)
      f(17,:,:,:) = cshift(f(17,:,:,:), shift = 1, dim = 2)
      f(18,:,:,:) = cshift(f(18,:,:,:), shift = 1, dim = 2)
      f(21,:,:,:) = cshift(f(21,:,:,:), shift = 1, dim = 2)
      f(22,:,:,:) = cshift(f(22,:,:,:), shift = 1, dim = 2)
      f(25,:,:,:) = cshift(f(25,:,:,:), shift = 1, dim = 2)
      f(26,:,:,:) = cshift(f(26,:,:,:), shift = 1, dim = 2)

! z+ move Dir 5,11,12,15,16   
      f(5,:,:,:) = cshift(f(5,:,:,:), shift = -1, dim = 3)
      f(11,:,:,:) = cshift(f(11,:,:,:), shift = -1, dim = 3)
      f(13,:,:,:) = cshift(f(13,:,:,:), shift = -1, dim = 3)
      f(15,:,:,:) = cshift(f(15,:,:,:), shift = -1, dim = 3)
      f(17,:,:,:) = cshift(f(17,:,:,:), shift = -1, dim = 3)
      f(19,:,:,:) = cshift(f(19,:,:,:), shift = -1, dim = 3)
      f(21,:,:,:) = cshift(f(21,:,:,:), shift = -1, dim = 3)
      f(23,:,:,:) = cshift(f(23,:,:,:), shift = -1, dim = 3)
      f(25,:,:,:) = cshift(f(25,:,:,:), shift = -1, dim = 3)

! z- move Dir: 6,13,14,17,18   
      f(6,:,:,:) = cshift(f(6,:,:,:), shift = 1, dim = 3)
      f(12,:,:,:) = cshift(f(12,:,:,:), shift = 1, dim = 3)
      f(14,:,:,:) = cshift(f(14,:,:,:), shift = 1, dim = 3)
      f(16,:,:,:) = cshift(f(16,:,:,:), shift = 1, dim = 3)
      f(18,:,:,:) = cshift(f(18,:,:,:), shift = 1, dim = 3)
      f(20,:,:,:) = cshift(f(20,:,:,:), shift = 1, dim = 3)
      f(22,:,:,:) = cshift(f(22,:,:,:), shift = 1, dim = 3)
      f(24,:,:,:) = cshift(f(24,:,:,:), shift = 1, dim = 3)
      f(26,:,:,:) = cshift(f(26,:,:,:), shift = 1, dim = 3)

! x+ move Dir 1,7,9,11,13 
      f(1,:,:,:) = cshift(f(1,:,:,:), shift = -1, dim = 1)
      f(7,:,:,:) = cshift(f(7,:,:,:), shift = -1, dim = 1)
      f(8,:,:,:) = cshift(f(8,:,:,:), shift = -1, dim = 1)
      f(11,:,:,:) = cshift(f(11,:,:,:), shift = -1, dim = 1)
      f(12,:,:,:) = cshift(f(12,:,:,:), shift = -1, dim = 1)
      f(19,:,:,:) = cshift(f(19,:,:,:), shift = -1, dim = 1)
      f(20,:,:,:) = cshift(f(20,:,:,:), shift = -1, dim = 1)
      f(21,:,:,:) = cshift(f(21,:,:,:), shift = -1, dim = 1)
      f(22,:,:,:) = cshift(f(22,:,:,:), shift = -1, dim = 1)
    
! x- move Dir 2,8,10,12,14
      f(2,:,:,:) = cshift(f(2,:,:,:), shift = 1, dim = 1)
      f(9,:,:,:) = cshift(f(9,:,:,:), shift = 1, dim = 1)
      f(10,:,:,:) = cshift(f(10,:,:,:), shift = 1, dim = 1)
      f(13,:,:,:) = cshift(f(13,:,:,:), shift = 1, dim = 1)
      f(14,:,:,:) = cshift(f(14,:,:,:), shift = 1, dim = 1)
      f(23,:,:,:) = cshift(f(23,:,:,:), shift = 1, dim = 1)
      f(24,:,:,:) = cshift(f(24,:,:,:), shift = 1, dim = 1)
      f(25,:,:,:) = cshift(f(25,:,:,:), shift = 1, dim = 1)
      f(26,:,:,:) = cshift(f(26,:,:,:), shift = 1, dim = 1)


! Now, midlink bounce back
! bottom
      f(1,  1, :, :)  = tmpxB(2,  :, :)
      f(7,  1, :, :)  = tmpxB(10, :, :)
      f(8,  1, :, :)  = tmpxB(9,  :, :)
      f(11, 1, :, :)  = tmpxB(14, :, :)
      f(12, 1, :, :)  = tmpxB(13, :, :)
      f(19, 1, :, :)  = tmpxB(26, :, :)
      f(20, 1, :, :)  = tmpxB(25, :, :)
      f(21, 1, :, :)  = tmpxB(24, :, :)
      f(22, 1, :, :)  = tmpxB(23, :, :)
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

          fx9 = force_realx(ix, iy, iz)
          fy9 = force_realy(ix, iy, iz)
          fz9 = force_realz(ix, iy, iz)
          G3 = u9*fx9 + v9*fy9 + w9*fz9
          Fbar(0) = -0.5*G3/RT
  
          do ip=1,npop-1
            G1 = cix(ip)*fx9 + ciy(ip)*fy9 + ciz(ip)*fz9
            Fbar(ip) = 0.5*(G1 - G3)/RT
          enddo

          fneqRR = 0.0
          call comptfneqRR(rho9, u9, v9, w9, f9,Fbar, fneqRR)
!--------Perform collision--------
          do ip = 0, npop-1
            if (cix(ip).gt.0) then
            f(ip, ix, iy, iz) = feqRR(ip,rho9,u9,v9,w9) + fneqRR(ip) - Fbar(ip)*feqRR_4(ip,rho9,u9,v9,w9)
            end if
          enddo
      END DO
      END DO
      END DO
!top
      f(2,  lx, :, :)  = tmpxT(1,  :, :)
      f(9,  lx, :, :)  = tmpxT(8,  :, :) + 6.*tp(9)*u0
      f(10, lx, :, :)  = tmpxT(7,  :, :) - 6.*tp(10)*u0
      !  write(*,*) 'monitor', u0, tp(10), -6.*tp(10)*u0
      !  pause
      f(13, lx, :, :)  = tmpxT(12, :, :)
      f(14, lx, :, :)  = tmpxT(11, :, :)
      f(23, lx, :, :)  = tmpxT(22, :, :) + 6.*tp(23)*u0
      f(24, lx, :, :)  = tmpxT(21, :, :) + 6.*tp(24)*u0
      f(25, lx, :, :)  = tmpxT(20, :, :) - 6.*tp(25)*u0
      f(26, lx, :, :)  = tmpxT(19, :, :) - 6.*tp(26)*u0
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

          fx9 = force_realx(ix, iy, iz)
          fy9 = force_realy(ix, iy, iz)
          fz9 = force_realz(ix, iy, iz)
          G3 = u9*fx9 + v9*fy9 + w9*fz9
          Fbar(0) = -0.5*G3/RT
  
          do ip=1,npop-1
            G1 = cix(ip)*fx9 + ciy(ip)*fy9 + ciz(ip)*fz9
            Fbar(ip) = 0.5*(G1 - G3)/RT
          enddo
          fneqRR = 0.0
          call comptfneqRR(rho9, u9, v9, w9, f9,Fbar, fneqRR)
!--------Perform collision--------
          do ip = 0, npop-1
            if (cix(ip).lt.0) then
            f(ip, ix, iy, iz) = feqRR(ip,rho9,u9,v9,w9) + fneqRR(ip)- Fbar(ip)*feqRR_4(ip,rho9,u9,v9,w9)
            end if
          enddo
      END DO
      END DO
      END DO
! left
      f(3,  :, 1, :)  = tmpxL(4,  :, :)
      f(15, :, 1, :)  = tmpxL(18, :, :)
      f(16, :, 1, :)  = tmpxL(17, :, :)
      f(20, :, 1, :)  = tmpxL(25, :, :)
      f(7,  :, 1, :)  = tmpxL(10, :, :)
      f(19, :, 1, :)  = tmpxL(26, :, :)
      f(24, :, 1, :)  = tmpxL(21, :, :)
      f(9 , :, 1, :)  = tmpxL(8,  :, :)
      f(23, :, 1, :)  = tmpxL(22, :, :)
      DO iz = 1,lz
      DO iy = 1,1
      DO ix = 1,lx
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

          fx9 = force_realx(ix, iy, iz)
          fy9 = force_realy(ix, iy, iz)
          fz9 = force_realz(ix, iy, iz)
          G3 = u9*fx9 + v9*fy9 + w9*fz9
          Fbar(0) = -0.5*G3/RT
  
          do ip=1,npop-1
            G1 = cix(ip)*fx9 + ciy(ip)*fy9 + ciz(ip)*fz9
            Fbar(ip) = 0.5*(G1 - G3)/RT
          enddo
          fneqRR = 0.0
          call comptfneqRR(rho9, u9, v9, w9, f9,Fbar, fneqRR)
!--------Perform collision--------
          do ip = 0, npop-1
            if (ciy(ip).gt.0) then
            f(ip, ix, iy, iz) = feqRR(ip,rho9,u9,v9,w9) + fneqRR(ip)- Fbar(ip)*feqRR_4(ip,rho9,u9,v9,w9)
            end if
          enddo
      END DO
      END DO
      END DO
! right
      f(4,  :, ly, :)  = tmpxR(3,  :, :)
      f(18, :, ly, :)  = tmpxR(15, :, :)
      f(17, :, ly, :)  = tmpxR(16, :, :)
      f(22, :, ly, :)  = tmpxR(23, :, :)
      f(8 , :, ly, :)  = tmpxR(9,  :, :)
      f(21, :, ly, :)  = tmpxR(24, :, :)
      f(26, :, ly, :)  = tmpxR(19, :, :)
      f(10, :, ly, :)  = tmpxR(7 , :, :)
      f(25, :, ly, :)  = tmpxR(20, :, :)
      DO iz = 1,lz
      DO iy = ly,ly
      DO ix = 1,lx
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

          fx9 = force_realx(ix, iy, iz)
          fy9 = force_realy(ix, iy, iz)
          fz9 = force_realz(ix, iy, iz)
          G3 = u9*fx9 + v9*fy9 + w9*fz9
          Fbar(0) = -0.5*G3/RT
  
          do ip=1,npop-1
            G1 = cix(ip)*fx9 + ciy(ip)*fy9 + ciz(ip)*fz9
            Fbar(ip) = 0.5*(G1 - G3)/RT
          enddo
          fneqRR = 0.0
          call comptfneqRR(rho9, u9, v9, w9, f9,Fbar, fneqRR)
!--------Perform collision--------
          do ip = 0, npop-1
            if (ciy(ip).lt.0) then
            f(ip, ix, iy, iz) = feqRR(ip,rho9,u9,v9,w9) + fneqRR(ip)- Fbar(ip)*feqRR_4(ip,rho9,u9,v9,w9)
            end if
          enddo
      END DO
      END DO
      END DO

      end subroutine streaming
!==================================================================




      
