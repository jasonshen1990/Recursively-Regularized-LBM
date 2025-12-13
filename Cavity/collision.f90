!===========================================================================
      subroutine collision_RR    
      use var_inc
      implicit none 

      real, dimension(0:npop-1) :: f9, fneqRR
      integer ix, iy, iz, ip 
      real fx9,fy9,fz9,G1,G2,G3
      real*8 rho9, u9, v9, w9, feq, feqRR, feqRR_4
      real, dimension(0:npop-1) :: Fbar
      real fchange

      DO iz = 1,lz
      DO iy = 1,ly
      DO ix = 1,lx

        fx9 = force_realx(ix, iy, iz)
        fy9 = force_realy(ix, iy, iz)
        fz9 = force_realz(ix, iy, iz)

        rho9 = rho(ix, iy, iz)
        u9   = ux(ix, iy, iz)
        v9   = uy(ix, iy, iz)
        w9   = uz(ix, iy, iz)
        f9(:) = f(:, ix, iy, iz)
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
      
          fchange = (1.- 1./tau) * fneqRR(ip) !+ Fbar(ip)*feqRR_4(ip, rho9, u9, v9, w9)
          f(ip, ix, iy, iz) = feqRR(ip,rho9,u9,v9,w9) + fchange

        enddo

     END DO
     END DO
     END DO

      end subroutine collision_RR
!===================================================================
      subroutine macrovar 
      use var_inc
      implicit none 

      integer  iz,iy,ix, ip
      real  rho9, ux9, uy9, uz9
      real, dimension(0:npop-1) :: f9

      do iz = 1,lz
      do iy = 1,ly
      do ix = 1,lx
          f9 = f(:,ix,iy,iz)

          rho9 = f9(0)+f9(1)+f9(2)+f9(3)+f9(4)+f9(5)+f9(6)+f9(7)+f9(8)  &
                +f9(9)+f9(10)+f9(11)+f9(12)+f9(13)+f9(14)+f9(15)+f9(16) &
                +f9(17)+f9(18)+f9(19)+f9(20)+f9(21)+f9(22)+f9(23)+f9(24) &
                +f9(25)+f9(26)

          ux9 = f9(1)-f9(2)+f9(7)+f9(8)-f9(9)-f9(10)+f9(11)+f9(12)-f9(13) &
               -f9(14)+f9(19)+f9(20)+f9(21)+f9(22)-f9(23)-f9(24)-f9(25)-f9(26)


          uy9 = f9(3)-f9(4)+f9(7)-f9(8)+f9(9)-f9(10)+f9(15)+f9(16)-f9(17) &
               -f9(18)+f9(19)+f9(20)-f9(21)-f9(22)+f9(23)+f9(24)-f9(25)-f9(26)

          uz9 = f9(5)-f9(6)+f9(11)-f9(12)+f9(13)-f9(14)+f9(15)-f9(16)+   &
               f9(17)-f9(18)+f9(19)-f9(20)+f9(21)-f9(22)+f9(23)-f9(24) +  &
               f9(25)-f9(26)

          rho(ix,iy,iz) = rho9
          ux(ix,iy,iz) = (ux9 + force_realx(ix,iy,iz)/2.) / rho9
          uy(ix,iy,iz) = (uy9 + force_realy(ix,iy,iz)/2.) / rho9
          uz(ix,iy,iz) = (uz9 + force_realz(ix,iy,iz)/2.) / rho9

      enddo
      enddo
      enddo
     
      end subroutine macrovar
! ========================================================================
      SUBROUTINE FORCING
!      use mpi
      use var_inc
      implicit none
     
      force_realx(:,:,:) = 0.0
      force_realy(:,:,:) = force_in_y
      force_realz(:,:,:) = 0.0

      RETURN
      END SUBROUTINE FORCING
!==========================================================================
      real function feq(ip, rho9, u9, v9, w9)
      use var_inc
      implicit none
      integer ip
      real rho9, u9, v9, w9
      real eu, uv

      eu = ( cix(ip)*u9 + ciy(ip)*v9 + ciz(ip)*w9 ) / RT
      uv = ( u9*u9 + v9*v9 +w9*w9 ) / RT
      feq = tp(ip) * rho9 * ( 1. + eu + 0.5*(eu*eu-uv) )
      return
      end function feq
!==================================================================