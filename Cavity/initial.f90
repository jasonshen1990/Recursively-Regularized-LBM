!==================================================================

      SUBROUTINE initpop 
      use var_inc
      implicit none

      real, dimension(lx,ly,lz):: usqr, G
      integer ip
      real feq, rho9, u9, v9, w9
      integer ix, iy, iz

      usqr = ux*ux + uy*uy + uz*uz
      usqr = 1.5*usqr   

! initialise density
      rho = 1.0 

      DO iz = 1,lz
      DO iy = 1,ly
      DO ix = 1,lx
        rho9 = rho(ix, iy, iz)
        u9   = ux(ix, iy, iz)
        v9   = uy(ix, iy, iz)
        w9   = uz(ix, iy, iz)
!--------set to equilibrium--------
        do ip = 0, npop-1
          f(ip, ix, iy, iz) =  feq(ip, rho9, u9, v9, w9) 
        enddo

     END DO
     END DO
     END DO

      END SUBROUTINE initpop
!==================================================================
