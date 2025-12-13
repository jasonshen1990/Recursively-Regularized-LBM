subroutine Hermite_polynomial
    use var_inc
    implicit none

    integer i, j, k, ip

    H000 = 1.
    do ip = 0, npop-1
        H100(ip) = cix(ip) / cs
        H010(ip) = ciy(ip) / cs
        H001(ip) = ciz(ip) / cs
        H110(ip) = cix(ip) * ciy(ip) / RT
        H101(ip) = cix(ip) * ciz(ip) / RT
        H011(ip) = ciy(ip) * ciz(ip) / RT
        H200(ip) = cix(ip) * cix(ip) / RT - 1
        H020(ip) = ciy(ip) * ciy(ip) / RT - 1
        H002(ip) = ciz(ip) * ciz(ip) / RT - 1
        !third
        H210(ip) = cix(ip) * cix(ip) * ciy(ip) / cs3 - ciy(ip) / cs
        H201(ip) = cix(ip) * cix(ip) * ciz(ip) / cs3 - ciz(ip) / cs
        H120(ip) = ciy(ip) * ciy(ip) * cix(ip) / cs3 - cix(ip) / cs
        H021(ip) = ciy(ip) * ciy(ip) * ciz(ip) / cs3 - ciz(ip) / cs
        H102(ip) = ciz(ip) * ciz(ip) * cix(ip) / cs3 - cix(ip) / cs
        H012(ip) = ciz(ip) * ciz(ip) * ciy(ip) / cs3 - ciy(ip) / cs
        H111(ip) = cix(ip) * ciy(ip) * ciz(ip) / cs3
    end do

end subroutine Hermite_polynomial
!==========================================================================
      real function feqRR(ip, rho9, ux9, uy9, uz9)
      use var_inc
      implicit none
      integer ip
      real rho9, ux9, uy9, uz9
      real uu, vv, ww
      real uv, eu, sec, thir, four, fifsix
      real a110, a101, a011, a200, a020, a002  ! Second
      real a210, a201, a120, a021, a102, a012, a111  ! Third

      uu = ux9 / cs
      vv = uy9 / cs
      ww = uz9 / cs
      a110 = uu * vv; a101 = uu * ww; a011 = vv * ww
      a200 = uu * uu; a020 = vv * vv; a002 = ww * ww  ! Second
      a210 = a200 * vv
      a201 = a200 * ww
      a120 = a020 * uu
      a021 = a020 * ww
      a102 = a002 * uu
      a012 = a002 * vv
      a111 = a110 * ww ! Third

      eu = H100(ip)*uu + H010(ip)*vv + H001(ip)*ww
      sec = 2. * ( H110(ip) * a110 + H101(ip) * a101 + H011(ip) * a011 ) &
          + H200(ip) * a200 + H020(ip) * a020 + H002(ip) * a002
      thir = H210(ip) * a210 + H201(ip) * a201 + H120(ip) * a120 & 
           + H021(ip) * a021 + H102(ip) * a102 + H012(ip) * a012 &
           + 2. * H111(ip) * a111

      feqRR = 1. + eu + 0.5 * ( sec + thir )

      feqRR = tp(ip) * rho9 * feqRR
      return
      end function feqRR
!==================================================================
!==========================================================================
      subroutine comptfneqRR(rho9, ux9, uy9, uz9, f9, Fbar, fneq9)
      use var_inc
     
      implicit none
      integer ip
      real, dimension(0:npop-1) :: f9, fneq9
      real rho9, ux9, uy9, uz9
      real uu, vv, ww
      real uv, eu, sec, thir, four, fifsix
      real a110, a101, a011, a200, a020, a002  ! Second
      real a210, a201, a120, a021, a102, a012, a111  ! Third
      real a220, a202, a022, a211, a121, a112 ! Fourth
      real a221, a212, a122, a222  ! Fifth and sixth
      real ane110, ane101, ane011, ane200, ane020, ane002  ! Second
      real ane210, ane201, ane120, ane021, ane102, ane012, ane111  ! Third
      real ane220, ane202, ane022, ane211, ane121, ane112 ! Fourth
      real ane221, ane212, ane122, ane222  ! Fifth and sixth
      real feqRR,feqRR_4
      real, dimension(0:npop-1) :: Fbar
      real fx9,fy9,fz9,G1,G2,G3
      integer ix, iy, iz

      uu = ux9 / cs
      vv = uy9 / cs
      ww = uz9 / cs
      a110 = rho9 * uu * vv; a101 = rho9 * uu * ww; a011 = rho9 * vv * ww
      a200 = rho9 * uu * uu; a020 = rho9 * vv * vv; a002 = rho9 * ww * ww  ! Second

      ane110 = 0.0; ane101 = 0.0; ane011 = 0.0
      ane200 = 0.0; ane020 = 0.0; ane002 = 0.0
      do ip = 0, npop-1
          ane110 = ane110 + H110(ip) * ( f9(ip) + Fbar(ip)*feqRR_4(ip, rho9, ux9, uy9, uz9) )
          ane101 = ane101 + H101(ip) * ( f9(ip) + Fbar(ip)*feqRR_4(ip, rho9, ux9, uy9, uz9) )
          ane011 = ane011 + H011(ip) * ( f9(ip) + Fbar(ip)*feqRR_4(ip, rho9, ux9, uy9, uz9) )
          ane200 = ane200 + H200(ip) * ( f9(ip) + Fbar(ip)*feqRR_4(ip, rho9, ux9, uy9, uz9) )
          ane020 = ane020 + H020(ip) * ( f9(ip) + Fbar(ip)*feqRR_4(ip, rho9, ux9, uy9, uz9) )
          ane002 = ane002 + H002(ip) * ( f9(ip) + Fbar(ip)*feqRR_4(ip, rho9, ux9, uy9, uz9) )
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
      ane111 = ane110 * ww + ane101 * vv + ane011 * uu  ! Third

    !   eu = H100(ip)*uu + H010(ip)*vv + H001(ip)*ww
    do ip = 0, npop-1
      sec = 2. * ( H110(ip) * ane110 + H101(ip) * ane101 + H011(ip) * ane011 ) &
          + H200(ip) * ane200 + H020(ip) * ane020 + H002(ip) * ane002
      thir = H210(ip) * ane210 + H201(ip) * ane201 + H120(ip) * ane120 &
           + H021(ip) * ane021 + H102(ip) * ane102 + H012(ip) * ane012 &
           + 2. * H111(ip) * ane111
      fneq9(ip) = tp(ip) * ( 0.5 * ( sec + thir ))
    end do
      return
      end subroutine comptfneqRR
!==================================================================
      real function feqRR_4(ip, rho9, ux9, uy9, uz9)
      use var_inc
      implicit none
      integer ip
      real rho9, ux9, uy9, uz9
      real uu, vv, ww
      real uv, eu, sec, thir, four, fifsix
      real a110, a101, a011, a200, a020, a002  ! Second
      real a210, a201, a120, a021, a102, a012, a111  ! Third

      uu = ux9 / cs
      vv = uy9 / cs
      ww = uz9 / cs
      a110 = uu * vv; a101 = uu * ww; a011 = vv * ww
      a200 = uu * uu; a020 = vv * vv; a002 = ww * ww  ! Second
      a210 = a200 * vv
      a201 = a200 * ww
      a120 = a020 * uu
      a021 = a020 * ww
      a102 = a002 * uu
      a012 = a002 * vv
      a111 = a110 * ww ! Third

      eu = H100(ip)*uu + H010(ip)*vv + H001(ip)*ww
      sec = 2. * ( H110(ip) * a110 + H101(ip) * a101 + H011(ip) * a011 ) &
          + H200(ip) * a200 + H020(ip) * a020 + H002(ip) * a002
      thir = H210(ip) * a210 + H201(ip) * a201 + H120(ip) * a120 & 
           + H021(ip) * a021 + H102(ip) * a102 + H012(ip) * a012 &
           + 2. * H111(ip) * a111

      feqRR_4 = 1. + eu + 0.5 * ( sec + thir )

      feqRR_4 = tp(ip) * rho9 * feqRR_4
      return
      end function feqRR_4