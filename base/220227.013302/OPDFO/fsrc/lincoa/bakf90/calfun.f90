!*==calfun.f90  processed by SPAG 7.50RE at 23:22 on 25 May 2021
      subroutine CALFUN(N, X, F)
      implicit none
!*--CALFUN4
!*** Start of declarations inserted by SPAG
      real * 8 del1, del2, del3, del4, F, FMAx, temp, v12, v13,  &
     &       v14, v23, v24, v34, X, zero
      integer N
!*** End of declarations inserted by SPAG
      common FMAx
      dimension X(*)
      zero = 0.0D0
      F = FMAx
      v12 = X(1) * X(5) - X(4) * X(2)
      v13 = X(1) * X(8) - X(7) * X(2)
      v14 = X(1) * X(11) - X(10) * X(2)
      v23 = X(4) * X(8) - X(7) * X(5)
      v24 = X(4) * X(11) - X(10) * X(5)
      v34 = X(7) * X(11) - X(10) * X(8)
      del1 = v23 * X(12) - v24 * X(9) + v34 * X(6)
      if (del1 > zero) then
          del2 = -v34 * X(3) - v13 * X(12) + v14 * X(9)
          if (del2 > zero) then
              del3 = -v14 * X(6) + v24 * X(3) + v12 * X(12)
              if (del3 > zero) then
                  del4 = -v12 * X(9) + v13 * X(6) - v23 * X(3)
                  if (del4 > zero) then
                      temp = (del1 + del2 + del3 + del4)**3 / (del1 * del2 * del3 * del4)
                      F = DMIN1(temp / 6.0D0, FMAx)
                  end if
              end if
          end if
      end if
      end subroutine CALFUN
