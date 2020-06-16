      subroutine initialize(n, npt, rhobeg, x, xbase, xpt, f, fval,     &
     & xopt, fopt, kopt, bmat, zmat, gq, hq, pq, nf, info, ftarget)

      use pdfomod, only : rp, zero, one, half, is_nan, is_posinf
      implicit none

      integer, intent(in) :: n, npt
      integer, intent(out) :: info, kopt, nf
      real(kind = rp), intent(in) :: rhobeg, ftarget
      real(kind = rp), intent(out) :: xbase(n), xpt(npt, n), f,         &
     & fval(npt), xopt(n), fopt, bmat(npt + n, n), zmat(npt, npt-n-1)
      real(kind = rp), intent(inout) :: x(n)
      real(kind = rp), intent(out) :: gq(n), hq((n*(n + 1))/2), pq(npt)

      integer :: ih, ipt, itemp, jpt
      real(kind = rp) :: fbeg, rhosq, reciq, recip, temp, xipt, xjpt
      logical :: evaluated(npt)

      ! Set the initial elements of XPT, BMAT, HQ, PQ and ZMAT to zero.
      xbase = x
      xpt = zero
      bmat = zero
      zmat = zero
      gq = zero
      hq = zero
      pq = zero

      ! Begin the initialization procedure. The coordinates of the
      ! displacement of the next initial interpolation point from XBASE
      ! are set in XPT(NF, .).
      rhosq = rhobeg*rhobeg
      recip = one/rhosq
      reciq = sqrt(half)/rhosq

      ! At return,
      ! INFO = 0: initialization finishes normally
      ! INFO = 1: return because f <= ftarget
      ! INFO = -1: return because x contains NaN
      ! INFO = -2: return because f is either NaN or positive infinity
      ! INFO = -3: return because the model contains NaN
      info = 0

      ! EVALUATED is a boolean array indicating whether the function
      ! value of the i-th interpolation point has been evaluated.
      ! We need it for a portable counting of the number of function
      ! evaluations, especially if the loop is conducted asynchronously.
      ! However, the loop here is not fully parallelizable if NPT>2N+1,
      ! because the definition XPT(2N+2:end, :) depends on FVAL(1:2N+1).
      evaluated = .false.
      do nf = 1, npt
          ! Set XPT(NF, :)
          if (nf <= 2*n + 1) then
              if (nf >= 2 .and. nf <= n + 1) then
                  xpt(nf, nf - 1) = rhobeg
              else if (nf > n + 1) then
                  xpt(nf, nf - n - 1) = -rhobeg
              end if
          else
              itemp = (nf - n - 2)/n
              jpt = nf - (itemp + 1)*n - 1
              ipt = jpt + itemp
              if (ipt > n) then
                  itemp = jpt
                  jpt = ipt - n
                  ipt = itemp
              end if
              xipt = rhobeg
              if (fval(ipt + n + 1) < fval(ipt + 1)) xipt = -xipt
              xjpt = rhobeg
              if (fval(jpt + n + 1) < fval(jpt + 1)) xjpt = -xjpt
              xpt(nf, ipt) = xipt
              xpt(nf, jpt) = xjpt
          end if

          ! Function evaluation at XPT(NF, :)
          x = xpt(nf, :) + xbase
          if (any(is_nan(x))) then
              f = sum(x)  ! Set F to NaN. It is necessary.
          else
              call calfun(n, x, f)
              evaluated(nf) = .true.
          end if
          fval(nf) = f

          ! Set FOPT, KOPT, and XOPT
          if (nf == 1) then
              fbeg = f
              fopt = f
              kopt = 1
          else if (f < fopt) then
              fopt = f
              kopt = nf
          end if
          xopt = xpt(kopt, :)

          ! Set INFO.
          ! The following should be done after setting fopt and xopt.
          if (any(is_nan(x))) then
              info = -1
              exit
          elseif (f <= ftarget) then
              info = 1
              exit
          elseif (is_posinf(f) .or. is_nan(f)) then
              info = -2
              exit
          end if

          ! Set the coefficients of the initial Lagrange functions
          ! (i.e., bmat and zmat) and the initial quadratic model (i.e.,
          ! gq, hq, pq). This is reached starting from the second
          ! function evaluation, namely when NF > 1.
          if (nf <= 2*n + 1) then
              ! Set the nonzero initial elements of BMAT and the
              ! quadratic model in the cases when NF <= 2*N + 1.
              if (nf >= 2 .and. nf <= n+1) then
                  gq(nf - 1) = (f - fbeg)/rhobeg
                  if (npt < nf + n) then
                      bmat(1, nf - 1) = -one/rhobeg
                      bmat(nf, nf - 1) = one/rhobeg
                      bmat(npt + nf - 1, nf - 1) = -half*rhosq
                  end if
              else if (nf > n + 1) then
                  bmat(nf - n, nf - n - 1) = half/rhobeg
                  bmat(nf, nf - n - 1) = -half/rhobeg
                  zmat(1, nf - n - 1) = -reciq - reciq
                  zmat(nf - n, nf - n - 1) = reciq
                  zmat(nf, nf - n - 1) = reciq
                  ih = ((nf - n - 1)*(nf - n))/2
                  temp = (fbeg - f)/rhobeg
                  hq(ih) = (gq(nf - n - 1) - temp)/rhobeg
                  gq(nf - n - 1) = half*(gq(nf - n - 1) + temp)
              end if
          else
              ! When NF > 2*N+1, set the off-diagonal second derivatives
              ! of the Lagrange functions and the quadratic model.
              ih = (ipt*(ipt - 1))/2 + jpt
              if (xipt < zero) ipt = ipt + n
              if (xjpt < zero) jpt = jpt + n
              zmat(1, nf - n - 1) = recip
              zmat(nf, nf - n - 1) = recip
              zmat(ipt + 1, nf - n - 1) = -recip
              zmat(jpt + 1, nf - n - 1) = -recip
              hq(ih) = (fbeg-fval(ipt+1)-fval(jpt+1)+f)/(xipt*xjpt)
          end if
          ! Check whether NaN occurs in the coefficients. 
          ! Do this only when NF > 1.
          if (nf > 1 .and. (any(is_nan(bmat)) .or. any(is_nan(zmat)).or.
     &     any(is_nan(gq)) .or. any(is_nan(hq)) .or. any(is_nan(pq))))  &
     &     then
              info = -3
              exit
          end if
      end do

      ! If the do loop is conducted seqentially, then the exit value of
      ! NF is UPPER_BOUND + 1.
      ! But this value is NOT portable: it depends on how the do loop
      ! is coducted (seqentially or an parallel) at the backend.
      nf = count(evaluated)

      end subroutine initialize
