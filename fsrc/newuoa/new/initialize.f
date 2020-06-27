      subroutine initialize(n, npt, rhobeg, x, xbase, xpt, f, fval,     &
     & xopt, fopt, kopt, bmat, zmat, gq, hq, pq, nf, info, ftarget)

      use consts, only : rp, zero, one, half
      use infnan
      implicit none

      integer, intent(in) :: n, npt
      integer, intent(out) :: info, kopt, nf
      real(kind = rp), intent(in) :: rhobeg, x(n), ftarget
      real(kind = rp), intent(out) :: xbase(n), xpt(n, npt), f,         &
     & fval(npt), xopt(n), fopt, bmat(npt + n, n), zmat(npt, npt-n-1)
      real(kind = rp), intent(out) :: gq(n), hq((n*(n + 1))/2), pq(npt)

      integer :: ih, ip, ipt(npt), itemp, jp, jpt(npt), npt_fun, npt_mod
      real(kind = rp) :: fbeg, rhosq, reciq, recip, temp, xip, xjp,     &
     & xtemp(n)
      logical :: evaluated(npt)


      ! Set the initial elements of XPT, BMAT, HQ, PQ and ZMAT to zero.
      xbase = x
      xpt = zero
      bmat = zero
      zmat = zero
      gq = zero
      hq = zero
      pq = zero  ! We will not update PQ. It is ZERO at return.

      ! At return,
      ! INFO = 0: initialization finishes normally
      ! INFO = 1: return because f <= ftarget
      ! INFO = -1: return because x contains NaN
      ! INFO = -2: return because f is either NaN or positive infinity
      ! INFO = -3: return because the model contains NaN
      info = 0

      ! Begin the initialization procedure. The coordinates of the
      ! displacement of the next initial interpolation point from XBASE
      ! are set in XPT(NF, .).
      rhosq = rhobeg*rhobeg
      recip = one/rhosq
      reciq = sqrt(half)/rhosq

      ! EVALUATED is a boolean array indicating whether the function
      ! value of the i-th interpolation point has been evaluated.
      ! We need it for a portable counting of the number of function
      ! evaluations, especially if the loop is conducted asynchronously.
      ! However, the loop here is not fully parallelizable if NPT>2N+1,
      ! because the definition XPT(2N+2:end, :) depends on FVAL(1:2N+1).
      evaluated = .false.

      ! NPT_FUN and NPT_MOD equal NPT, unless it turns out necessary to
      ! return due to abnormality (NaN or Inf occurs, or F < FTARGET).
      npt_fun = npt
      npt_mod = npt
      

      ! Evaluate F at the starting point X.
      if (any(is_nan(x))) then  ! X contains NaN. Return immediately.
          f = sum(x)  ! Set F to NaN. It is necessary.
          info = -1
          npt_fun = 0 
          npt_mod = 0 
      else
          call calfun(n, x, f)
          evaluated(1) = .true.
          fval(1) = f
          fbeg = f
          ! Check whether to exit.
          if (f <= ftarget) then
              info = 1
              npt_fun = 0 
              npt_mod = 0 
          end if
          if (is_posinf(f) .or. is_nan(f)) then
              info = -2
              npt_fun = 0 
              npt_mod = 0 
          end if
      end if
      
      ! Set XPT, FVAL, KOPT, FOPT, and XOPT.
      do nf = 2, npt_fun
          ! Set XPT(NF, :)
          if (nf <= 2*n + 1) then
              if (nf <= n + 1) then
                  xpt(nf - 1, nf) = rhobeg
              else
                  xpt(nf - n - 1, nf) = -rhobeg
              end if
          else
              itemp = (nf - n - 2)/n
              jp = nf - (itemp + 1)*n - 1
              ip = jp + itemp
              if (ip > n) then
                  itemp = jp
                  jp = ip - n
                  ip = itemp
              end if
              ipt(nf) = ip
              jpt(nf) = jp
              ! SIGN(A, B) = |A| if B > 0 and -|A| if B < 0. Note that:
              ! 1. When B = 0, the result depends on the compiler, but 
              !    it is normally |A|. 
              ! 2. SIGN(A, B) is defined in terms of |A| but not
              !    directly in terms of A; use with caution when A is
              !    possibly negative.
              xip = sign(rhobeg, (fval(ip + n + 1) - fval(ip + 1)))
              xjp = sign(rhobeg, (fval(jp + n + 1) - fval(jp + 1)))
              xpt(ip, nf) = xip
              xpt(jp, nf) = xjp
          end if

          ! Function evaluation at XPT(NF, :)
          xtemp = xpt(:, nf) + xbase
          if (any(is_nan(xtemp))) then
              f = sum(xtemp)  ! Set F to NaN. It is necessary.
              info = -1
              npt_mod = 0 
              exit
          end if
          call calfun(n, xtemp, f)
          evaluated(nf) = .true.
          fval(nf) = f

          ! Check whether to exit.
          if (f <= ftarget) then
              info = 1
              npt_mod = 0 
              exit
          end if
          if (is_posinf(f) .or. is_nan(f)) then
              info = -2
              npt_mod = 0 
              exit
          end if
      end do

      ! Set GQ, HQ, BMAT, and ZMAT.
      do nf = 2, npt_mod
          ! Set the coefficients of the initial Lagrange functions
          ! (i.e., bmat and zmat) and the initial quadratic model (i.e.,
          ! gq, hq, pq). This is reached starting from the second
          ! function evaluation, namely when NF > 1.
          f = fval(nf)
          if (nf <= 2*n + 1) then
              ! Set the nonzero initial elements of BMAT and the
              ! quadratic model in the cases when NF <= 2*N + 1.
              if (nf <= n + 1) then
                  gq(nf - 1) = (f - fbeg)/rhobeg
                  if (npt < nf + n) then
                      bmat(1, nf - 1) = -one/rhobeg
                      bmat(nf, nf - 1) = one/rhobeg
                      bmat(npt + nf - 1, nf - 1) = -half*rhosq
                  end if
              else 
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

              ! IP, JP, XIP, and XJP will be used below.
              ip = ipt(nf)
              jp = jpt(nf)
              xip = xpt(ip, nf)
              xjp = xpt(jp, nf)

              ih = (ip*(ip - 1))/2 + jp
              if (xip < zero) then 
                  ip = ip + n
              end if
              if (xjp < zero) then
                  jp = jp + n
              end if
              zmat(1, nf - n - 1) = recip
              zmat(nf, nf - n - 1) = recip
              zmat(ip + 1, nf - n - 1) = -recip
              zmat(jp + 1, nf - n - 1) = -recip
              hq(ih) = (fbeg - fval(ip+1) - fval(jp+1) + f)/(xip*xjp)
          end if

          ! Check whether NaN occurs in the coefficients. 
          ! Note that PQ is always ZERO --- the initial HESSIAN only has
          ! the explicit part.
          if ((any(is_nan(bmat)) .or. any(is_nan(zmat)).or.             &
     &     any(is_nan(gq)) .or. any(is_nan(hq)))) then
              info = -3
              exit
          end if
      end do

      ! If the do loop is conducted seqentially, then the exit value of
      ! NF is UPPER_BOUND + 1.
      ! But this value is NOT portable: it depends on how the do loop
      ! is coducted (seqentially or an parallel) at the backend.
      nf = count(evaluated)

      if (nf == 0) then
      ! In this case, the starting X contains NaN. F was set to NaN.
          kopt = 1
          fopt = f
          xopt = zero
      else
          kopt = minloc(fval, dim = 1, mask = evaluated)
          fopt = fval(kopt)
          xopt = xpt(:, kopt)
      end if

      return

      end subroutine initialize
