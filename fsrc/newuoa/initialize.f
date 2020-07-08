      subroutine initialize(n, npt, rhobeg, x, xbase, xpt, f, fval,     &
     & xopt, fopt, kopt, bmat, zmat, gq, hq, pq, nf, info, ftarget)

      use consts, only : rp, zero, one, half
      use infnan
      implicit none

      integer, intent(in) :: n, npt
      integer, intent(out) :: info, kopt, nf
      real(RP), intent(in) :: rhobeg, x(n), ftarget
      real(RP), intent(out) :: xbase(n), xpt(n, npt), f, fval(npt),     &
     & xopt(n), fopt, bmat(n, npt + n), zmat(npt, npt-n-1)
      real(RP), intent(out) :: gq(n), hq(n, n), pq(npt)

      integer :: k, ih, ip, ipt(npt), itemp, jp, jpt(npt), npt1, npt2
      real(RP) :: fbeg, fip, fjp, rhosq, reciq, recip, xip, xjp,xtemp(n)
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
      ! are set in XPT(:, .).
      rhosq = rhobeg*rhobeg
      recip = one/rhosq
      reciq = sqrt(half)/rhosq

      ! EVALUATED is a boolean array indicating whether the function
      ! value of the i-th interpolation point has been evaluated.
      ! We need it for a portable counting of the number of function
      ! evaluations, especially if the loop is conducted asynchronously.
      ! However, the loop here is not fully parallelizable if NPT>2N+1,
      ! because the definition XPT(;, 2N+2:end) depends on FVAL(1:2N+1).
      evaluated = .false.

      ! NPT1 and NPT2 equal NPT, unless it turns out necessary to
      ! return due to abnormality (NaN or Inf occurs, or F < FTARGET).
      npt1 = npt
      npt2 = npt
      
      ! When K > 2*N + 1, the IPT(K) and JPT(K) entries of XPT(:, K)
      ! will be RHOBEG or -RHOBEG, depending on the values of F(IPT(K))
      ! and F(JPT(K)). Consequently, the Hessian of the quadratic model
      ! will get a possibly nonzero (IPT(K), JPT(K)) entry.
      ! We initilize IPT and JPT to 1 in case the initialization of XPT
      ! aborts before finishing, which will leave IPT and JPT undefined.
      ipt = 1
      jpt = 1

      ! Evaluate F at the starting point X.
      if (any(is_nan(x))) then  ! X contains NaN. Return immediately.
          f = sum(x)  ! Set F to NaN. It is necessary.
          info = -1
          npt1 = 0 
          npt2 = 0
      else
          call calfun(n, x, f)
          evaluated(1) = .true.
          fval(1) = f
          fbeg = f
          ! Check whether to exit.
          if (f <= ftarget) then
              info = 1
              npt1 = 0 
              npt2 = 0
          end if
          if (is_posinf(f) .or. is_nan(f)) then
              info = -2
              npt1 = 0 
              npt2 = 0
          end if
      end if
      
      ! Set XPT, FVAL, KOPT, FOPT, and XOPT.
      
      ! Set XPT(:, 2 : 2*N + 1). 
      do k = 2, min(npt1, n + 1)
          xpt(k - 1, k) = rhobeg
      end do
      do k = n+2, min(npt1, 2*n + 1)
          xpt(k - n - 1, k) = -rhobeg
      end do
       
      ! Set FVAL(2 : 2*N + 1) by evaluating F. Totally parallelizable.
      do k = 2, min(npt1, 2*n + 1)
          xtemp = xpt(:, k) + xbase
          if (any(is_nan(xtemp))) then
              f = sum(xtemp)  ! Set F to NaN. It is necessary.
              info = -1
              npt2 = 0 
              exit
          end if
          call calfun(n, xtemp, f)
          evaluated(k) = .true.
          fval(k) = f

          ! Check whether to exit.
          if (f <= ftarget) then
              info = 1
              npt2 = 0 
              exit
          end if
          if (is_posinf(f) .or. is_nan(f)) then
              info = -2
              npt2 = 0 
              exit
          end if
      end do

      ! Set XPT(:, 2*N + 2 : NPT). It depends on FVAL(2 : 2*N + 1).
      do k = 2*n + 2, npt2
          itemp = (k - n - 2)/n
          jp = k - (itemp + 1)*n - 1
          ip = jp + itemp
          if (ip > n) then
              itemp = jp
              jp = ip - n
              ip = itemp
          end if
          ipt(k) = ip
          jpt(k) = jp
          ! XIP and XJP are the only places that depend on FVAL(2:2*N+1)
          ! SIGN(A, B) = |A| if B > 0 and -|A| if B < 0. Note that:
          ! 1. When B = 0, the result depends on the compiler, but 
          !    it is normally |A|. 
          ! 2. SIGN(A, B) is defined in terms of |A| but not
          !    directly in terms of A; use with caution when A is
          !    possibly negative.
          ! 3. MATLAB has also a SIGN function, but the definition is
          !    SIGN(A) = 1 if A > 0, 0 if A = 0, and -1 if A < 0.
          xpt(ip, k) = sign(rhobeg, fval(ip + n + 1) - fval(ip + 1)) 
          xpt(jp, k) = sign(rhobeg, fval(jp + n + 1) - fval(jp + 1))
      end do

      ! Set FVAL(2*N + 2 : NPT) by evaluating F. Totally parallelizable.
      do k = 2*n + 2, npt2
          xtemp = xpt(:, k) + xbase
          if (any(is_nan(xtemp))) then
              f = sum(xtemp)  ! Set F to NaN. It is necessary.
              info = -1
              exit
          end if
          call calfun(n, xtemp, f)
          evaluated(k) = .true.
          fval(k) = f

          ! Check whether to exit.
          if (f <= ftarget) then
              info = 1
              exit
          end if
          if (is_posinf(f) .or. is_nan(f)) then
              info = -2
              exit
          end if
      end do

      ! Set NF, KOPT, FOPT, and XOPT.
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

      ! We return immediately if NF < NPT, because it indicates that 
      ! something abnormality ocurred during the initialization of XPT
      ! and FVAL, causing the initialization to abort.  
      if (nf < npt) then
          return
      end if

      ! Set GQ and HQ.

      ! Set GQ by forward difference.
      gq(1 : n) = (fval(2 : n + 1) - fbeg)/rhobeg
      ! If possible, revise GQ to central difference. 
      k = min(npt - n - 1, n)
      gq(1 : k) = half*(gq(1 : k) + (fbeg - fval(n+2 : n+1+k))/rhobeg)

      ! Set the diagonal of HQ  by 2nd-order central finite difference.
      do k = 1, min(npt - n - 1, n) 
          hq(k, k) = ((fval(k + 1) - fbeg)/rhobeg -                     &
     &     (fbeg - fval(k + n + 1))/rhobeg)/rhobeg
      end do
      ! When NPT > 2*N + 1, set the off-diagonal entries of HQ.
      do k = 2*n + 2, npt 
          ! IP, JP, XIP, and XJP will be used below.
          ip = ipt(k)
          jp = jpt(k)
          xip = xpt(ip, k)
          xjp = xpt(jp, k)
          ih = (ip*(ip - 1))/2 + jp
          if (xip < zero) then 
              fip = fval(ip + n + 1)
          else
              fip = fval(ip + 1)
          end if
          if (xjp < zero) then
              fjp = fval(jp + n + 1)
          else
              fjp = fval(jp + 1)
          end if
          hq(ip, jp) = (fbeg - fip - fjp + fval(k))/(xip*xjp)
          hq(jp, ip) = hq(ip, jp)
      end do

      ! Set BMAT and ZMAT. They depend on XPT but not FVAL.

      ! Set the nonzero initial elements of BMAT. 
      ! When NPT >= 2*N + 1, this defines BMAT completely; 
      ! When NPT <= 2*N, this defines the first NPT-N-1 rows of BMAT.
      do k = 1, min(npt - n - 1, n)
          bmat(k, k + 1) = half/rhobeg
          bmat(k, n + k + 1) = -half/rhobeg
      end do

      ! When NPT <= 2*N, set the NPT - N to N rows of BMAT. 
      do k = npt - n, n 
          bmat(k, 1) = -one/rhobeg
          bmat(k, k + 1) = one/rhobeg
          bmat(k, npt + k) = -half*rhosq
      end do
                
      ! Set the nonzero initial elements of ZMAT. 
      ! When NPT <= 2*N + 1, this defines ZMAT completely; 
      ! When NPT > 2*N + 1, this defines the first N columns of ZMAT.
      do k = 1, min(npt - n - 1, n) 
          zmat(1, k) = - reciq - reciq
          zmat(k + 1, k) = reciq
          zmat(k + n + 1, k) = reciq
      end do

      ! When NPT > 2*N+1, set the N + 1 to NPT - N - 1 columns of ZMAT.
      do k = n + 1, npt - n - 1
          ! IP, JP, XIP, and XJP will be used below.
          ip = ipt(k + n + 1)
          jp = jpt(k + n + 1)
          xip = xpt(ip, k + n + 1)
          xjp = xpt(jp, k + n + 1)
          if (xip < zero) then 
              ip = ip + n
          end if
          if (xjp < zero) then
              jp = jp + n
          end if
          zmat(1, k) = recip
          zmat(k + n + 1, k) = recip
          zmat(ip + 1, k) = -recip
          zmat(jp + 1, k) = -recip
      end do

!      ! Check whether NaN occurs in the coefficients. 
!      if ((any(is_nan(bmat)) .or. any(is_nan(zmat)).or.                 &
!     &     any(is_nan(gq)) .or. any(is_nan(hq)))) then
!          info = -3
!      end if

      return

      end subroutine initialize
