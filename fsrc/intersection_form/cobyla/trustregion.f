!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of trustregion.f90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 10-Aug-2021.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      module trustregion_mod

      contains

      subroutine trstlp(n, m, A, b, rho, d, ifull)

! Generic modules
      use consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, TENTH, EPS, H&
     &UGENUM, DEBUGGING, SRNLEN
      use info_mod, only : FTARGET_ACHIEVED, MAXFUN_REACHED, TRSUBP_FAIL&
     &ED, SMALL_TR_RADIUS, NAN_X, NAN_INF_F
      use infnan_mod, only : is_nan, is_posinf, is_finite
      use debug_mod, only : errstop
      use output_mod, only : retmssg, rhomssg, fmssg
      use lina_mod, only : inprod, matprod, eye, planerot, isminor

      implicit none

      integer(IK), intent(in) :: n
      integer(IK), intent(in) :: m
      real(RP), intent(in) :: A(:, :)
      ! (N, M+1)
      real(RP), intent(in) :: b(:)
      ! M+1
      real(RP), intent(in) :: rho
      real(RP), intent(inout) :: d(:)
      ! N
      integer(IK), intent(out) :: ifull

      real(RP) :: hypt
      real(RP) :: z(n, n)
      real(RP) :: zdota(n)
      real(RP) :: sdirn(n)
      real(RP) :: vmultc(size(A, 2))
      real(RP) :: vmultd(size(A, 2))
      real(RP) :: cgrad(n)
      real(RP) :: cgz(n)
      real(RP) :: cgzabs(n)
      real(RP) :: cgzk
      real(RP) :: cgzkabs
      real(RP) :: dnew(n)

      real(RP) :: alpha
      real(RP) :: beta
      real(RP) :: dd
      real(RP) :: grot(2, 2)
      real(RP) :: optnew
      real(RP) :: optold
      real(RP) :: ratio
      real(RP) :: cstrv
      real(RP) :: cvold
      real(RP) :: sd
      real(RP) :: sp
      real(RP) :: spabs
      real(RP) :: ss
      real(RP) :: step
      real(RP) :: stpful
      real(RP) :: summ
      real(RP) :: summabs
      real(RP) :: summd
      real(RP) :: temp
      real(RP) :: tempa
      real(RP) :: tmpv(size(A, 2))
      real(RP) :: tmpvabs(size(A, 2))
      real(RP) :: vsave
      real(RP) :: zdotw
      real(RP) :: zdvabs
      real(RP) :: zdwabs
      real(RP) :: dold(size(d))
      ! N
      integer(IK) :: i
      integer(IK) :: iact(size(A, 2))
      integer(IK) :: icon
      integer(IK) :: icount
      integer(IK) :: isave
      integer(IK) :: iter
      integer(IK) :: j
      integer(IK) :: k
      integer(IK) :: kk
      integer(IK) :: kl
      integer(IK) :: kw
      integer(IK) :: maxiter
      integer(IK) :: mcon
      integer(IK) :: nact
      integer(IK) :: nactx
      integer(IK) :: stage


! This subroutine calculates an N-component vector D by the following two stages. In the first
! stage, D is set to the shortest vector that minimizes the greatest violation of the constraints
!       dot_product(A(1:N, K), D) >= B(K),  K= 1, 2, 3, ..., M,
! subject to the Euclidean length of D being at most RHO. If its length is strictly less than RHO,
! then we use the resultant freedom in D to minimize the objective function
!       dot_product(-A(1:N, M+1), D)
! subject to no increase in any greatest constraint violation. This notation allows the gradient of
! the objective function to be regarded as the gradient of a constraint. Therefore the two stages
! are distinguished by MCON == M and MCON > M respectively.
!
! It is possible that a degeneracy may prevent D from attaining the target length RHO. Then the
! value IFULL=0 would be set, but usually IFULL=1 on return.
!
! In general NACT is the number of constraints in the active set and IACT(1),...,IACT(NACT) are
! their indices, while the remainder of IACT contains a permutation of the remaining constraint
! indices.
!
! Further, Z is an orthogonal matrix whose first NACT columns can be regarded as the result of
! Gram-Schmidt applied to the active constraint gradients. For J=1, 2, ..., NACT, the number
! ZDOTA(J) is the scalar product of the J-th column of Z with the gradient of the J-th active
! constraint. D is the current vector of variables and here the residuals of the active constraints
! should be zero. Further, the active constraints have nonnegative Lagrange multipliers that are
! held at the beginning of VMULTC. The remainder of this vector holds the residuals of the inactive
! constraints at d, the ordering of the components of VMULTC being in agreement with the permutation
! of the indices of the constraints that is in IACT. All these residuals are nonnegative, which is
! achieved by the shift CSTRV that makes the least residual zero.


!??????????????????????????????? NACT <= min(M, N)?????????????????????????????????????????????????

      ifull = 1
      mcon = m
      stage = 1
      nact = 0
      cstrv = ZERO

! Initialize Z and some other variables. The value of CSTRV will be appropriate to D=0, while ICON
! will be the index of a most violated constraint if CSTRV is positive. Usually during the first
! stage the vector SDIRN gives a search direction that reduces all the active constraint violations
! by one simultaneously.
      z = eye(n, n)
      !!!! What is size (number of columns) of Z and related variables (CGZ, ZDOTA, ZDOTW)???????
      d = ZERO
      cstrv = maxval([b(1:m), ZERO])
      icon = maxloc(b(1:m), dim=1)
!iact(1:m) = [(k, k=1, m)]
      iact = [(k, k=1, m + 1)]
      ! What is the size of IACT? M or M + 1?
      vmultc(1:m) = cstrv - b(1:m)

      if (.not. (cstrv <= zero)) then
          call trstlp_sub(iact(1:m), ifull, stage, nact, A(:, 1:m), b(1:&
     &m), rho, cstrv, d, vmultc(1:m), z)
!-------------------------------------------------------------------------------------------------------!
!call trstlp_sub(iact(1:m), stage, nact, A(:, 1:m), b(1:m), rho, d, vmultc(1:m), z) ! Is this enough????
! Decide IFULL by ||D||
! Decide CSTRV by A and b
! Decide ZDOTA by Z and A
!-------------------------------------------------------------------------------------------------------!
      else
          ifull = 0
      end if
      mcon = m + 1
      stage = 2
      icon = mcon
      iact(mcon) = mcon
      vmultc(mcon) = ZERO
      optold = ZERO
      icount = 0_IK

      if (ifull == 0) then
          call trstlp_sub(iact, ifull, stage, nact, A, b, rho, cstrv, d,&
     & vmultc, z)
      end if

      end subroutine trstlp


      subroutine trstlp_sub(iact, ifull, stage, nact, A, b, rho, cstrv, &
     &d, vmultc, z)

! Generic modules
      use consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, TENTH, EPS, H&
     &UGENUM, DEBUGGING, SRNLEN
      use info_mod, only : FTARGET_ACHIEVED, MAXFUN_REACHED, TRSUBP_FAIL&
     &ED, SMALL_TR_RADIUS, NAN_X, NAN_INF_F
      use infnan_mod, only : is_nan, is_posinf, is_finite
      use debug_mod, only : errstop
      use output_mod, only : retmssg, rhomssg, fmssg
      use lina_mod, only : inprod, matprod, eye, planerot, isminor

      implicit none

!integer(IK), intent(in) :: n
!integer(IK), intent(in) :: m
      integer(IK), intent(in) :: stage
      real(RP), intent(in) :: A(:, :)
      !(n, m+1)
      real(RP), intent(in) :: b(:)
      real(RP), intent(in) :: rho
      real(RP), intent(inout) :: d(:)
      real(RP), intent(inout) :: cstrv
      real(RP), intent(inout) :: vmultc(:)
      real(RP), intent(inout) :: z(:, :)


      integer(IK), intent(inout) :: iact(:)
      integer(IK), intent(out) :: ifull
      integer(IK), intent(inout) :: nact


      real(RP) :: hypt
      real(RP) :: dnew(size(d))
      real(RP) :: cgrad(size(d))
      real(RP) :: cgz(size(d))
      real(RP) :: cgzabs(size(d))
      real(RP) :: cgzk
      real(RP) :: cgzkabs

      real(RP) :: alpha
      real(RP) :: beta
      real(RP) :: dd
      real(RP) :: grot(2, 2)
      real(RP) :: optnew
      real(RP) :: optold
      real(RP) :: ratio
      real(RP) :: sdirn(size(d))
      real(RP) :: cvold
      real(RP) :: sd
      real(RP) :: sp
      real(RP) :: spabs
      real(RP) :: ss
      real(RP) :: step
      real(RP) :: stpful
      real(RP) :: temp
      real(RP) :: tmpv(size(A, 2))
      real(RP) :: tmpvabs(size(A, 2))
      real(RP) :: vsave
      real(RP) :: zdota(size(z, 2))
      real(RP) :: zdotw
      real(RP) :: zdvabs
      real(RP) :: zdwabs
      real(RP) :: dold(size(d))
      ! N
      real(RP) :: vmultd(size(A, 2))
      ! Is this necessary?????
      integer(IK) :: i
      integer(IK) :: icon
      integer(IK) :: icount
      integer(IK) :: isave
      integer(IK) :: iter
      integer(IK) :: j
      integer(IK) :: k
      integer(IK) :: kk
      integer(IK) :: kl
      integer(IK) :: kw
      integer(IK) :: maxiter
      integer(IK) :: mcon
      integer(IK) :: m
      integer(IK) :: n
      integer(IK) :: nactx


      n = size(A, 1)
      mcon = size(A, 2)
      if (stage == 1) then
          m = mcon
      else
          m = mcon - 1
      end if

      if (stage == 1) then
          icon = maxloc(b(1:mcon), dim=1)
      else
          icon = mcon
! In Powell's code, stage 2 uses the ZDOTA calculated by stage 1. Here we re-calculate ZDOTA so
! that we do not need to pass it from stage 1 to stage 2 in order to reduce the coupling.
          zdota(1:nact) = [(inprod(z(:, i), A(:, iact(i))), i=1, nact)]
          temp = maxval([b(1:m) - matprod(d, A(:, 1:m)), ZERO])
          temp = abs(temp - cstrv) / max(1.0D0, abs(cstrv))
          if (temp > 1.0D2 * epsilon(1.0D0)) then
              open (unit=11, file='fort', status='old', position='append&
     &', action='write')
              write (11, *) 'cstrv', cstrv, temp
              write (11, *) - int(log10(temp))
              close (11)
          end if
      end if

      icount = 0_IK
      sdirn = ZERO
      optold = ZERO
      ifull = 1

! Powell's code can encounter infinite cycling, which did happen when testing the following CUTEst
! problems: DANWOODLS, GAUSS1LS, GAUSS2LS, GAUSS3LS, KOEBHELB, TAX13322, TAXR13322.
! Indeed, in all these cases, Inf and NaN appear in D due to extremely large values in A (up to
! 10^219). To resolve this, we set the maximal number of iterations to MAXITER, and terminate when
! NaN or Inf occurs in D.
!maxiter = min(10000_IK, 100_IK * max(m, n))
      maxiter = min(10000_IK, 100_IK * max(m, n))
      do iter = 1, maxiter
          if (is_finite(sum(abs(d)))) then
              dold = d
          else
              d = dold
              !! IFULL = ?????
              exit
          end if

          if (stage == 1) then
              optnew = cstrv
          else
              optnew = -inprod(d, A(:, mcon))
          end if

! End the current stage of the calculation if 3 consecutive iterations have either failed to
! reduce the best calculated value of the objective function or to increase the number of active
! constraints since the best value was calculated. This strategy prevents cycling, but there is
! a remote possibility that it will cause premature termination.
          if (icount == 0_IK .or. optnew < optold) then
              optold = optnew
              nactx = nact
              icount = 3_IK
          else if (nact > nactx) then
              nactx = nact
              icount = 3_IK
          else
              icount = icount - 1_IK
              if (icount == 0_IK) exit
          end if

! If ICON exceeds NACT, then we add the constraint with index IACT(ICON) to the active set.
! Apply Givens rotations so that the last N-NACT-1 columns of Z are orthogonal to the gradient
! of the new constraint, a scalar product being set to zero if its nonzero value could be due to
! computer rounding errors.
          if (icon > nact) then
              kk = iact(icon)
              ! IACT(ICON) will be updated later. So KK is necessary.
              cgrad = A(:, kk)
              cgz = matprod(cgrad, z)
              cgzabs = matprod(abs(cgrad), abs(z))
              where (isminor(cgz, cgzabs))
              ! In MATLAB: CGZ(ISMINOR(CGZ, CGZABS)) = ZERO
                  cgz = ZERO
              end where
              do k = n - 1, nact + 1, -1
! Apply a 2D Givens rotation to Z(:, [K,K+1]) from the right to zero CGRAD'*Z(:, K+1) out.
                  if (abs(cgz(k + 1)) > ZERO) then
! Powell wrote CGZ(K + 1) /= ZERO instead of ABS(CGZ(K + 1)) > ZERO. The two
! conditions differ if CGZ(K + 1) is NaN.
                      grot = planerot(cgz([k, k + 1]))
                      z(:, [k, k + 1]) = matprod(z(:, [k, k + 1]), trans&
     &pose(grot))
                      cgz(k) = sqrt(cgz(k)**2 + cgz(k + 1)**2)
                  end if
              end do

              if (nact < n .and. abs(cgz(nact + 1)) > ZERO) then
! Add the new constraint if this can be done without a deletion from the active set.
! Powell wrote "CGZ(NACT + 1) /= ZERO" instead of "ABS(CGZ(NACT + 1)) > ZERO", the two
! conditions differ if CGZ(NACT + 1) is NaN.
                  nact = nact + 1
                  zdota(nact) = cgz(nact)
                  vmultc(icon) = vmultc(nact)
                  vmultc(nact) = ZERO
              else
! The next instruction is reached if a deletion has to be made from the active set in
! order to make room for the new active constraint, because the new constraint gradient
! is a linear combination of the gradients of the old active constraints. Set the
! elements of VMULTD to the multipliers of the linear combination.
!----------- Powell wrote the following comment, but IOUT is never defined ------------!
! Further, set IOUT to the index of the constraint to be deleted, but branch if no
! suitable index can be found.
!--------------------------------------------------------------------------------------!
                  do k = nact, 1, -1
                  ! NACT >= N >= 1
! The following part used to lead to a memory error when VMULTD(K) is accessed (by why not Z(I, K) and others?????)
                      cgzk = inprod(cgrad, z(:, k))
                      cgzkabs = inprod(abs(cgrad), abs(z(:, k)))
                      if (isminor(cgzk, cgzkabs)) then
                          vmultd(k) = ZERO
                      else
                          vmultd(k) = cgzk / zdota(k)
                      end if
                      if (k >= 2) then
                          cgrad = cgrad - vmultd(k) * A(:, iact(k))
                      end if
                  end do

                  ratio = minval(vmultc(1:nact) / vmultd(1:nact), mask=(&
     &vmultd(1:nact) > ZERO .and. iact(1:nact) <= m))
                  if (ratio < ZERO .or. .not. any(vmultd(1:nact) > ZERO &
     &.and. iact(1:nact) <= m)) exit

! Revise the Lagrange multipliers and reorder the active constraints so that the one to
! be replaced is at the end of the list. Also calculate the new value of ZDOTA(NACT) and
! branch if it is not acceptable.
                  vmultc(1:nact) = max(ZERO, vmultc(1:nact) - ratio * vm&
     &ultd(1:nact))

!-------Powell's code includes the following, but it is IMPOSSIBLE TO REACH------------!
!if (icon < nact) then
!    do k = icon, nact - 1
!        hypt = sqrt(zdota(k + 1)**2 + inprod(z(:, k), A(:, iact(k + 1)))**2)
!        grot = planerot([zdota(k + 1), inprod(z(:, k), A(:, iact(k + 1)))])
!        z(:, [k, k + 1]) = matprod(z(:, [k + 1, k]), transpose(grot))
!        zdota([k, k + 1]) = [hypt, (zdota(k + 1) / hypt) * zdota(k)]
!    end do
!    iact(icon:nact) = [iact(icon + 1:nact), iact(icon)]
!    vmultc(icon:nact) = [vmultc(icon + 1:nact), vmultc(icon)]
!end if
!--------------------------------------------------------------------------------------!

                  if (inprod(z(:, nact), A(:, kk)) == ZERO) exit
                  ! Why?????
                  zdota(nact) = inprod(z(:, nact), A(:, kk))
                  ! Why ?????
                  vmultc(icon) = ZERO
                  vmultc(nact) = ratio
              end if

! Update IACT and ensure that the objective function continues to be treated as the last
! active constraint if stage 2 is in progress.
              iact([icon, nact]) = iact([nact, icon])
              if (stage == 2 .and. kk /= mcon) then
                  hypt = sqrt(zdota(nact)**2 + inprod(z(:, nact - 1), A(&
     &:, kk))**2)
                  ! What if HYPT = 0???
                  grot = planerot([zdota(nact), inprod(z(:, nact - 1), A&
     &(:, kk))])
                  z(:, [nact - 1, nact]) = matprod(z(:, [nact, nact - 1]&
     &), transpose(grot))
                  zdota([nact - 1, nact]) = [hypt, (zdota(nact) / hypt) &
     &* zdota(nact - 1)]
                  iact([nact - 1, nact]) = [kk, iact(nact - 1)]
                  vmultc([nact - 1, nact]) = vmultc([nact, nact - 1])
              end if

! If stage 1 is in progress, then set SDIRN to the direction of the next change to the
! current vector of variables.
              if (stage == 1) then
                  sdirn = sdirn - ((inprod(sdirn, A(:, iact(nact))) - ON&
     &E) / zdota(nact)) * z(:, nact)
              else
                  sdirn = (ONE / zdota(nact)) * z(:, nact)
              end if

          else

! Delete the constraint that has the index IACT(ICON) from the active set.
              if (icon < nact) then
                  do k = icon, nact - 1
                      hypt = sqrt(zdota(k + 1)**2 + inprod(z(:, k), A(:,&
     & iact(k + 1)))**2)
                      ! What if HYPT = 0 ???
                      grot = planerot([zdota(k + 1), inprod(z(:, k), A(:&
     &, iact(k + 1)))])
                      z(:, [k, k + 1]) = matprod(z(:, [k + 1, k]), trans&
     &pose(grot))
                      zdota([k, k + 1]) = [hypt, (zdota(k + 1) / hypt) *&
     & zdota(k)]
                  end do
                  iact(icon:nact) = [iact(icon + 1:nact), iact(icon)]
                  vmultc(icon:nact) = [vmultc(icon + 1:nact), vmultc(ico&
     &n)]
              end if
              nact = nact - 1

! If stage 1 is in progress, then set SDIRN to the direction of the next change to the
! current vector of variables.
              if (stage == 1) then
                  sdirn = sdirn - inprod(sdirn, z(:, nact + 1)) * z(:, n&
     &act + 1)
              else
                  sdirn = (ONE / zdota(nact)) * z(:, nact)
              end if
          end if

! Calculate the step to the boundary of the trust region or take the step  that reduces CSTRV to
! zero. The two statements below that include the factor EPS  prevent some harmless underflows
! that occurred in a test calculation (here, EPS is the machine epsilon; Powell's original code
! used 1.0E-6, and Powell's code was written in SINGLE PRECISION). Further, we skip the step if
! it could be zero within a reasonable tolerance for computer rounding errors.
          dd = rho**2 - sum(d**2, mask=(abs(d) >= EPS * rho))
          sd = inprod(sdirn, d)
          ss = inprod(sdirn, sdirn)
          if (dd <= ZERO) exit
          temp = sqrt(ss * dd)
          if (abs(sd) >= EPS * temp) then
              temp = sqrt(ss * dd + sd * sd)
          end if
          stpful = dd / (temp + sd)
          step = stpful
          if (stage == 1) then
              if (isminor(cstrv, step)) then
                  exit
              end if
              step = min(step, cstrv)
          end if

! Set DNEW to the new variables if STEP is the steplength, and reduce CSTRV to the corresponding
! maximum residual if stage 1 is being done. Because DNEW will be changed during the calculation
! of some Lagrange multipliers, it will be restored to the following value later.
          dnew = d + step * sdirn
          if (stage == 1) then
              cvold = cstrv
              cstrv = maxval([b(iact(1:nact)) - matprod(dnew, A(:, iact(&
     &1:nact))), ZERO])
          end if

! Set VMULTD to the VMULTC vector that would occur if D became DNEW. A device is included to
! force VMULTD(K)=ZERO if deviations from this value can be attributed to computer rounding
! errors. First calculate the new Lagrange multipliers.
          do k = nact, 1, -1
! What if NACT = 0? Is it possible? Powell's code will carry out the loop for one time.
              zdotw = inprod(z(:, k), dnew)
              zdwabs = inprod(abs(z(:, k)), abs(dnew))
! Powell's original code set ZDOTW = 0 when ISMINOR(ZDOTW, ZDWABS) = TRUE, and then take
! VMULTD(K) = ZDOTW/ZDOTA, which is NaN if ZDOTW = 0 = ZDOTA. The following code avoids NaN.
              if (isminor(zdotw, zdwabs)) then
                  vmultd(k) = ZERO
              else
                  vmultd(k) = zdotw / zdota(k)
              end if
              if (k >= 2) then
                  dnew = dnew - vmultd(k) * A(:, iact(k))
              end if
          end do
          if (stage == 2) then
              vmultd(nact) = max(ZERO, vmultd(nact))
          end if

! Complete VMULTD by finding the new constraint residuals.
          dnew = d + step * sdirn
          tmpv = matprod(dnew, A(:, iact)) - b(iact) + cstrv
          ! Indeed, only TMPV(nact+1:mcon) is needed.
          tmpvabs = matprod(abs(dnew), abs(A(:, iact))) + abs(b(iact)) +&
     & cstrv
          where (isminor(tmpv, tmpvabs))
              tmpv = ZERO
          end where
          vmultd(nact + 1:mcon) = tmpv(nact + 1:mcon)

! Calculate the fraction of the step from D to DNEW that will be taken.
          tmpv = vmultc / (vmultc - vmultd)
          !
          ratio = min(ONE, minval(tmpv(1:mcon), mask=(vmultd(1:mcon) < Z&
     &ERO)))
          if (ratio < ONE) then
              icon = minloc(tmpv(1:mcon), mask=(vmultd(1:mcon) < ZERO), &
     &dim=1)
          else
              icon = 0
          end if

! Update d, VMULTC and CSTRV.
          d = (ONE - ratio) * d + ratio * dnew
          vmultc(1:mcon) = max(ZERO, (ONE - ratio) * vmultc(1:mcon) + ra&
     &tio * vmultd(1:mcon))
          if (stage == 1) then
!cstrv = cvold + ratio * (cstrv - cvold)
              cstrv = (ONE - ratio) * cvold + ratio * cstrv
          end if

          if (icon == 0) then
              if (step >= stpful) then
              ! Indeed, STEP == STPFUL.
                  return
              else
                  exit
              end if
          end if
      end do
      ifull = 0
      end subroutine trstlp_sub

      end module trustregion_mod

! 1. Check the update of ZDOTA; is it always <Z(:, I), A(:, IACT(I))> ?
! 2. Check KK
! 3. Check CSTRV. Is it always max([b-A*x, 0])?
! 4. What is the objective of trstlp_sub?