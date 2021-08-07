!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of geometry.f90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 07-Aug-2021.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      module geometry_mod

      contains

      function goodgeo(factor_alpha, factor_beta, rho, sim, simi) result&
     &(good_geo)

      use consts_mod, only : IK, RP, ONE, DEBUGGING, SRNLEN
      use debug_mod, only: errstop, verisize

      implicit none

! Inputs
      real(RP), intent(in) :: sim(:, :)
      real(RP), intent(in) :: simi(:, :)
      real(RP), intent(in) :: factor_alpha
      real(RP), intent(in) :: factor_beta
      real(RP), intent(in) :: rho

! Output
      logical :: good_geo

! Local variables
      integer(IK) :: n
      real(RP) :: parsig
      real(RP) :: pareta
      real(RP) :: vsig(size(sim, 1))
      real(RP) :: veta(size(sim, 1))
      character(len=SRNLEN) :: srname='goodgeo'

! Get and verify the sizes
      n = size(sim, 1)
      if (DEBUGGING) then
          if (n < 1) then
              call errstop(srname, 'SIZE(SIM, 1) < 1')
          end if
          call verisize(sim, n, n+1)
          call verisize(simi, n, n)
      end if

! Calculate the values of sigma and eta.
      parsig = factor_alpha * rho
      pareta = factor_beta * rho
! VSIG(J) (J=1, .., N) is The Euclidean distance from vertex J to the opposite face of
! the current simplex. But what about vertex N+1?
      vsig = ONE / sqrt(sum(simi**2, dim=2))
      veta = sqrt(sum(sim(:, 1:n)**2, dim=1))
      good_geo = all(vsig >= parsig) .and. all(veta <= pareta)

      end function goodgeo


      function setdrop_tr(actrem, d, factor_alpha, factor_delta, rho, si&
     &m, simi) result(jdrop)

      use consts_mod, only : IK, RP, ZERO, ONE, DEBUGGING, SRNLEN
      use lina_mod, only : matprod, inprod
      use debug_mod, only: errstop, verisize

      implicit none

! Inputs
      real(RP), intent(in) :: actrem
      real(RP), intent(in) :: d(:)
      real(RP), intent(in) :: factor_alpha
      real(RP), intent(in) :: factor_delta
      real(RP), intent(in) :: rho
      real(RP), intent(in) :: sim(:, :)
      real(RP), intent(in) :: simi(:, :)

! Output
      integer(IK) :: jdrop

! Local variables
      integer(IK) :: n
      real(RP) :: distx(size(sim, 1))
      real(RP) :: edgmax
      real(RP) :: parsig
      real(RP) :: ratio
      real(RP) :: sigbar(size(sim, 1))
      real(RP) :: simid(size(sim, 1))
      real(RP) :: vsig(size(sim, 1))
      character(len=SRNLEN) :: srname="setdrop_tr"

! Get and verify the sizes
      n = size(sim, 1)
      if (DEBUGGING) then
          if (n < 1) then
              call errstop(srname, 'SIZE(SIM, 1) < 1')
          end if
          call verisize(d, n)
          call verisize(sim, n, n+1)
          call verisize(simi, n, n)
      end if

      jdrop = 0_IK
      if (actrem <= ZERO) then
          ratio = ONE
      else
          ratio = ZERO
      end if
      simid = matprod(simi, d)
      if (maxval(abs(simid)) > ratio) then
          jdrop = int(maxloc(abs(simid), dim=1), kind(jdrop))
      end if

      if (actrem > ZERO) then
          distx = sqrt(sum((spread(d, dim=2, ncopies=n) - sim(:, 1:n))**&
     &2, dim=1))
      else
          distx = sqrt(sum(sim(:, 1:n)**2, dim=1))
      end if

      edgmax = factor_delta * rho
      parsig = factor_alpha * rho
      vsig = ONE / sqrt(sum(simi**2, dim=2))
      sigbar = abs(simid) * vsig
      if (any(distx > edgmax .and. (sigbar >= parsig .or. sigbar >= vsig&
     &))) then
          jdrop = int(maxloc(distx, mask=(sigbar >= parsig .or. sigbar >&
     &= vsig), dim=1), kind(jdrop))
      end if

      end function setdrop_tr



      function setdrop_geo(factor_alpha, factor_beta, rho, sim, simi) re&
     &sult(jdrop)

      use consts_mod, only : IK, RP, ONE, DEBUGGING, SRNLEN
      use debug_mod, only: errstop, verisize

      implicit none

! Inputs
      real(RP), intent(in) :: sim(:, :)
      real(RP), intent(in) :: simi(:, :)
      real(RP), intent(in) :: factor_alpha
      real(RP), intent(in) :: factor_beta
      real(RP), intent(in) :: rho

! Output
      integer(IK) :: jdrop

! Local variables
      integer(IK) :: n
      real(RP) :: parsig
      real(RP) :: pareta
      real(RP) :: vsig(size(sim, 1))
      real(RP) :: veta(size(sim, 1))
      character(len=SRNLEN) :: srname='setdrop_geo'

! Get and verify the sizes
      n = size(sim, 1)
      if (DEBUGGING) then
          if (n < 1) then
              call errstop(srname, 'SIZE(SIM, 1) < 1')
          end if
          call verisize(sim, n, n+1)
          call verisize(simi, n, n)
      end if

! Calculate the values of sigma and eta.
      parsig = factor_alpha * rho
      pareta = factor_beta * rho
! VSIG(J) (J=1, .., N) is The Euclidean distance from vertex J to the opposite face of
! the current simplex. But what about vertex N+1?
      vsig = ONE / sqrt(sum(simi**2, dim=2))
      veta = sqrt(sum(sim(:, 1:n)**2, dim=1))

! Decide which vertex to drop from the simplex. It will be replaced by a new point to improve
! acceptability of the simplex. See equations (15) and (16) of the COBYLA paper.
      if (maxval(veta) > pareta) then
          jdrop = int(maxloc(veta, dim=1), kind(jdrop))
      else
          jdrop = int(minloc(vsig, dim=1), kind(jdrop))
      end if

      end function setdrop_geo


      function geostep(jdrop, cpen, datmat, factor_gamma, rho, simi) res&
     &ult(d)

      use consts_mod, only : IK, RP, ZERO, ONE, TWO, DEBUGGING, SRNLEN
      use lina_mod, only : matprod, inprod
      use debug_mod, only: errstop, verisize

      implicit none

! Inputs
      integer(IK), intent(in) :: jdrop
      real(RP), intent(in) :: simi(:, :)
      real(RP), intent(in) :: factor_gamma
      real(RP), intent(in) :: cpen
      real(RP), intent(in) :: datmat(:, :)
      real(RP), intent(in) :: rho

! Output
      real(RP) :: d(size(simi, 1))

! Local variables
      integer(IK) :: m
      integer(IK) :: n
      real(RP) :: cvmaxp
      real(RP) :: cvmaxm
      real(RP) :: vsig(size(simi, 1))
      real(RP) :: A(size(simi, 1), size(datmat, 1) - 1)
      character(len=SRNLEN) :: srname='geostep'

! Get and verify the sizes
      m = size(datmat, 1) - 2
      n = size(datmat, 2) - 1
      if (DEBUGGING) then
          if (m < 0 .or. n < 1) then
              call errstop(srname, 'SIZE(DATMAT) is invalid')
          end if
          call verisize(simi, n, n)
      end if

! VSIG(J) (J=1, .., N) is The Euclidean distance from vertex J to the opposite face of
! the current simplex. But what about vertex N+1?
      vsig = ONE / sqrt(sum(simi**2, dim=2))
      d = factor_gamma * rho * vsig(jdrop) * simi(jdrop, :)
! Calculate the coefficients of the linear approximations to the objective and constraint functions,
! placing minus the objective function gradient after the constraint gradients in the array A.
! When __USE_INTRINSIC_ALGEBRA__ = 1, the following code may not produce the same result as
! Powell's, because the intrinsic MATMUL behaves differently from a naive triple loop in
! finite-precision arithmetic.
! Is it more reasonable to save A transpose instead of A? Better name for A?
      A = transpose(matprod(datmat(1:m + 1, 1:n) - spread(datmat(1:m + 1&
     &, n + 1), dim=2, ncopies=n), simi))
      A(:, m+1) = - A(:, m+1)
      cvmaxp = maxval([ZERO, -matprod(d, A(:, 1:m)) - datmat(1:m, n + 1)&
     &])
      cvmaxm = maxval([ZERO, matprod(d, A(:, 1:m)) - datmat(1:m, n + 1)]&
     &)
      if (TWO * inprod(d, A(:, m + 1)) < cpen * (cvmaxp - cvmaxm)) then
          d = -d
      end if

      end function geostep


      end module geometry_mod