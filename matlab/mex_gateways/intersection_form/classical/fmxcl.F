!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of fmxcl.F90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 09-Aug-2020.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! FMXCL_MOD is a module that provides the following procedures:
!
! fmxAllocate
! fmxReadMPtr
! fmxWriteMPtr
!
! They are analogous to the subroutines with the same name in FMXAPI_MOD.
! The difference is that the procesures in FMXAPI_MOD deal with
! reading/writing of INTEGER(IK) and REAL(RP), while the procesures
! here read/write integers and reals with the "classical" kinds:
!
! classical kind for INTEGER: IK_CL = IK_DFT, e.g., kind(0);
! classical kind for REAL: RP_CL = DP, e.g., kind(0.0D0).
!
! N.B.:
! 0. IK_CL and RP_CL can be changed easily if needed (see below),
!    but we decide not to support RP_CL = REAL128.
! 1. These procedures are needed when interfacing the "classical mode"
!    of Powell's Fortran code with MEX.
! 2. Why cannot we emerge the fmxReadMPtr and fmxWriteMPtr here with
!    those in FMXAPI_MOD? Because we cannot decide at compilation time
!    whether IK = IK_CL and whether RP = RP_CL (we can if we
!    include ppf.h, but we do not want to do that). Consequently, it is
!    undeciable whether the interfaces of fmxReadMPtr should or not
!    include the subroutines read_rscalar_cl, read_iscalar, etc. A
!    solution is to let fmxReadMPtr include all possible combinations of
!    integer kind and real kind (write_rscalar_sp, write_rscalar_dp,
!    write_rscalar_qp, write_iscalar_int16_sp, write_iscalar_int16_dp,
!    ...), but there are too many combinations!
!    The same applies to fmxWriteMPtr.
! 3. We decide to name the procedures in exactly the same way as in
!    FMXAPI_MOD so that the classical and normal modes of the MEX
!    gateways can have exactly the same I/O code, although they use
!    fmxReadMPtr/fmxWriteMPtr from different modules.

! Coded by Zaikun ZHANG in July 2020.


#include "fintrf.h"


      module fmxcl_mod

      use consts_mod, only : IK_CL => IK_DFT, RP_CL => DP
      use fmxapi_mod, only : convTol, mwOne, notComplex
      use fmxapi_mod, only : mxGetM, mxGetN, mexErrMsgIdAndTxt
      use fmxapi_mod, only : mxCreateDoubleScalar
      use fmxapi_mod, only : mxCreateDoubleMatrix
      use fmxapi_mod, only : mxCopyPtrToReal8, mxCopyReal8ToPtr
      use fmxapi_mod, only : fmxVerifyClassShape, fmxGetDble
      implicit none
      private
      public :: IK_CL
      public :: RP_CL
      public :: fmxAllocate
      public :: fmxReadMPtr
      public :: fmxWriteMPtr

      real(RP_CL), parameter :: ONE = 1.0_RP_CL

      interface fmxAllocate
! fmxAllocate allocates the space for a vector/matrix
          module procedure alloc_rvector_cl_sp, alloc_rmatrix_cl_sp
          module procedure alloc_rvector_cl_dp, alloc_rmatrix_cl_dp
      end interface fmxAllocate

      interface fmxReadMPtr
! fmxReadMPtr reads the numeric data associated with an mwPointer.
! It verifies the class and shape of the data before reading, and
! converts the data to REAL(RP_CL) if necessary.
          module procedure read_rscalar_cl, read_rvector_cl, read_rmatri&
     &x_cl
          module procedure read_iscalar_cl
      end interface fmxReadMPtr

      interface fmxWriteMPtr
! fmxWriteMPtr associates numeric data with an mwPointer. It converts
! the data to REAL(DP) if necessary, and allocates space if the
! data is a vector or matrix. Therefore, it is necessary to call
! mxDestroyArray when the usage of the vector/matrix terminates.
          module procedure write_rscalar_cl, write_rmatrix_cl, write_rve&
     &ctor_cl
          module procedure write_iscalar_cl
      end interface fmxWriteMPtr


      contains


      subroutine alloc_rvector_cl_sp(x, n)
! ALLOC_RVECTOR_CL_SP allocates the space for an allocatable
! single-precision vector X, whose size is N after allocation.
      use consts_mod, only : SP, MSSGLEN
      implicit none

! Input
      integer(IK_CL), intent(in) :: n

! Output
      real(SP), allocatable, intent(out) :: x(:)

! Intermediate variable
      integer :: alloc_status
      character(len = MSSGLEN) :: eid, mssg

! Deallocate X if X is allocated. This is necessary.
      if (allocated(x)) then
          deallocate (x)
      end if

! Allocate memory for X
      allocate (x(n), stat = alloc_status)
      if (alloc_status /= 0) then
          eid = 'FMXAPI:AllocateFailed'
          mssg = 'ALLOC_RVECTOR_CL: Memory allocation fails.'
          call mexErrMsgIdAndTxt(trim(eid), trim(mssg))
      end if

! Set X = 0; otherwise, compilers may complain.
      x = 0.0_SP

      end subroutine alloc_rvector_cl_sp


      subroutine alloc_rmatrix_cl_sp(x, m, n)
! ALLOC_RMATRIX_CL_SP allocates the space for a single-precision
! matrix X, whose size is (M, N) after allocation.
      use consts_mod, only : SP, MSSGLEN
      implicit none

! Input
      integer(IK_CL), intent(in) :: m, n

! Output
      real(SP), allocatable, intent(out) :: x(:, :)

! Intermediate variable
      integer :: alloc_status
      character(len = MSSGLEN) :: eid, mssg

! Deallocate X if X is allocated. This is necessary.
      if (allocated(x)) then
          deallocate (x)
      end if

! Allocate memory for X
      allocate (x(m, n), stat = alloc_status)
      if (alloc_status /= 0) then
          eid = 'FMXAPI:AllocateFailed'
          mssg = 'ALLOC_RMATRIX_CL: Memory allocation fails.'
          call mexErrMsgIdAndTxt(trim(eid), trim(mssg))
      end if

! Set X = 0; otherwise, compilers may complain.
      x = 0.0_SP
      end subroutine alloc_rmatrix_cl_sp


      subroutine alloc_rvector_cl_dp(x, n)
! ALLOC_RVECTOR_CL_DP allocates the space for an allocatable
! single-precision vector X, whose size is N after allocation.
      use consts_mod, only : DP, MSSGLEN
      implicit none

! Input
      integer(IK_CL), intent(in) :: n

! Output
      real(DP), allocatable, intent(out) :: x(:)

! Intermediate variable
      integer :: alloc_status
      character(len = MSSGLEN) :: eid, mssg

! Deallocate X if X is allocated. This is necessary.
      if (allocated(x)) then
          deallocate (x)
      end if

! Allocate memory for X
      allocate (x(n), stat = alloc_status)
      if (alloc_status /= 0) then
          eid = 'FMXAPI:AllocateFailed'
          mssg = 'ALLOC_RVECTOR_CL: Memory allocation fails.'
          call mexErrMsgIdAndTxt(trim(eid), trim(mssg))
      end if

! Set X = 0; otherwise, compilers may complain.
      x = 0.0_DP

      end subroutine alloc_rvector_cl_dp


      subroutine alloc_rmatrix_cl_dp(x, m, n)
! ALLOC_RMATRIX_CL_DP allocates the space for a single-precision
! matrix X, whose size is (M, N) after allocation.
      use consts_mod, only : DP, MSSGLEN
      implicit none

! Input
      integer(IK_CL), intent(in) :: m, n

! Output
      real(DP), allocatable, intent(out) :: x(:, :)

! Intermediate variable
      integer :: alloc_status
      character(len = MSSGLEN) :: eid, mssg

! Deallocate X if X is allocated. This is necessary.
      if (allocated(x)) then
          deallocate (x)
      end if

! Allocate memory for X
      allocate (x(m, n), stat = alloc_status)
      if (alloc_status /= 0) then
          eid = 'FMXAPI:AllocateFailed'
          mssg = 'ALLOC_RMATRIX_CL: Memory allocation fails.'
          call mexErrMsgIdAndTxt(trim(eid), trim(mssg))
      end if

! Set X = 0; otherwise, compilers may complain.
      x = 0.0_DP
      end subroutine alloc_rmatrix_cl_dp


      subroutine read_rscalar_cl(px, x)
! READ_RSCALAR_CL reads the double scalar associated with an mwPointer
! PX and saves the data in X, which is a REAL(RP_CL) scalar.
      use consts_mod, only : DP, MSSGLEN
      implicit none

! Input
      mwPointer, intent(in) :: px

! Output
      real(RP_CL), intent(out) :: x

! Intermediate variable
      real(DP) :: x_dp(1)
      character(len = MSSGLEN) :: eid, mssg

! Check input type and size
      call fmxVerifyClassShape(px, 'double', 'scalar')

! Read the input
      call mxCopyPtrToReal8(fmxGetDble(px), x_dp, mwOne)

! Convert the input to the type expected by the Fortran code
      x = real(x_dp(1), kind(x))
! Check whether the type conversion is proper
      if (kind(x) /= kind(x_dp)) then
          if (abs(x-x_dp(1)) > convTol*max(abs(x), ONE)) then
              eid = 'FMXAPI:ConversionError'
              mssg = 'READ_RSCALAR_CL: Large error occurs when convertin&
     &g REAL(DP) to REAL(RP_CL).'
              call mexErrMsgIdAndTxt(trim(eid), trim(mssg))
          end if
      end if
      end subroutine read_rscalar_cl


      subroutine read_rvector_cl(px, x)
! READ_RVECTOR_CL reads the double vector associated with an mwPointer
! PX and saves the data in X, which is a REAL(RP_CL) vector.
      use consts_mod, only : DP, MSSGLEN
      implicit none

! Input
      mwPointer, intent(in) :: px

! Output
      real(RP_CL), allocatable, intent(out) :: x(:)

! Intermediate variables
      real(DP), allocatable :: x_dp(:)
      integer(IK_CL) :: n
      mwSize :: n_mw
      character(len = MSSGLEN) :: eid, mssg

! Check input type and size
      call fmxVerifyClassShape(px, 'double', 'vector')

! Get size
      n_mw = int(mxGetM(px)*mxGetN(px), kind(n_mw))
      n = int(n_mw, kind(n))

! Allocate memory for X and X_DP
      call fmxAllocate(x, n)
      call fmxAllocate(x_dp, n)

! Copy input to X_DP
      call mxCopyPtrToReal8(fmxGetDble(px), x_dp, n_mw)

! Convert X_DP to the type expected by the Fortran code
      x = real(x_dp, kind(x))
! Check whether the type conversion is proper
      if (kind(x) /= kind(x_dp)) then
          if (maxval(abs(x-x_dp)) > convTol*max(maxval(abs(x)), ONE)) th&
     &en
              eid = 'FMXAPI:ConversionError'
              mssg = 'READ_RVECTOR_CL: Large error occurs when convertin&
     &g REAL(DP) to REAL(RP_CL).'
              call mexErrMsgIdAndTxt(trim(eid), trim(mssg))
          end if
      end if

! Deallocate X_DP
      deallocate(x_dp)

      end subroutine read_rvector_cl


      subroutine read_rmatrix_cl(px, x)
! READ_MATRIX_CL reads the double matrix associated with an mwPointer
! PX and saves the data in X, which is a REAL(RP_CL) matrix.
      use consts_mod, only : DP, MSSGLEN
      implicit none

! Input
      mwPointer, intent(in) :: px

! Output
      real(RP_CL), allocatable, intent(out) :: x(:, :)

! Intermediate variables
      real(DP), allocatable :: x_dp(:, :)
      integer(IK_CL) :: m, n
      mwSize :: xsize
      character(len = MSSGLEN) :: eid, mssg

! Check input type and size
      call fmxVerifyClassShape(px, 'double', 'matrix')

! Get size
      m = int(mxGetM(px), kind(m))
      n = int(mxGetN(px), kind(n))
      xsize = int(m*n, kind(xsize))


! Allocate memory for X and X_DP
      call fmxAllocate(x, m, n)
      call fmxAllocate(x_dp, m, n)

! Copy input to X_DP
      call mxCopyPtrToReal8(fmxGetDble(px), x_dp, xsize)

! Convert X_DP to the type expected by the Fortran code
      x = real(x_dp, kind(x))
! Check whether the type conversion is proper
      if (kind(x) /= kind(x_dp)) then
          if (maxval(abs(x-x_dp)) > convTol*max(maxval(abs(x)), ONE)) th&
     &en
              eid = 'FMXAPI:ConversionError'
              mssg = 'READ_MATRIX_CL: Large error occurs when converting&
     & REAL(DP) to REAL(RP_CL).'
              call mexErrMsgIdAndTxt(trim(eid), trim(mssg))
          end if
      end if

! Deallocate X_DP
      deallocate(x_dp)

      end subroutine read_rmatrix_cl


      subroutine read_iscalar_cl(px, x)
! READ_ISCALAR_CL reads a MEX input X that is a double scalar with an
! integer value. Such a value will be passed to the Fortran code as an
! integer but passed by MEX as a double.
      use consts_mod, only : DP, MSSGLEN
      implicit none

! Input
      mwPointer, intent(in) :: px

! Output
      integer(IK_CL), intent(out) :: x

! Intermediate variable
      real(DP) :: x_dp(1)
      character(len = MSSGLEN) :: eid, mssg

! Check input type and size
      call fmxVerifyClassShape(px, 'double', 'scalar')

! Read the input
      call mxCopyPtrToReal8(fmxGetDble(px), x_dp, mwOne)

! Convert the input to the type expected by the Fortran code
      x = int(x_dp(1), kind(x))

! Check whether the type conversion is proper
      if (abs(x-x_dp(1)) > 0.5_DP) then
          eid = 'FMXAPI:ConversionError'
          mssg = 'READ_ISCALAR_CL: Large error occurs when converting RE&
     &AL(DP) to INTEGER(IK_CL).'
          call mexErrMsgIdAndTxt(trim(eid), trim(mssg))
      end if
      end subroutine read_iscalar_cl


      subroutine write_rscalar_cl(x, px)
! WRITE_RSCALAR_CL associates a REAL(RP_CL) scalar X with an mwPointer PX,
! after which X can be passed to MATLAB either as an output of
! mexFunction or an input of mexCallMATLAB.
      use consts_mod, only : DP, MSSGLEN
      implicit none

! Input
      real(RP_CL), intent(in) :: x

! Output
      mwPointer, intent(out) :: px

! Intermediate variable
      real(DP) :: x_dp
      character(len = MSSGLEN) :: eid, mssg

! Convert X to REAL(DP), which is expected by mxCopyReal8ToPtr
      x_dp = real(x, kind(x_dp))
      if (kind(x_dp) /= kind(x)) then
! Check whether the type conversion is proper
          if (abs(x-x_dp) > convTol*max(abs(x), ONE)) then
              eid = 'FMXAPI:ConversionError'
              mssg = 'WRITE_RSCALAR_CL: Large error occurs when converti&
     &ng REAL(RP_CL) to REAL(DP).'
              call mexErrMsgIdAndTxt(trim(eid), trim(mssg))
          end if
      end if

      px = mxCreateDoubleScalar(x_dp)

      end subroutine write_rscalar_cl


      subroutine write_rvector_cl(x, px, rowcol)
! WRITE_RVECTOR_CL associates a REAL(RP_CL) vector X with an
! mwPointer PX, after which X can be passed to MATLAB either as an output
! of mexFunction or an input of mexCallMATLAB. If ROWCOL = 'row', then
! the vector is passed as a row vector, otherwise, it will be a column vector.
      use consts_mod, only : DP, MSSGLEN
      implicit none

! Input
      real(RP_CL), intent(in) :: x(:)
      character(len = *), intent(in), optional :: rowcol

! Output
      mwPointer, intent(out) :: px

! Intermediate variable
      real(DP) :: x_dp(size(x))
      integer(IK_CL) :: n
      mwSize :: n_mw
      logical :: row
      character(len = MSSGLEN) :: eid, mssg

! Get size of X
      n_mw = int(size(x), kind(n_mw))
      n = int(n_mw, kind(n))

! Convert X to REAL(DP), which is expected by mxCopyReal8ToPtr
      x_dp = real(x, kind(x_dp))
! Check whether the type conversion is proper
      if (kind(x) /= kind(x_dp)) then
          if (maxval(abs(x-x_dp)) > convTol*max(maxval(abs(x)), ONE)) th&
     &en
              eid = 'FMXAPI:ConversionError'
              mssg = 'WRITE_RVECTOR_CL: Large error occurs when converti&
     &ng REAL(RP_CL) to REAL(DP).'
              call mexErrMsgIdAndTxt(trim(eid), trim(mssg))
          end if
      end if

      row = .false.
      if (present(rowcol)) then
          if (rowcol == 'row' .or. rowcol == 'ROW' .or. rowcol =='Row') &
     &then
              row = .true.
          end if
      end if
! Create a MATLAB matrix using the data in X_DP
      if (row) then
          px = mxCreateDoubleMatrix(mwOne, n_mw, notComplex)
      else
          px = mxCreateDoubleMatrix(n_mw, mwOne, notComplex)
      end if
      call mxCopyReal8ToPtr(x_dp, fmxGetDble(px), n_mw)

      end subroutine write_rvector_cl


      subroutine write_rmatrix_cl(x, px)
! WRITE_MATRIX_CL associates a REAL(RP_CL) matrix X with an mwPointer
! PX, after which X can be passed to MATLAB either as an output of
! mexFunction or an input of mexCallMATLAB.
      use consts_mod, only : DP, MSSGLEN
      implicit none

! Input
      real(RP_CL), intent(in) :: x(:, :)

! Output
      mwPointer, intent(out) :: px

! Intermediate variable
      real(DP) :: x_dp(size(x, 1), size(x, 2))
      integer(IK_CL) :: m, n
      mwSize :: m_mw, n_mw
      character(len = MSSGLEN) :: eid, mssg

! Get size of X
      m = int(size(x, 1), kind(m))
      n = int(size(x, 2), kind(n))
      m_mw = int(m, kind(m_mw))
      n_mw = int(n, kind(n_mw))

! Convert X to REAL(DP), which is expected by mxCopyReal8ToPtr
      x_dp = real(x, kind(x_dp))
! Check whether the type conversion is proper
      if(kind(x) /= kind(x_dp)) then
          if (maxval(abs(x-x_dp)) > convTol*max(maxval(abs(x)), ONE)) th&
     &en
              eid = 'FMXAPI:ConversionError'
              mssg = 'WRITE_MATRIX_CL: Large error occurs when convertin&
     &g REAL(RP_CL) to REAL(DP).'
              call mexErrMsgIdAndTxt(trim(eid), trim(mssg))
          end if
      end if

! Create a MATLAB matrix using the data in X_DP
      px = mxCreateDoubleMatrix(m_mw, n_mw, notComplex)
      call mxCopyReal8ToPtr(x_dp, fmxGetDble(px), m_mw*n_mw)

      end subroutine write_rmatrix_cl


      subroutine write_iscalar_cl(x, px)
! WRITE_RSCALAR_CL associates an INTEGER(IK_CL) scalar X with an
! mwPointer PX, after which X can be passed to MATLAB either as an output
! of mexFunction or an input of mexCallMATLAB.
      use consts_mod, only : DP, MSSGLEN
      implicit none

! Input
      integer(IK_CL), intent(in) :: x

! Output
      mwPointer, intent(out) :: px

! Intermediate variable
      real(DP) :: x_dp
      character(len = MSSGLEN) :: eid, mssg

      x_dp = real(x, kind(x_dp))

      if (abs(x - x_dp) > 0.5_DP) then
          eid = 'FMXAPI:ConversionError'
          mssg = 'WRITE_ISCALAR_CL: Large error occurs when converting I&
     &NTEGER(IK_CL) to REAL(DP).'
          call mexErrMsgIdAndTxt(trim(eid), trim(mssg))
      end if

      px = mxCreateDoubleScalar(x_dp)

      end subroutine write_iscalar_cl


      end module fmxcl_mod