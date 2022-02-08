! GETHUGE subroutine
!
! Authors:
!     Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
!     and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
!     Department of Applied Mathematics,
!     The Hong Kong Polytechnic University.
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! This function returns in huge_value the upper bound for the Fortran type data_type. If the given type is unknown,
! this function will return -huge(0.0D0) in huge_value.

subroutine gethuge(data_type, huge_value)
use pdfoconst
implicit none
! The largest length of the input string. It is 20 because the input should be 'integer', 'float', 'single', 'double', 'fun',
! 'function', 'con', 'constraint'
integer, parameter :: maxlen = 20
character(len=maxlen), intent(in) :: data_type ! the input string, which specifies the data type
double precision, intent(out) :: huge_value

integer, parameter :: intZero = 0 ! integer 0 of the default kind
real, parameter :: floatZero = 0.0 ! floating-point 0 of the default kind
real(kind=sp), parameter :: singleZero = 0.0 ! single-precision floating-point 0
double precision, parameter :: doubleZero = 0.0d0 ! double-precision floating-point 0

integer, parameter :: intHuge = huge(intZero)
real, parameter :: floatHuge = huge(floatZero)
real(kind=sp), parameter :: singleHuge = huge(singleZero)
double precision, parameter :: doubleHuge = huge(doubleZero)

if (data_type .eq. 'integer' .or. data_type .eq. 'Integer' .or. data_type .eq. 'INTEGER') then
    huge_value = real(intHuge, dp)
elseif (data_type .eq. 'float' .or. data_type .eq. 'real' .or. data_type .eq. 'Float' .or. data_type .eq. 'Real' .or. &
                data_type .eq. 'FLOAT' .or. data_type .eq. 'REAL') then
    huge_value = real(floatHuge, dp)
elseif (data_type .eq. 'single' .or. data_type .eq. 'Single' .or. data_type .eq. 'SINGLE') then
    huge_value = real(singleHuge, dp)
elseif (data_type .eq. 'double' .or. data_type .eq. 'Double' .or. data_type .eq. 'DOUBLE') then
    huge_value = doubleHuge
elseif (data_type .eq. 'fun' .or. data_type .eq. 'Fun' .or. data_type .eq. 'FUN' .or. data_type .eq. 'function' .or. &
                data_type .eq. 'Function' .or. data_type .eq. 'FUNCTION') then
    huge_value = real(HUGEFUN, dp) ! HUGEFUN is defined in pdfoconst
elseif (data_type .eq. 'con' .or. data_type .eq. 'Con' .or. data_type .eq. 'CON' .or. data_type .eq. 'constraint' .or. &
                data_type .eq. 'Constraint' .or. data_type .eq. 'CONSTRAINT') then
    huge_value = real(HUGECON, dp) ! HUGECON is defined in pdfoconst
else
    huge_value = -doubleHuge
endif

return
end subroutine gethuge
