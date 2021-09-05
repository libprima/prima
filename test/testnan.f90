program testnan
use, intrinsic :: ieee_arithmetic
write (*, *) kind(.true.)
write (*, *) ieee_is_nan(1.0)
write (*, *) kind(ieee_is_nan(1.0))
end program testnan
