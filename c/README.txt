This is the C interface for calling the Fortran version of PRIMA.

Contributed by Julien Schueller ( schueller@phimeca.com ) in September 2023.

## Usage
See https://github.com/libprima/prima#c for the usage.

## Developer notes
Since C does not provide a direct avenue for optional function arguments, we use the following values to indicate that a value is not present:

| type | value |
| ---- | ----- |
| integer | 0 |
| double  | NAN |
| array   | NULL |

The rationale for using 0 to indicate that an integer argument is not present lies in the fact that most of the time the code expects any integer arguments to be positive integers. See https://github.com/libprima/prima/pull/135 for a more detailed discussion.

For doubles, we set them to NAN but the check for NAN happens in the Fortran code using the `is_nan` function from the `infnan` module. See `rhobeg` and `rhoend` for an example.
