To develop, `cd` into the `src` directory and run

```pip install --editable .```

This will install prima locally in an editable fashion. From there you can run the examples/cobyla/cobyla_example.py (from any directory) and go from there.

## Style notes

- In cases where there are Python or numpy functions that can achieve the same aim (i.e. sum vs np.sum) we prefer the numpy functions.
  - Rationale:
    - In some cases we must use the numpy function because it has functionality the Python function doesn't, i.e. np.sum has an axis argument, whereas sum does not.
    - In light of this, for the sake of consistency we prefer to use the numpy functions everywhere.
- Most of the comments are copied from Fortran verbatim, except in cases where they need to modified due to specifics of the Python language. In these cases a note will be made of the difference between Fortran and Python
  - Rationale:
      - The main purpose of this is to keep the Python and Fortran codebases as similar as possible.
- For determining the dimensions of an array, we exclusively use `np.size` instead of `np.shape` or `some_array.shape` or `len`
  - Rationale:
    - Fortran uses `SIZE` so this helps us to be as consistent with the Fortran code as possible.

## A note on Fortran's `max` and `maxval` and their Python equivalents

| Fortran | Python |
|---------|--------|
| `max`   | `np.maximum` |
| `maxval` | `np.max` or `max` |

Fortran's `max` and numpy's `maximum` accept two arguments, either of which can be a scalar or an array,
and returns an elementwise maximum of the two arguments. In the case of a scalar and an array argument it
returns an elementwise maximum of the scalar and each element of the array.

Fortran's `maxval` takes a single argument and returns the largest item in the list, similar to numpy's `max`.
Python's built-in `max` can either take a single argument and return the largest element, or it can take a
variable number of arguments and return the largest of all of them.

Per the style notes, prefer `np.max` over the built-in `max`.

For the larger of two numbers Fortran uses `max`, and we prefer `np.maximum` over the Python `max` for consistenchy in the translation.

This note applies to `min` and `minval` as well.


## A note on indices

Consider the following Fortran code

```
do i=0:5
  print *, *
end do
```

It can be easily and automatically translated to Python as

```
for i in range(0, 6):
  print(i)
```

Now consider the following similar loop

```
do i=1:5
  print *, some_array(i)
end do
```

This can be translated to Python as

```
for i in range(1, 6):
  print(some_array[i-1])
```

This leads to awkward Python code, since the more pythonic code would range from 0 to 5, and the indexing would be `some_array[i]`. In order to make the Python code more usable, we will attempt to write more "pythonic" code, even though that makes the translation a little bit more difficult.