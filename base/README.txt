This directory contains "the base versions" of the package. For tests and developments only.

N.B.: When profiling the solvers against the "base" version, e.g.,

```matlab
profile('newuoa', 'base');
```

the "base" version is the one in the `devbase` directory. This directory is a symlink pointing to a
recent base version that we use as a benchmark for the development of the current version. It is not
necessarily the latest one.
