This directory contains the "archiva" versions of the package. For tests and developments only.

N.B.: When profiling the solvers against the "archiva" version, e.g.,

```matlab
profile('newuoa', 'archiva');
```

the "archiva" version is the one in the `dev_arch` directory. This directory is a symlink pointing to a
recent archiva version that we use as a benchmark for the development of the current version. It is not
necessarily the latest one.
