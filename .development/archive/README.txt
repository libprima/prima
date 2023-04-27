This directory contains "the archive versions" of the package. For tests and developments only.

N.B.: When profiling the solvers against the "archive" version, e.g.,

```matlab
profile('newuoa', 'archive');
```

the "archive" version is the one in the `dev_arch` directory. This directory is a symlink pointing to a
recent archive version that we use as a benchmark for the development of the current version. It is not
necessarily the latest one.
