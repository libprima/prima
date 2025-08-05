This benchmarks tests how much difference the RESCUE of BOBYQA or IDZ of NEWUOA/LINCOA makes.

Tests: for SOLVER in {'bobyqa', 'newuoa', 'lincoa'},

options = struct();
options.classical = true or false;  % Test the classical version or not
prof(SOLVER, 'norma', options);
prof(SOLVER, 'norma', 'all', options);
prof(SOLVER, 'norma', DIM_RANGE, options);
