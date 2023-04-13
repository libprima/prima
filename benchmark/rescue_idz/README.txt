This benchmarks tests how much difference the RESCUE of BOBYQA or IDZ of NEWUOA/LINCOA makes.

Tests: for SOLVER in {'bobyqa', 'newuoa', 'lincoa'},

options = struct();
options.classical = true or false;  % Test the classical version or not
prof(SOLVER, 'last', options);
prof(SOLVER, 'last', 'all', options);
prof(SOLVER, 'last', DIM_RANGE, options);
