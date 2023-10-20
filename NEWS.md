# New C API

The solver and problem definition options moved to a new `prima_options` struct
that must be initialized before use with `prima_init_options`,
and the results move to `prima_results`:
```
prima_options options;  
prima_init_options(&options);  
options.iprint = PRIMA_MSG_EXIT;  
options.rhoend= 1e-3;  
options.maxfun = 200*n;  
prima_results results;  
const int rc = prima_bobyqa(&fun, n, x, &options, &results);
prima_free_options(&options);  
prima_free_results(&results);
```
Both should be freed after use with their dedicated free function.
