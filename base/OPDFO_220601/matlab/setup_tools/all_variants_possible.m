function [variant_list, default_variant] = all_variants_possible()
%ALL_VARIANTS_POSSIBLE returns a cell array containing the names of all the possible variants of
% the Fortran solvers in this package; `default_variant`, if requested, is the name of the default
% variant.

variant_list = {'modern', 'classical'};

default_variant = 'modern';

return
