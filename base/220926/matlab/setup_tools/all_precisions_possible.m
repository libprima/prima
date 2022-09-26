function [precision_list, default_precision] = all_precisions_possible()
%ALL_VARIANTS_POSSIBLE returns a cell array containing the names of all the possible precisions of
% the Fortran solvers in this package; `default_precision`, if requested, is the name of the default
% precision.

precision_list = {'double', 'single', 'quadruple'};

default_precision = 'double';

return
