function hpa = half_precision_available()
%HALF_PRECISION_SUPPORTED checks whether half precision is supported

hpa = ismac_silicon();
