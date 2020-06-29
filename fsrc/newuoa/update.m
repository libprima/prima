function [bmat, zmat, idz] = update(bmat, zmat, idz, vlag, beta, knew)

n = size(bmat, 1);
npt = size(bmat, 2) - n;
