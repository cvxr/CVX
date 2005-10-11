function [xir,xjc] = sdpasplit(xcol,maxnnz)
% ------------------------------------------------------------
% Allocate enough space
% ------------------------------------------------------------
 n = length(xcol);
 if nargin < 2
   maxnnz = n;
 end
 xir = zeros(maxnnz,1);
 xjc = zeros(maxnnz+1,1);
% ------------------------------------------------------------
% Find all nonempty blocks, and store their block start in xjc.
% ------------------------------------------------------------
 knz = 0;
 k = -1;
 for j = 1:n
   if xcol(j) > k
     knz = knz + 1;
     k = xcol(j);
     xir(knz) = k;
     xjc(knz) = j;
   end
 end
% ------------------------------------------------------------
% Close last block in xjc and shrink xjc, xir to correct size
% ------------------------------------------------------------
 xjc(knz+1) = n+1;
 xjc = xjc(1:knz+1);
 xir = xir(1:knz);
