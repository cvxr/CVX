%%********************************************************************* 
%% chebymat: compute the Chebyshev polynomial of 
%%           deg m of a matrix B.
%%        
%% [blk,Avec,C,b,X0,y0,Z0,objval,p] = chebymat(B,m,feas,solve);
%%
%% B = a square matrix.
%% m = degree of polynomial. 
%% feas  = 1 if want feasible starting point
%%       = 0 if otherwise.
%% solve = 0 just to initialize
%%       = 1 if want to solve using sdp.m
%%       = 2 if want to solve using sdphlf.m
%%
%% p = Chebyshev polynomial of B in Matlab format.
%%
%% SDPT3: version 3.0 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last modified: 2 Feb 01
%%********************************************************************* 

 function [blk,Avec,C,b,X0,y0,Z0,objval,p] = chebymat(B,m,feas,solve);

   if (nargin <= 3); solve = 0; end;
   if (nargin <= 2); feas  = 0;  end; 

   N = length(B); 
   if (m >= N); error('degree >= size of B'); end;

   A = cell(1,m+1); 
   [V,H,R] = orthbasis(B,m); 
   A{1} = V{m+1};
   for k = 2:m+1;  A{k} = V{k-1};  end;
   [blk,Avec,C,b,X0,y0,Z0,objval,xx] = norm_min(A,feas,solve);
%%
   if (solve) 
      y  = -xx; 
      h  = R(1:m+1,1:m+1) \ [zeros(m,1); 1]; 
      x1 = R(1:m,1:m)*(h(m+1)*y(1:m) + h(1:m));
      p  = [1; -x1(m:-1:1)]; 
      objval = norm(polyvalm(p,B));
   else 
      objval = []; p = [];
   end 
%%********************************************************************* 
