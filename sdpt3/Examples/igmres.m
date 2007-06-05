%%********************************************************************
%% igmres: compute ideal GMRES polynomial of deg m
%%         of a matrix B.
%%
%%    min || p(B) ||_2 
%%    where p is a polynomial of degree <= m and p(0) = 1. 
%%
%%    B = a general square matrix.
%%----------------------------------------------------------- 
%%
%% [blk,Avec,C,b,X0,y0,Z0,objval,p] = igmres(B,m,feas,solve); 
%%
%% B = a square matrix.
%% m = degree of polynomial. 
%% feas  = 1 if want feasible starting point
%%       = 0 if otherwise.
%% solve = 0 just to initialize
%%       = 1 if want to solve the problem.
%%
%% p = ideal GMRES polynomial in Matlab format.
%%
%% SDPT3: version 3.0 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last modified: 2 Feb 01
%%********************************************************************

 function [blk,Avec,C,b,X0,y0,Z0,objval,p] = igmres(B,m,feas,solve); 
 
   if (nargin <= 3); solve = 0; end;
   if (nargin <= 2); feas  = 0; end; 

   N = length(B); 
   if (m >= N); error('degree >= size of B'); end;

   A = cell(1,m+1); 
   [V,H,R] = orthbasis(B,m-1);
   A{1} = eye(N);
   for k = [1:m];  A{k+1} = B*V{k};  end; 
   [blk,Avec,C,b,X0,y0,Z0,objval,xx] = norm_min(A,feas,solve); 
   if (solve)
      y = xx; 
      x1 = R(1:m,1:m)*y(1:m);
      p = [x1(m:-1:1); 1];  
   else 
      objval = []; p = [];
   end
%%********************************************************************
