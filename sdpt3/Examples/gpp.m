%%*******************************************************
%% gpp: pgraph partitioning problem.
%%
%% (primal problem) min   Tr C*X
%%                  s.t.  Tr(ee'*X) = alpha,      
%%                        diag(X) = e, 
%%                      
%%   C = -(diag(B*e)-B). 
%%
%% (dual problem)   max   alpha*y1 + e'*y 
%%                  s.t. y1*e*e^T + diag(y) + Z = C. 
%%-------------------------------------------------------
%%
%% [blk,Avec,C,b,X0,y0,Z0,objval,X] = gpp(B,alpha,feas,solve);
%%
%% B: weighted adjacency matrix of a graph with n nodes.
%% alpha: any real number,
%%        for alpha in (0,n^2), primal problem is strictly feasible.
%%                  in {0,n^2}, primal problem is feasible but not strictly.
%%                  outside [0,n^2], primal problem is infeasible.
%%        [default = 1]. 
%% feas  = 1 if want feasible starting point
%%       = 0 if otherwise.
%% solve = 0 just to initialize
%%       = 1 if want to solve the problem. 
%%
%% See graph.m --- generate random adjacency matrix. 
%%
%% SDPT3: version 3.0 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last modified: 2 Feb 01
%%*******************************************************

 function [blk,Avec,C,b,X0,y0,Z0,objval,X] = gpp(B,alpha,feas,solve);

   if nargin < 2; alpha = 1; end;
   if nargin < 3; feas  = 0; end; 
   if nargin < 4; solve = 0; end; 
   if ~isreal(B); error('only real B allowed'); end; 
 
   n = length(B); e = ones(n,1); 
   C = -(spdiags(B*e,0,n,n)-B); 
   b = [alpha; e];
   blk{1,1} = 's';  blk{1,2} = n;
  
   A = cell(1,n+1);
   A{1} = e*e'; 
   for k = 1:n; A{k+1} = sparse(k,k,1,n,n); end; 

   Avec = svec(blk,A,ones(size(blk,1),1)); 
   if (feas == 1);
      error('feas = 1, this option is not avaliable'); 
   elseif (feas == 0);  
      [X0,y0,Z0] = infeaspt(blk,Avec,C,b); 
   end
   if (solve) 
      [obj,X,y,Z] = sqlp(blk,Avec,C,b,[],X0,y0,Z0);
      objval = obj(1);
   else
      objval = []; X = [];
   end      
%%=======================================================








