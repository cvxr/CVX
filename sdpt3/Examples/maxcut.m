%%*******************************************************
%% maxcut: MAXCUT problem.
%%
%% (primal problem) min  Tr C*X
%%                  s.t.  diag(X) = b; 
%%                              
%% Here, b = e,  C = -(diag(B*e)-B)/4. 
%%
%% (dual problem)   max  b'*y
%%                  s.t. diag(y) + Z = C. 
%%-------------------------------------------------------
%% [blk,Avec,C,b,X0,y0,Z0,objval,X] = maxcut(B,feas,solve);
%%
%% B: weighted adjacency matrix of a graph.
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

 function [blk,Avec,C,b,X0,y0,Z0,objval,X] = maxcut(B,feas,solve);

   if nargin < 2; feas = 0; end; 
   if nargin < 3; solve = 0; end; 
   if ~isreal(B); error('only real B allowed'); end; 
 
   n = length(B); e = ones(n,1); 
   C{1} = -(spdiags(B*e,0,n,n)-B)/4; 
   b = e;
   blk{1,1} = 's';  blk{1,2} = n;
  
   A = cell(1,n);
   for k = 1:n; A{k} = sparse(k,k,1,n,n); end; 

   Avec = svec(blk,A,ones(size(blk,1),1)); 
   if (feas)
      y0 = -1.1*abs(C{1})*e;
      Z0{1} = C{1} - spdiags(y0,0,n,n);
      X0{1} = spdiags(b,0,n,n); 
   else
      [X0,y0,Z0] = infeaspt(blk,Avec,C,b); 
   end 
   if (solve)
      [obj,X,y,Z] = sqlp(blk,Avec,C,b,[],X0,y0,Z0);
      objval = obj(1); 
   else
      objval = []; X = [];
   end      
%%*******************************************************
