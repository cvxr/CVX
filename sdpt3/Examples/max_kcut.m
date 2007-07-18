%%*******************************************************
%% max_kcut: MAX k-CUT problem.
%%
%% (primal problem) min  <C,X>
%%                  s.t. diag(X) = b
%%                       Xij >= -1/(K-1) for all i ~= j
%%                       X positive semidefinite
%%
%% Here, b = e,  C = -(1-1/K)/2* (diag(B*e)-B). 
%%-------------------------------------------------------
%% [blk,Avec,C,b,objval,X] = max_kcut(B,K,solve);
%%
%% B: weighted adjacency matrix of a graph.
%% solve = 0 just to initialize
%%       = 1 if want to solve the problem. 
%%
%% See graph.m --- generate random adjacency matrix. 
%%
%% SDPT3: version 4.0 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last modified: 18 May 07
%%*******************************************************

 function [blk,Avec,C,b,objval,X] = max_kcut(B,K,solve);

   if (nargin < 3); solve = 0; end; 
   if ~isreal(B); error('only real B allowed'); end; 
 
   n = length(B); e = ones(n,1); 
   n2 = n*(n-1)/2; 
   C{1} = -(1-1/K)/2*(spdiags(B*e,0,n,n)-B); 
   b = e;
   blk{1,1} = 's';  blk{1,2} = n;
   blk{2,1} = 'l';  blk{2,2} = n2; 
  
   A = cell(1,n);
   for j = 1:n; A{j} = sparse(j,j,1,n,n); end; 
   Avec = svec(blk(1,:),A,1); 
   tmp = speye(n*(n+1)/2); 
   idx = cumsum([1:n]); 
   Atmp = tmp(:,setdiff([1:n*(n+1)/2],idx)); 
   Avec{1,1} = [Avec{1,1}, Atmp/sqrt(2)]; 
   Avec{2,1} = [sparse(n2,n), -speye(n2,n2)]; 
   b = [b; -1/(K-1)*ones(n2,1)]; 
   C{2,1} = zeros(n2,1); 
%%
   if (solve)
      [obj,X,y,Z] = sqlp(blk,Avec,C,b);
      objval = obj(1); 
   else
      objval = []; X = [];
   end      
%%*******************************************************
