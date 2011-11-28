%%***********************************************************
%% etp: Education testing problem.
%%
%% (dual problem) maximize     e'*d
%%                  subject to   B - diag(d) >= 0
%%                               d >= 0
%%
%% (primal problem) minimize     Tr B*X
%%                  subject to   X >= 0
%%                               diag(X) >= e
%%
%% Ref: M.T. Chu, J.W. Wright, IMA J. of Numerical Anal.,
%%      15 (1995), pp. 141--160.
%%-----------------------------------------------------------
%% [blk,Avec,C,b,X0,y0,Z0,objval,d] = etp(B,feas,solve);
%%
%% B = nxn positive definite.
%% feas  = 1 if want feasible starting point
%%       = 0 if otherwise.
%% solve = 0 just to initialize
%%       = 1 if want to solve the problem.
%%
%% SDPT3: version 3.0 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last modified: 2 Feb 01
%%***********************************************************

   function [blk,Avec,C,b,X0,y0,Z0,objval,d] = etp(B,feas,solve);

   if (nargin < 2); feas = 0; end
   if (nargin < 3); solve = 0; end 
   if isempty(feas); feas = 0; end
   if (~isreal(B))
      error('only real B allowed');
   elseif (norm(B-B','fro') > 1e-13);
      error(' B must be symmetric'); 
   end;
%%
%% validate B
%%
   n = length(B);  
   d = eig(B); d = real(d);
   if (min(d) < 0); 
      error('B must be positive def'); 
   end;
%%
%%
   blk{1,1} = 's'; blk{1,2} = n; 
   blk{2,1} = 'l'; blk{2,2} = n;  
   b = ones(n,1);  
   C{1,1} = B; 
   C{2,1} = zeros(n,1); 

   A = cell(2,n); 
   for k = 1:n 
       A{1,k} = sparse(k,k,1,n,n); 
       A{2,k} = [zeros(k-1,1); -1; zeros(n-k,1)]; 
   end;  

   Avec = svec(blk,A,ones(size(blk,1),1)); 
   if (feas == 1); 
      y0 = 0.9*min(d)*ones(n,1);
      Z0 = ops(C,'-',Atyfun(blk,Avec,[],[],y0));      
      X0{1,1} = 1.1*eye(n); 
      X0{2,1} = 0.1*ones(n,1); 
   elseif (feas == 0);   
      [X0,y0,Z0] = infeaspt(blk,Avec,C,b); 
   end;
   if (solve)
      [obj,X,y,Z] = sqlp(blk,Avec,C,b,[],X0,y0,Z0);
      objval = obj(2); 
      d = y;
   else
      objval = []; d = [];
   end
%%===========================================================

