%%*********************************************************************
%% norm_min: matrix 2-norm minimization problem 
%%
%% (primal) minimize || B_0 + x_1*B_1 + ... + x_m*B_m ||_2
%%
%% (dual)   maximize    - Tr B_0*Q
%%          subject to  Tr B_i*Q = 0
%%                      sum of singular values of Q is less than one.
%% 
%% The matrices B_i are of size p x q.
%%
%% The problem is equivalent to the following SDP:
%% 
%%    minimize    t 
%%                [ t*I     B(x) ] 
%%    subject to  [              ]  >= 0 
%%                [ B(x)^T  t*I  ]  
%%
%% with B(x) = B_0 + x_1*B_1 + ... + x_m*B_m.
%%
%% Adapted from Vandenberghe and Boyd, SIAM Review 96.        
%%-------------------------------------------------------------------- 
%%
%%  [blk,Avec,C,b,X0,y0,Z0,objval,x] = norm_min(B,feas,solve); 
%%
%%  Input:   B = cell array such that B{k} = B_k.
%%           feas  = 1 if want feasible starting point
%%                 = 0 if otherwise.
%%           solve = 0 just to initialize
%%                 = 1 if want to solve the problem. 
%%
%% SDPT3: version 3.0 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last modified: 2 Feb 01
%%*********************************************************************

   function [blk,Avec,C,b,X0,y0,Z0,objval,x] = norm_min(B,feas,solve);   

   if (nargin < 3); solve = 0; end; 
   if (nargin < 2); feas = 0; end; 
%%
%% dimensions 
%%
   m = size(B,2)-1; 
   [p,q] = size(B{1}); isrealB = 1;  
   for k = 1:m; 
       [p1,q1] = size(B{k}); 
       if any([p q]-[p1 q1]); 
          error('size of B_k is not the same for all k'); 
       end 
       isrealB = min(isrealB, isreal(B{k})); 
   end
   cmp = 1-isrealB; 
%%
%% A_i = [0 B_i; B_i' 0],  i=1,...,m+1;  A_{m+1} = I
%%     
   if (~cmp) 
       n = p+q;    
       blk{1,1} = 's'; blk{1,2} = n;
       b = [zeros(m,1); -1];
       C = -[zeros(p), B{1}; B{1}', zeros(q)]; 
       A = cell(1,m+1); 
       for k = 1:m
           Bk = B{k+1};  
           A{k} = [zeros(p), Bk; Bk', zeros(q)];
       end
       A{m+1} = -speye(n,n); 
       Avec = svec(blk,A,ones(size(blk,1),1)); 
   elseif (cmp)
       n = p+q;    
       blktmp{1,1} = 's'; blktmp{1,2} = n;
       b = [zeros(2*m,1); -1];
       Ctmp{1} = -[zeros(p),  B{1}; B{1}', zeros(q)]; 
       A = cell(1,2*m+1); ii = sqrt(-1); 
       for k = 1:m
           Bk = B{k+1};
           A{k}   = [zeros(p), Bk; Bk', zeros(q)];
           A{k+m} = [zeros(p), ii*Bk; -ii*Bk', zeros(q)];
       end
       A{2*m+1} = -speye(n,n);
       [blk,Avec,C,b] = convertcmpsdp(blktmp,A,Ctmp,b); 
   end  
%%
%%
   if (feas == 1); 
      X0{1} = eye(n)/n;
      y0 = 1.1*ops(C,'norm')*[zeros(length(b)-1,1); 1];
      Z0 = ops(C,'-',Atyfun(blk,Avec,[],[],y0));
   elseif (feas == 0)
      [X0,y0,Z0] = infeaspt(blk,Avec,C,b); 
   end 
   if (solve) 
      [obj,X,y,Z] = sqlp(blk,Avec,C,b,[],X0,y0,Z0);
      objval = -mean(obj); 
      if ~cmp; x = y(1:m); 
      else;    x = y(1:m) +sqrt(-1)*y(m+1:2*m); 
      end;
   else
      objval = []; x = [];
   end
%%*********************************************************************
