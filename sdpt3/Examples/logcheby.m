%%***************************************************************
%% logcheby: logarithmic Chebyshev approximation problem. 
%%           SDP formulation. 
%%    minimize   t 
%%      x,t
%%    such that  1/t <= (x'*B(i,:))/f(i) <= t,  i = 1:p.
%% 
%%    B = pxm matrix,  f = p-vector,  p > m.
%%
%%    Ref: Vanderberghe & Boyd. 
%%--------------------------------------------------------------
%% [blk,Avec,C,b,X0,y0,Z0,obj,x] = logcheby(B,f,feas,solve);
%%
%% B = pxm matrix  (p > m).
%% f = px1 vector,  [B(:,j)./f must be positive for each j]  
%% feas  = 1 if want feasible starting point
%%       = 0 if otherwise.
%% solve = 0 if just want initialization. 
%%       = 1 if want to solve the problem. 
%%
%% SDPT3: version 3.0 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last modified: 2 Feb 01
%%***************************************************************

  function [blk,Avec,C,b,X0,y0,Z0,objval,x] = logcheby(B,f,feas,solve);
 
     if (nargin < 4); solve = 0; end; 
     if (nargin < 3); feas = 1;  end; 
     if isempty(feas); feas = 1; end; 

     if ~isreal(B); error('B must be real'); end;

     [p,m] = size(B); 
     blk{1,1} = 's'; blk{1,2} = 2*ones(1,p); 
     blk{2,1} = 'l'; blk{2,2} = p; 

     E = zeros(p,m); 
     for j = [1:m];  E(:,j) = B(:,j)./f;  end;
     if any(E < 1e-10);
        error(' B(:,j)./f must have all entry positive');    
     end;
     for i = [1:p];  beta(i) = sum(E(i,:)); end; 
%%
     temp = zeros(2*p+1,1);
     temp(2:2:2*p) = ones(p,1); 
     C{1,1} = spdiags(temp,1,2*p,2*p); C{1,1} = C{1,1} + C{1,1}'; 
     C{2,1} = zeros(p,1); 
     b  = [zeros(m,1); 1];

     A = cell(2,m+1); 
     temp = zeros(2*p,1);
     for k = 1:m
         temp(1:2:2*p-1) = -E(:,k);
         A{1,k} = spdiags(temp,0,2*p,2*p); 
         A{2,k} = E(:,k);
     end;
     temp = zeros(2*p,1);
     temp(2:2:2*p) = ones(p,1); 
     A{1,m+1} = spdiags(temp,0,2*p,2*p); 
     A{2,m+1} = ones(p,1);
%%
    Avec = svec(blk,A,ones(size(blk,1),1)); 
    if (feas == 1); 
        X0{1,1} = speye(2*p)/(2*p);
        X0{2,1} = ones(p,1)/(2*p); 
        y0 = [ones(m,1); -1.1*max(beta)]/min(beta);
        Z0 = ops(C,'-',Atyfun(blk,Avec,[],[],y0));
    elseif (feas == 0);
        [X0,y0,Z0] = infeaspt(blk,Avec,C,b);
    end;  
    if (solve)
       [obj,X,y,Z] = sqlp(blk,Avec,C,b,[],X0,y0,Z0);
       objval = -mean(obj);
       x = y(1:m);
    else
       objval = []; x = [];
    end 
%%***************************************************************
