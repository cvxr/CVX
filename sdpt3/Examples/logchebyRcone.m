%%***************************************************************
%% logchebyRcone: logarithmic Chebyshev approximation problem. 
%%                Rotated cone formulation. 
%%    minimize   t 
%%      x,t
%%    such that  1/t <= (B(i,:)*x)/f(i) <= t,  i = 1:p.
%% 
%%    B = pxm matrix,  f = p-vector,  p > m.
%%
%%    e.g. p = 20; m = 5; B = rand(p,m); f = rand(p,1); 
%%--------------------------------------------------------------
%% [blk,At,C,b,objval,x] = logchebyRcone(B,f,solve);
%%
%% B = pxm matrix  (p > m).
%% f = px1 vector,  [B(:,j)./f must be positive for each j]  
%% solve = 0 if just want initialization. 
%%       = 1 if want to solve the problem. 
%%
%% SDPT3: version 3.0 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last modified: 2 Feb 01
%%***************************************************************

  function [blk,At,C,b,objval,x] = logchebyRcone(B,f,solve);
 
     if (nargin < 3); solve = 0; end; 
     if ~isreal(B); error('B must be real'); end;

     [p,m] = size(B); 
     E = zeros(p,m); 
     for j = [1:m];  E(:,j) = B(:,j)./f;  end;
     if any(E < 1e-10);
        error(' B(:,j)./f must have all entry positive');    
     end
%%
     blk{1,1} = 'l'; blk{1,2} = p;
     blk{2,1} = 'l'; blk{2,2} = p;
     blk{3,1} = 'r'; blk{3,2} = 3;  
     blk{4,1} = 'u'; blk{4,2} = 1; 
     At{1,1} = [B, -f, zeros(p,2)]; 
     At{2,1} = [-B, zeros(p,1), f, zeros(p,1)]; 
     At{3,1} = [sparse(3,m), -speye(3,3)];
     At{4,1} = [sparse(1,m+2), 1]; 
     C{1,1} = zeros(p,1); 
     C{2,1} = zeros(p,1);
     C{3,1} = zeros(3,1);
     C{4,1} = sqrt(2); 
     b = [zeros(m,1); -1; 0; 0]; 
%%
%% convert rotated cone to socp cone
%%
    [blk,At,C,b] = convertRcone(blk,At,C,b);
%%
    if (solve)
       [obj,X,y,Z] = sqlp(blk,At,C,b);
       objval = -mean(obj);
       x = y(1:m);
    else
       objval = []; x = [];
    end 
%%***************************************************************
