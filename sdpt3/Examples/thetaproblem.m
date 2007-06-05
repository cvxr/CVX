%%*******************************************************
%% thetaproblem: Lovasz theta number. 
%%
%% (primal) min  Tr C*X
%%          s.t. X(i,j) = 0 if (i,j) is an edge of G, 
%%               Tr(X) = 1.                         
%%  b = e1, 
%%  C = -ones(n), 
%%  A1 = eye(n), Ak = ei*ej' + ej*ei', if (i,j) is an edge. 
%%-------------------------------------------------------
%%
%% [blk,Avec,C,b,X0,y0,Z0,objval,X] = thetaproblem(G,feas,solve);
%%
%% G: adjacency matrix of a graph.
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

 function [blk,Avec,C,b,X0,y0,Z0,objval,X] = thetaproblem(G,feas,solve);

    if nargin < 2; feas  = 0; end; 
    if nargin < 3; solve = 0; end; 
    if ~isreal(G); error('only real G allowed');  end; 

    n = length(G); 
    m = sum(sum(triu(G,1))) + 1; 
    e1 = [1 zeros(1,m-1)]'; 
    C{1} = -ones(n); 
    b = e1; 
    blk{1,1} = 's';  blk{1,2} = n; 
  
    A = cell(1,m); A{1} = speye(n); 
    cnt = 2; 
    for i = 1:n 
        idx = find(G(i,i+1:n)); 
        idx = idx+i;      %% adjust index.  
        for j = 1:length(idx) 
	    A{1,cnt} = sparse([i idx(j)],[idx(j) i],[1 1],n,n);    
            cnt  = cnt + 1; 
        end
    end  

    Avec = svec(blk,A,ones(size(blk,1),1)); 
    if (feas == 1)
       y0 = -2*n*e1;
       Z0 = 2*n*eye(n)+C{1}; 
       X0 = eye(n)/n;
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

