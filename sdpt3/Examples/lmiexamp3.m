%%***********************************************************
%% lmiexamp3: generate SDP data for the following LMI problem
%%
%%  max  eta
%%  s.t. [A*P+P*A'  P*B'] +  eta*[0    0]  <= [-G  0]
%%       [B*P       0   ]        [0    I]     [ 0  0]
%%        P >= 0 
%%  P and eta are variables, P symmetric.
%%
%%  Ref: Body et al, Linear matrix inequalities in system and 
%%       control theory,  p. 10. 
%%***********************************************************
%% Here is an example on how to use this function to 
%% find an optimal P. 
%%
%% A = [-1  0  0; 0 -2  0; 1  1 -1];
%% B = [1  3  5; 2 4 6]; 
%% G = ones(3,3); 
%%
%% [blk,Avec,C,b] = lmiexamp3(A,B,G);
%% [obj,X,y,Z] = sqlp(blk,Avec,C,b);
%% n = size(A,2); N = n*(n+1)/2; 
%% blktmp{1,1} = 's'; blktmp{1,2} = n; 
%% P = smat(blktmp,y(1:N)); 
%%
%% SDPT3: version 3.0 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last modified: 20 Apr 02
%%***********************************************************

   function [blk,Avec,C,b] = lmiexamp3(A,B,G); 
%%
   [m,n] = size(A);  
   [m2,n2] = size(B); 
   if (n ~= n2); error('lmiexamp3: A, B not compatible'); end; 
%%  
   blk{1,1} = 's'; blk{1,2} = m + m2; 
   I = speye(n);  
   Avec(1,1) = lmifun2(A,I,I,B);  
   tmp =  [sparse(m,n+m2);  sparse(m2,n) speye(m2,m2)];
   Avec{1} = [Avec{1} svec(blk,tmp,1)]; 
%%
   C{1,1} = [-G  sparse(m,m2); sparse(m2,n+m2)];
%%
   N = n*(n+1)/2; 
   b = [zeros(N,1); 1]; 
%%**********************************************************

