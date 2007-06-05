%%***********************************************************
%% lmiexamp2: generate SDP data for the following LMI problem
%%
%%  max  -Tr(P)
%%  s.t. A1*P + P*A1' + B*diag(d)*B' <= 0
%%       A2*P + P*A2' + B*diag(d)*B' <= 0
%%                                -d <= 0 
%%                            sum(d)  = 1
%%***********************************************************
%% Here is an example on how to use this function to 
%% find an optimal P. 
%%
%% A1 = [-1  0  0; 0 -2  0; 1  1 -1];
%% A2 = A1 + 0.1*A1'; 
%% B  = [1 2; 3 4; 5 6]; 
%%
%% [blk,Avec,C,b] = lmiexamp2(A1,A2,B);
%% [obj,X,y,Z] = sqlp(blk,Avec,C,b);
%% n = size(A1,2); N = n*(n+1)/2; dlen = size(B,2); 
%% P = smat(blk(1,:),y(1:N)); 
%% d = y(N+[1:dlen]); 
%%
%% SDPT3: version 3.0 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last modified: 20 Apr 02
%%***********************************************************

   function [blk,Avec,C,b] = lmiexamp2(A1,A2,B); 

%%
   n1 = size(A1,2);  n2 = size(A2,2);  
   if (n1 ~= n2); error('lmiexamp2: A1, A2 not compatible'); end; 
%%   
   n = n1; 
   I = speye(n); 
   blk{1,1}  = 's'; blk{1,2} = n;
   Avec{1,1} = lmifun(A1,I,B);
   C{1,1} = sparse(n,n); 
%%
   blk{2,1}  = 's'; blk{2,2} = n;
   Avec{2,1} = lmifun(A2,I,B);
   C{2,1} = sparse(n,n); 
%%
%% and constraints:  -d <= 0
%%
   N = n*(n+1)/2; dlen = size(B,2); 
   blk{3,1}  = 'l'; blk{3,2} = dlen;
   Avec{3,1} = [sparse(dlen,N)    -speye(dlen,dlen)]; 
   C{3,1} = zeros(dlen,1);  
%%
%% add in the constraint:  sum(d) = 1
%%
   blk{4,1}  = 'u'; blk{4,2} = 1;
   Avec{4,1} = sparse([zeros(1,N) ones(1,dlen)]); 
   C{4,1} = 1; 
%%
   blktmp{1,1} = 's'; blktmp{1,2} = n;   
   b = [-svec(blktmp,I); zeros(dlen,1)]; 
%%**********************************************************

