%%***********************************************************
%% lmiexamp1: generate SDP data for the following LMI problem
%%
%%  max  -eta
%%  s.t. B*P + P*B'  <= 0
%%        -P         <= -I 
%%         P - eta*I <= 0
%%         P(1,1)     = 1
%%***********************************************************
%% Here is an example on how to use this function to 
%% find an optimal P. 
%%
%% B = [-1  0  0; 5 -2  0; 1  1 -1];
%% [blk,At,C,b] = lmiexamp1(B);
%% [obj,X,y,Z] = sqlp(blk,At,C,b);
%% P = smat(blk(1,:),y);  
%%
%% SDPT3: version 3.0 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last modified: 20 Apr 02
%%***********************************************************

   function [blk,At,C,b] = lmiexamp1(B); 

   n = length(B); n2 = n*(n+1)/2; 
   I  = speye(n); 
   z0 = sparse(n2,1); 
   blktmp{1,1} = 's'; blktmp{1,2} = n;
%%
   blk{1,1} = 's'; blk{1,2} = n;
   blk{2,1} = 's'; blk{2,2} = n;
   blk{3,1} = 's'; blk{3,2} = n;
   blk{4,1} = 'u'; blk{4,2} = 1; 
%%
   At{1,1} = [lmifun(B,I),     z0];
   At{2,1} = [lmifun(-I/2,I),  z0]; 
   At{3,1} = [lmifun(I/2,I),   svec(blktmp,-I,1)]; 
   At{4,1} = sparse([1, zeros(1,n2)]); 
%%   
   C{1,1} = sparse(n,n); 
   C{2,1} = -speye(n); 
   C{3,1} = sparse(n,n); 
   C{4,1} = 1; 
%%
   b = [zeros(n2,1); -1]; 
%%**********************************************************
