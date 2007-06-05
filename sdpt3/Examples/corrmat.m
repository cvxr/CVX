%%****************************************************
%% corrmat: find the nearest correlation matrix
%%          to a given symmetric matrix. 
%%
%%  [blk,At,C,b] = corrmat(H); 
%%
%%  min || X - H ||_F
%%  s.t.  diag(X) = e, X psd.  
%%
%%  min   t 
%%  s.t.  diag(X) = e
%%        svec(X) + y = svec(H), 
%%        X psd, ||y||_2 <= t. 
%%
%% SDPT3: version 3.0 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last modified: 18 Dec 01
%%****************************************************

  function [blk,At,C,b] = corrmat(H); 

   if size(H,1) ~= size(H,2); error('corrmat: matrix must be square'); end;
   n = length(H); 
   n2 = n*(n+1)/2; 

   blk{1,1} = 's'; blk{1,2} = n;
   for k=1:n; AA{1,k} = spconvert([k,k,1; n,n,0]); end; 
   matrepdiag = svec(blk(1,:),AA); 
   At{1,1} = [matrepdiag{1},  speye(n2)]; 

   blk{2,1} = 'q'; blk{2,2} = n2+1;
   At{2,1} = [sparse(n,n2+1); sparse(n2,1), speye(n2)]; 
   b = [ones(n,1); svec(blk(1,:),H)];
   C{1,1} = sparse(n,n); C{2,1} = [1; zeros(n2,1)]; 
%%****************************************************
