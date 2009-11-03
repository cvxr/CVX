%%*******************************************************
%%  Doptdesign: D-optimal experiment design.
%%
%%  max  log(det(sum_{i=1}^p lambda_i v_i v_i^T)) 
%%   
%%  s.t. lambda_i >= 0, sum_{i=1}^p lambda_i = 1. 
%%
%%  V:       nxp matrix with n <= p. 
%%  lambda:  lambda_i is the fraction of the experiments 
%%           allocated to test vector v_i.
%%  S: = V*diag(\lambda)*V'.      
%%
%% SDPT3: version 3.0 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last modified: 2 Feb 01
%%******************************************************* 

   function [blk,At,C,b,OPTIONS,lambda] = Doptdesign(V,solve);

   if (nargin == 1); solve = 0; end

   [n,p] = size(V);  
   if (n > p); 
     error(' size(V,1) > size(V,2)'); 
   end
%%
%% form At, C, b
%%
   b = zeros(p,1); 

   blk{1,1} = 's'; blk{1,2} = n; 
   F = cell(1,p); 
   for k = 1:p
      F{1,k} = -V(:,k)*V(:,k)';
   end
   At(1) = svec(blk(1,:),F,1);
   C{1,1} = sparse(n,n);
 
   blk{2,1} = 'l'; blk{2,2} = p; 
   At{2,1} = -speye(p,p); 
   C{2,1} = zeros(p,1); 
 
   blk{3,1} = 'u'; blk{3,2} = 1;
   At{3,1} = ones(1,p); 
   C{3,1} = 1; 

   OPTIONS.parbarrier{1,1} = 1; 
   OPTIONS.parbarrier{2,1} = 0; 
   OPTIONS.parbarrier{3,1} = 0; 
%%
   if (solve)
      [obj,X,y,Z] = sqlp(blk,At,C,b,OPTIONS);
      lambda = y; 
   else
      lambda = [];
   end
%%*******************************************************
