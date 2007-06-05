%%*******************************************************************
%% schurmat_lblk: compute A*D*A'
%%
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%*******************************************************************

   function [schur,UU,EE] = schurmat_lblk(blk,At,par,schur,UU,EE,p,dd);

   global idxdenAl
   
   iter = par.iter; 
   n = sum(blk{p,2});  

   if (iter==1) 
      idxdenAl{p} = checkdense(At{p}'); 
   end
   ddsch = dd{p}; 
   if ~isempty(idxdenAl{p}); 
      idxden = idxdenAl{p}; 
      len = length(idxden); 
      Ad = At{p}(idxden,:)' *spdiags(sqrt(ddsch(idxden)),0,len,len); 
      UU = [UU, Ad];
      if isempty(EE)
         count = 0; 
      else
         count = max(max(EE(:,1)),max(EE(:,2))); 
      end
      tmp = count + [1:len]'; 
      EE = [EE; [tmp, tmp, -ones(len,1)] ]; 
      ddsch(idxden) = zeros(len,1); 
   end
   schurtmp = At{p}' *spdiags(ddsch,0,n,n) *At{p}; 
   schur = schur + schurtmp;
%%*******************************************************************
