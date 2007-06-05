%%*****************************************************************************
%% findepcols: find dependent columns of a matrix
%%              
%%  [idxB,W] = findepcols(A);
%%
%% AB = A(:,idxB), where idxB are the indicies of independent columns. 
%% AN = A(:,idxN) = AB*W, where idxN are the indices of dependent columns. 
%%
%% SDPT3: version 3.0 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 7 May 01
%%*****************************************************************************

   function  [idxB,W,dd] = findepcols(A);
%%
%% find independent columns of A
%%
      [m,n2] = size(A);
      idxB = [1:n2]'; W = []; dd = [];
      numdenrow = 0;
      BB = A'*A; 
      L = sparcholfun(symbcholfun(BB),BB);
      dd(L.perm) = abs(L.d);
      tol = 1e-13;  
      idxN = find(dd < tol);
      idxB = setdiff([1:n2]',idxN);
      if ~isempty(idxN)
         A2 = A(:,idxB);
         BB2 = A2'*A2; 
         L = sparcholfun(symbcholfun(BB2),BB2);
         dd2(L.perm) = abs(L.d); 
         idxN2 = find(dd2 < tol);
         idxN = [idxN idxB(idxN2)];
         %%
         fprintf('\n number of nearly dependent columns = %1.0d',length(idxN)); 
         idxB = setdiff([1:n2]',idxN);
         AB = A(:,idxB); AN = A(:,idxN);     
         if (nargout >= 2) 
            W = findcoeff(AB,AN);
         end
      else
         AB = A; AN = []; 
      end
%%
%%*****************************************************************************
%% findcoeff: find W so that AN = AB*W. 
%%
%% W = findcoeff(AB,AN);
%% 
%% SDPT3: version 3.0
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last modified: 2 Feb 01
%%*****************************************************************************

   function W = findcoeff(AB,AN);

   tol = 1e-8; 
   [m,n] = size(AB); 
%%
   matlabversion = sscanf(version,'%f');
   matlabversion = matlabversion(1);
   if (matlabversion >= 6.5)
      options = 2; 
   else
      options = 1;
   end
   if (options == 1)
      tmp = AB'*AB;
      S = sparcholfun(symbcholfun(tmp),tmp); 
      rhs = AB'*AN;
      for k=1:size(AN,2)
         W(:,k) = bwblkslvfun(S, fwblkslvfun(S,rhs(:,k)) ./ S.d);
      end
   else 
      [L,U,P,Q] = lu(sparse(AB));    
      rhs  = P*AN;
      Lhat = L(1:n,:); 
      W = Q*( U \ (Lhat \ rhs(1:n,:))); 
   end
   nnorm = norm(AN-AB*W,'fro')/max(1,norm(AN,'fro'));
   if (nnorm > tol) 
      fprintf('\n basis rows may be identified incorrectly, %3.1e',nnorm); 
   end
%%*****************************************************************************
