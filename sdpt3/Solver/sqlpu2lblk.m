%%***************************************************************************
%% sqlpu2lblk: decide whether to convert ublk to lblk
%%
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 10 Jul 2007
%%***************************************************************************

   function [blk,At,C,X,Z,ublk2lblk,ublkidx] = sqlpu2lblk(blk,At,C,X,Z,par,convertlen); 

   cachesize = 128; 
%%
   ublk2lblk = zeros(size(blk,1),1);  
   ublkidx = cell(size(blk,1),2); 

   for p = 1:size(blk,1) 
      pblk = blk(p,:); n0 = sum(pblk{2}); 
      if strcmp(pblk{1},'u') & (pblk{2} > 0) 
         ublk2lblk(p) = 1; 
         if (pblk{2} > convertlen); return; end
         AAt = At{p}*At{p}';
         mexschurfun(AAt,1e-15*max(1,diag(AAt)));      
         [L.R,indef,L.perm] = chol(AAt,'vector');
         L.d = full(diag(L.R)).^2; 
         if (~indef) & (max(L.d)/min(L.d) < 1e6) 
            ublk2lblk(p) = 0; 
            msg = '*** no conversion for ublk'; 
            if (par.printlevel); fprintf(' %s',msg); end
         else
            dd(L.perm,1) = abs(L.d); 
            idxN = find(dd < 1e-11*mean(L.d));
            idxB = setdiff([1:n0]',idxN);
            ddB  = dd(idxB);
            ddN  = dd(idxN);
            if ~isempty(ddN) & ~isempty(ddB) & (min(ddB)/max(ddN) < 10) 
               idxN = []; idxB = [1:n0]'; 
            end
            ublkidx{p,1} = n0; ublkidx{p,2} = idxN; 
            if ~isempty(idxN)
   	       restol = 1e-8;
               [W,resnorm] = findcoeff(At{p}',idxB,idxN);
               resnorm(2) = norm(C{p}(idxN) - W'*C{p}(idxB));
               if (max(resnorm) < restol)
                  feasible = 1; 
                  blk{p,2} = length(idxB); 
                  Atmp  = At{p}'; 
                  At{p} = Atmp(:,idxB)'; 
                  C{p}  = C{p}(idxB);  
                  X{p} = X{p}(idxB); Z{p} = Z{p}(idxB); 
                  msg = 'removed dependent columns in constraint matrix for ublk'; 
                  if (par.printlevel); fprintf('\n %s\n',msg); end
               end 
            end
         end
      end
   end
%%***************************************************************************
%%***************************************************************************
%% findcoeff: 
%%
%% [W,resnorm] = findcoeff(A,idXB,idXN);
%% 
%% idXB = indices of independent columns of A. 
%% idxN = indices of   dependent columns of A.
%% 
%% AB = A(:,idxB); AN = A(:,idxN) = AB*W
%%
%% SDPT3: version 3.0
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last modified: 2 Feb 01
%%***************************************************************************

   function [W,resnorm] = findcoeff(A,idxB,idxN);

   AB = A(:,idxB);
   AN = A(:,idxN);
   [m,n] = size(AB); 
%%
%%-----------------------------------------
%% find W so that AN = AB*W
%%-----------------------------------------
%% 
   [L,U,P,Q] = lu(sparse(AB));    
   rhs  = P*AN;
   Lhat = L(1:n,:); 
   W = Q*( U \ (Lhat \ rhs(1:n,:))); 
   resnorm = norm(AN-AB*W,'fro')/max(1,norm(AN,'fro'));
%%***************************************************************************
