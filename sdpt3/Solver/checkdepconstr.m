%%*****************************************************************************
%% checkdepconst: compute AAt to determine if the 
%%             constraint matrices Ak are linearly independent. 
%%              
%% [At,b,y,idxB,neardepconstr,feasible,AAt] = checkdepconstr(blk,At,b,y,rmdepconstr);
%%
%% rmdepconstr = 1, if want to remove dependent columns in At
%%             = 0, otherwise.
%% 
%% idxB = indices of linearly independent columns of original At.
%% neardepconstr = 1 if there is nearly dependent columns in At
%%            = 0, otherwise.
%% Note: the definition of "nearly dependent" is dependent on the 
%%       threshold used to determine the small diagonal elements in 
%%       the LDLt factorization of A*At. 
%% 
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%*****************************************************************************

   function  [At,b,y,idxB,neardepconstr,feasible,AAt] = ...
              checkdepconstr(blk,At,b,y,rmdepconstr);
   
   global existlowrank printlevel
   global Lsymb nnzmatold existspcholsymb
%%  
%% compute AAt
%%
   m = length(b);  
   AAt = sparse(m,m);  numdencol = 0; UU = []; 
   for p = 1:size(blk,1)
      pblk = blk(p,:); 
      if strcmp(pblk{1},'s'); 
         m1 = size(At{p,1},2); m2 = m - m1;
         if (m2 > 0) 
            if (m2 > 0.3*m); AAt = full(AAt); end
            dd = At{p,3}; 
            lenn = size(dd,1);
            idxD = [0; find(diff(dd(:,1))); lenn];
            ss = [0, cumsum(pblk{3})]; 
            if (existlowrank)
               AAt(1:m1,1:m1) = AAt(1:m1,1:m1) + At{p,1}'*At{p,1};  
               for k = 1:m2
                  idx = [ss(k)+1 : ss(k+1)];
                  idx2 = [idxD(k)+1: idxD(k+1)];
                  ii = dd(idx2,2)-ss(k); %% undo cumulative indexing
                  jj = dd(idx2,3)-ss(k);
                  len = pblk{3}(k); 
                  Dk = spconvert([ii,jj,dd(idx2,4); len,len,0]);
                  tmp = svec(pblk,At{p,2}(:,idx)*Dk*At{p,2}(:,idx)'); 
                  tmp2 = AAt(1:m1,m1+k) + (tmp'*At{p,1})'; 
                  AAt(1:m1,m1+k) = tmp2;
                  AAt(m1+k,1:m1) = tmp2';
               end
            end
            DD = spconvert([dd(:,2:4); sum(pblk{3}),sum(pblk{3}),0]);
            VtVD = (At{p,2}'*At{p,2})*DD; 
            VtVD2 = VtVD'.*VtVD; 
            for k = 1:m2
               idx0 = [ss(k)+1 : ss(k+1)]; 
               %%tmp = VtVD(idx0,:)' .* VtVD(:,idx0);
               tmp = VtVD2(:,idx0);
               tmp = tmp*ones(length(idx0),1); 
               tmp3 = AAt(m1+[1:m2],m1+k) + mexqops(pblk{3},tmp,ones(length(tmp),1),1);
               AAt(m1+[1:m2],m1+k) = tmp3; 
            end
         else
            AAt = AAt + abs(At{p,1})'*abs(At{p,1}); 
         end
      else 
         decolidx = checkdense(At{p,1}'); 
         if ~isempty(decolidx); 
            n2 = size(At{p,1},1); 
            dd = ones(n2,1); 
            len= length(decolidx); 
            dd(decolidx) = zeros(len,1);  
            AAt = AAt + (abs(At{p,1})' *spdiags(dd,0,n2,n2)) *abs(At{p,1});
            tmp = At{p,1}(decolidx,:)'; 
            UU = [UU, tmp]; 
            numdencol = numdencol + len; 
         else
            AAt = AAt + abs(At{p,1})'*abs(At{p,1});  
         end
      end
   end
   if (numdencol > 0) & (printlevel)
      fprintf('\n number of dense column in A = %d',numdencol); 
   end
   numdencol = size(UU,2); 
%%
%% 
%%
   feasible = 1; neardepconstr = 0; 
   if ~issparse(AAt); AAt = sparse(AAt); end 
   nnzmatold = mexnnz(AAt);
   rho = 1e-15;
   diagAAt = diag(AAt);  
   mexschurfun(AAt,rho*max(diagAAt,1));
   [L.R,indef,L.perm] = chol(AAt,'vector'); 
   L.d = full(diag(L.R)).^2; 
   if (indef) 
      msg = 'AAt is not pos. def.'; 
      idxB = [1:m]; 
      neardepconstr = 1; 
      if (printlevel); fprintf('\n checkdepconstr: %s',msg); end
      return; 
   end
%%
%% find independent rows of A
%%
   dd = zeros(m,1); 
   idxB = [1:m]';
   dd(L.perm) = abs(L.d); 
   idxN = find(dd < 1e-13*mean(L.d));
   ddB = dd(setdiff([1:m],idxN));
   ddN = dd(idxN);
   if ~isempty(ddN) & ~isempty(ddB) & (min(ddB)/max(ddN) < 10) 
      %% no clear separation of elements in dd
      %% do not label constraints as dependent
      idxN = []; 
   end
   if ~isempty(idxN)     
      neardepconstr = 1; 
      if (printlevel)
         fprintf('\n number of nearly dependent constraints = %1.0d',length(idxN)); 
      end
      if (numdencol==0)
         if (rmdepconstr)
            idxB = setdiff([1:m]',idxN);
            if (printlevel)
               fprintf('\n checkdepconstr: removing dependent constraints...');
            end
            [W,resnorm] = findcoeffsub(blk,At,idxB,idxN);
   	    tol = 1e-8;
            if (resnorm > sqrt(tol))
               idxB = [1:m]'; 
               neardepconstr = 0; 
               if (printlevel)
                  fprintf('\n checkdepconstr: basis rows cannot be reliably identified,'); 
                  fprintf(' abort removing nearly dependent constraints'); 
               end
               return; 
            end
            tmp = W'*b(idxB) - b(idxN);
            nnorm = norm(tmp)/max(1,norm(b)); 
            if (nnorm > tol) 
               feasible = 0; 
               if (printlevel)
                  fprintf('\n checkdepconstr: inconsistent constraints exist,');
                  fprintf(' problem is infeasible.');
               end
            else
               feasible = 1; 
               for p = 1:size(blk,1) 
                  At{p,1} = At{p,1}(:,idxB);
               end
               b = b(idxB);
               y = y(idxB); 
               AAt = AAt(idxB,idxB);               
            end
	 else
            if (printlevel)
               fprintf('\n To remove these constraints,');
               fprintf(' re-run sqlp.m with OPTIONS.rmdepconstr = 1.'); 
            end
         end
      else
         if (printlevel)
            fprintf('\n warning: the sparse part of AAt may be nearly singular.');
         end
      end
   end
%%*****************************************************************************
%% findcoeffsub: 
%%
%% [W,resnorm] = findcoeffsub(blk,At,idXB,idXN);
%% 
%% idXB = indices of independent columns of At. 
%% idxN = indices of   dependent columns of At.
%% 
%% AB = At(:,idxB); AN = At(:,idxN) = AB*W
%%
%% SDPT3: version 3.0
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last modified: 2 Feb 01
%%*****************************************************************************

   function [W,resnorm] = findcoeffsub(blk,At,idxB,idxN);

   AB = []; AN = [];
   for p = 1:size(blk,1) 
      AB = [AB; At{p,1}(:,idxB)];
      AN = [AN; At{p,1}(:,idxN)];
   end
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
%%*****************************************************************************
