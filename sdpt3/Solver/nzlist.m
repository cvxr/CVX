%%***********************************************************************
%% nzlist: find the combined list of non-zero elements
%%         of Aj, j = 1:k, for each k,
%%         assuming that the Aj's are permuted such that 
%%         A1 has the fewest nonzero elements, followed by A2, and so on.
%%
%% [isspA,nzlistA,nzlistAsum,isspAy,nzlistAy] = nzlist(blk,At,par)
%%
%%  isspA(p,k) = 1 if Apk is sparse, 0 if it is dense. 
%%  nzlistA = px2 cell array.
%%            nzlistA{p,1}(k) is the starting row index (in C convention) 
%%            in the 2-column matrix nzlistA{p,2} that
%%            stores the row and column index of the nonzero elements
%%            of Apk. 
%%            nzlistA{p,1}(k) = inf if nnz(Apk) exceeds given threshold. 
%%  nzlistAsum = px2 cell array.
%%            nzlistA{p,1}(k) is the starting row index (in C convention) 
%%            in the 2-column matrix nzlistA{p,2} that
%%            stores the row and column index of the nonzero elements
%%            of Apk that are not already present
%%            in the combined list from Ap1+...Ap,k-1.
%%  nzlistAy = px1 cell array.
%%            nzlistAy{p} is a 2-column matrix that stores the 
%%            row and column index of the nonzero elements of
%%            Ap,1+.... Ap,m. 
%%            nzlistAy{p} = inf  if the number of nonzero elements 
%%            exceeds a given threshold. 
%%
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%***********************************************************************

  function [isspA,nzlistA,nzlistAsum,isspAy,nzlistAy] = nzlist(blk,At,par)

  spdensity   = par.spdensity;
  smallblkdim = par.smallblkdim; 
  m = par.numcolAt; 
%%
  numblk = size(blk,1); 
  isspA = zeros(numblk,m);   
  nzlistA = cell(numblk,2);  nzlistAsum = cell(numblk,2);   
  isspAy = zeros(numblk,1);  nzlistAy = cell(numblk,1); 
%%
  for p = 1:size(blk,1)
     pblk = blk(p,:); 
     if strcmp(pblk{1},'s') & ((max(pblk{2}) > smallblkdim) | (length(pblk{2}) <= 10))
        numblk = length(pblk{2}); 
        n = sum(pblk{2});
        n2 = sum(pblk{2}.*pblk{2});  
        if (numblk == 1) 
           nztol = spdensity*n;
           nztol2 = spdensity*n2/2; 
        else
           nztol = spdensity*n/2; 
           nztol2 = spdensity*n2/4; 
        end
        nzlist1 = zeros(1,m+1); nzlist2 = []; 
        nzlist3 = zeros(1,m+1); nzlist4 = []; breakyes = zeros(1,2); 
        Asum = sparse(n,n); 
%%
        m1 = size(At{p,1},2); 
        for k = 1:m1
           Ak = mexsmat(blk,At,1,p,k); 
           nnzAk = nnz(Ak); 
           isspA(p,k) = (nnzAk < spdensity*n2) | (numblk > 1); 
           if ~all(breakyes)
              [I,J] = find(abs(Ak) > 0);  
              %%
              %% nonzero elements of Ak.
              %%
   	      if (breakyes(1) == 0); 
                 if (nnzAk <= nztol) 
                    idx = find(I<=J); 
                    nzlist1(k+1) = nzlist1(k)+length(idx);      
                    nzlist2 = [nzlist2; [I(idx), J(idx)] ]; 
                 else
                    nzlist1(k+1:m+1) = inf*ones(1,m-k+1);
                    breakyes(1) = 1; 
                 end
              end
              %%   
	      %% nonzero elements of |A1|+...+|Ak|.
              %%
	      if (breakyes(2) == 0)
                 nztmp = zeros(length(I),1); 
                 for t = 1:length(I); 
                    i=I(t); j=J(t); nztmp(t)=Asum(i,j);           
                 end
                 %% find new nonzero positions when Ak is added to Asum. 
                 idx = find(nztmp == 0);   
                 nzlist3(k+1) = nzlist3(k) + length(idx); 
                 if (nzlist3(k+1) < nztol2); 
                    nzlist4 = [ nzlist4; [I(idx), J(idx)] ];
                 else
                    nzlist3(k+1:m+1) = inf*ones(1,m-k+1);  
                    breakyes(2) = 1;
                 end
                 Asum = Asum+abs(Ak);
              end
           end
        end
	if (numblk == 1) 
   	   isspAy(p,1) = (nzlist1(m+1) < inf) | (nzlist3(m+1) < inf);
        else
   	   isspAy(p,1) = 1;
        end
	nzlistA{p,1} = nzlist1;    nzlistA{p,2} = nzlist2;          
    	nzlistAsum{p,1} = nzlist3; nzlistAsum{p,2} = nzlist4;         
        %%
        %% nonzero elements of (A1*y1+...Am*ym). 
        %% 
        if (nzlist3(m+1) < inf); 
           if (length(pblk) > 2) 
              m2 = length(pblk{3});
              len = sum(pblk{3}); 
              DD = spconvert([At{p,3}(:,2:4); len, len, 0]);
              Asum = Asum + abs(At{p,2}*DD*At{p,2}');
           end
           [I,J] = find(Asum > 0); 
           if (length(I) < nztol2)
              nzlistAy{p} = [I, J];
           else
              nzlistAy{p} = inf; 
           end
        else
           nzlistAy{p} = inf; 
        end
     end
  end
%%***********************************************************************
