%%*****************************************************************
%% compute lower and upper bounds for the exact primal
%%         optimal value. 
%%
%%  LB <= true optimal dual value = true optimal primal value <= UB.
%%
%%*****************************************************************

  function [LB,UB] = SDPvalBounds(blk,At,C,b,X,y,mu); 

  if (nargin < 7); mu = 1.1; end
  Aty = Atyfun(blk,At,[],[],y);
  Znew = ops(C,'-',Aty); 
%%
  eigX = cell(size(blk,1),1); 
  for p = 1:size(blk,1)
     pblk = blk(p,:);
     if strcmp(pblk{1},'s')
        eigX{p} = eig(full(X{p})); 
     elseif strcmp(pblk{1},'l')
        eigX{p} = X{p}; 
     end
  end
%%
%% compute lower bound
%%
  pert = 0; 
  for p = 1:size(blk,1)
     pblk = blk(p,:);
     if strcmp(pblk{1},'s')
        eigtmp = eig(full(Znew{p})); 
        idx  = find(eigtmp < 0);
        Xbar = mu*max(eigX{p}); 
     elseif strcmp(pblk{1},'l')
        eigtmp = Znew{p};
        idx  = find(eigtmp < 0);
        Xbar = mu*max(eigX{p}); 
     end      
     numneg = length(idx); 
     if (numneg) 
        mineig = min(eigtmp(idx));
        pert = pert + Xbar*sum(eigtmp(idx));
        %%fprintf('\n numneg = %3.0d,  mineigZnew = %- 3.2e',numneg,mineig);
     end
  end
  LB0 = b'*y; 
  LB  = b'*y + pert; 
  fprintf('\n <b,y> = %-10.9e, LB = %-10.9e\n',LB0,LB); 
%%
%% compute upper bound
%%
   Xbar = X;
   %% construct Xbar that is positive semidefinite 
   for p = 1:size(blk,1)
     pblk = blk(p,:);
     n = sum(pblk{2}); 
     eigXp = eigX{p}; 
     if strcmp(pblk{1},'s')
        Xbar{p} = Xbar{p} - min(eigXp)*speye(n,n);
     elseif strcmp(pblk{1},'l')
        Xbar{p} = Xbar{p} - min(eigXp)*ones(n,1);
     end
  end
  Rp  = b-AXfun(blk,At,[],Xbar);
  UB  = blktrace(blk,C,Xbar) + mu*abs(y)'*abs(Rp);
  UB0 = blktrace(blk,C,X);
  fprintf('\n <C,X> = %-10.9e, UB = %-10.9e\n',UB0,UB);
%%*****************************************************************

