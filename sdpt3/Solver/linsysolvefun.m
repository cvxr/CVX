%%*************************************************************************
%% linsysolvefun: Solve H*x = b
%%
%% x = linsysolvefun(L,b)
%% where L contains the triangular factors of H. 
%%
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%*************************************************************************
 
  function x = linsysolvefun(L,b)
 
  x = zeros(size(b)); 
  for k=1:size(b,2)
     if strcmp(L.matfct_options,'chol')
        x(L.perm,k) = mextriang(L.R, mextriang(L.R,b(L.perm,k),2) ,1); 
        %% x(L.perm,k) = L.R \ (b(L.perm,k)' / L.R)';
     elseif strcmp(L.matfct_options,'spchol')
        x(L.perm,k) = mextriangsp(L.Rt,mextriangsp(L.R,b(L.perm,k),2),1);
     elseif strcmp(L.matfct_options,'ldl')
	x(L.p,k) = ((L.D\ (L.L \ b(L.p,k)))' / L.L)';
     elseif strcmp(L.matfct_options,'spldl')
        btmp = b(:,k).*L.s;
        xtmp(L.p,1) = L.Lt\ (L.D\ (L.L \ btmp(L.p)));
	x(:,k) = xtmp.*L.s; 
     elseif strcmp(L.matfct_options,'lu')
        x(:,k) = L.U \ (L.L \ b(L.p,k));
     elseif strcmp(L.matfct_options,'splu') 
	btmp = b(:,k)./L.s; 
        x(L.q,k) = L.U \ (L.L \ (btmp(L.p)));
     end
  end
%%*************************************************************************
