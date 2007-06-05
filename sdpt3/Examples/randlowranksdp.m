%%******************************************************************
%% randlowranksdp.m : creates random feasible SDP problems where the 
%%                constraint matrices are low-rank matrices of the 
%%                form V*diag(d)*V'. 
%%
%% [blk2,At2,C,b,blk,At] = randlowranksdp(n,m1,m2,r);
%%
%% blk2,At2: data with low-rank structure coded.
%% blk, At:  data without taking low-rank structure into account.
%%
%% n = size of the sdp variable
%% m1 = number of general constraints
%% m2 = number of low-rank constraints
%% r = rank of each constraint matrix. 
%%
%%******************************************************************

  function [blk2,At2,C,b,blk,At] = randlowranksdp(n,m1,m2,r);

  blk = cell(1,2);
  if (m1 > 0) 
     [blk,At0,C,b0] = randsdp(n,[],[],m1);
  else 
     b0 = []; At0 = cell(1);
  end
  if (m2 > 0)
     if (nargout > 4)
        [blk2,At2,C,b2,blk,At1] = randlowranksdpfun(n,m2,r);
     else
        [blk2,At2,C,b2] = randlowranksdpfun(n,m2,r);
     end
  else 
     b2 = []; At1 = cell(1); blk2 = blk; 
  end
  b = [b0; b2];
  if (nargout > 4)
     At{1} = [At0{1}, At1{1}];
  end
  At2{1,1} = At0{1}; 
%%******************************************************************
%%******************************************************************
  function [blk2,At2,C,b,blk,At] = randlowranksdpfun(n,m,r);

  randn('state',0);

  blk{1,1} = 's'; blk{1,2} = n; 
  blk2{1,1} = 's'; blk2{1,2} = n; blk2{1,3} = r*ones(1,m); 
  %%
  %% construct data low rank structure
  %%
  At2 = cell(1,3); 
  b = zeros(m,1); 
  X0 = randn(n); X0 = X0*X0'; X0 = 0.5*(X0+X0');
  ss = [0,cumsum(blk2{1,3})];
  V = randn(n,m*r);  
  dd = [];
  for k = 1:length(blk2{1,3})
     idx = [ss(k)+1 : ss(k+1)];
     len = blk2{1,3}(k);
     Dk = randn(len,len); Dk = 0.5*(Dk+Dk');
     [ii,jj,vv] = find(Dk);
     numnz = length(ii);
     dd = [dd; k*ones(numnz,1),ii,jj,vv]; %% each row has the form [constr,i,j,val]
     tmp1 = X0*V(:,idx); 
     tmp2 = V(:,idx)*Dk;
     b(k) = sum(sum(tmp1.*tmp2));
  end
  At2{1,2} = V; At2{1,3} = dd;
  Z0 = X0; 
  y0 = randn(m,1);
  Aty = At2{1,2}*spdiags(mexexpand(blk2{1,3},y0),0,r*m,r*m)*At2{1,2}';
  C{1} = Z0 + norm(Aty,'fro')*speye(n,n) + Aty; 
  C{1} = 0.5*(C{1}+C{1}');
  %%
  %% construct data without exploiting low rank structure
  %%
  if (nargout > 4)
     idxD = [0; find(diff(dd(:,1))); size(dd,1)];
     for k = 1:m
        idx = [ss(k)+1 : ss(k+1)];
        Vk = At2{1,2}(:,idx);
        len = blk2{1,3}(k);
        idx2 = [idxD(k)+1:idxD(k+1)];
        Dk = spconvert([dd(idx2,2:4); len,len,0]);
        A{k} = Vk*Dk*Vk';
     end
     At = svec(blk,A);
  end
%%******************************************************************
