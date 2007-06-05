%%***************************************************
%% min sum_k bk*yk
%% s.t. sum yk*Hk  <= 0  
%%      y1 = 1
%% Hk = -hankel(ek) if   1 <= k <= n
%%    = -hankel(0,e(k-n+1)) if n+1 <= k <=2*n-1 
%%
%% [blk,At,C,b] = sdphankel(n);
%%***************************************************

   function [blk,At,C,b] = sdphankel(n);

   randn('seed',0); 
   tmp = randn(n,n); 
   tmp = tmp+tmp'; 
   X{1} = tmp + norm(tmp,'fro')*speye(n,n);
%%   
   for k = 1:n
      ek = zeros(n,1); ek(k) = -1;    
      AA{k} = sparse(hankel(ek));
   end
   zz = zeros(n,1);
   for k = n+1:2*n-1
      ek = zeros(n,1); ek(k-n+1) = -1;
      AA{k} = sparse(hankel(zz,ek)); 
   end
   blk{1,1} = 's'; blk{1,2} = n; 
   At = svec(blk,AA,1);  
   C{1} = spconvert([n n 0]);
   b = AXfun(blk,At,[],X); 
%%
   blk{2,1} = 'u'; blk{2,2} = 1; 
   ee = zeros(1,2*n-1); ee(1) = 1;
   At{2,1} = ee;
   C{2,1} = 1; 
%%***************************************************
