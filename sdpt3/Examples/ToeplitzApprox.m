%%*************************************************************************
%% ToeplitzApprox: find the nearest symmetric positive definite Toeplitz 
%%                 matrix to a given symmetric matrix F. 
%%
%%
%%  max -y(n+1)
%%  s.t.  T(y(1:n)) + y(n+1)*B   >= 0 
%%        [I   0  ] + sum_{k=1}^n y(k) [0          gam(k)*e_k ]  + y(n+1)*B >= 0 
%%        [0 -beta]                    [gam(k)*e_k'  -2q(k)   ] 
%%
%%  where B = diag([zeros(n,1); 1])
%%        q(1) = - Tr(F); q(k+1) = -sum of upper and lower kth diagonals of F
%%        gam(1) = sqrt(n); gam(k) = sqrt(2*(n-k+1)) for k=2:n 
%%        beta = norm(F,'fro')^2
%%*************************************************************************

 function [blk,At,C,b] =  ToeplitzApprox(F)

 n = length(F); 
 gam = sqrt([n, 2*(n-1:-1:1)]);
 q = zeros(n,1); 
 q(1) = -sum(diag(F));
 for k=1:n-1
    q(k+1) = -2*sum(diag(F,k)); 
 end
 beta = norm(F,'fro')^2;
 
 blk{1,1} = 's'; blk{1,2} = n+1;
 blk{2,1} = 's'; blk{2,2} = n+1; 
 
 b = [zeros(n,1); -1]; 
 C{1,1} = sparse(n+1,n+1); 
 C{2,1} = spdiags([ones(n,1); -beta],0,n+1,n+1); 
 
 Acell = cell(1,n+1);
 Acell{1} = -spdiags([ones(n,1); 0],0,n+1,n+1); 
 for k = 1:n-1
    tmp = -spdiags([ones(n,1); 0],k,n+1,n+1);  
    Acell{k+1} = tmp + tmp';
 end
 Acell{n+1} = -spconvert([n+1,n+1,1]); 
 At(1,1) = svec(blk(1,:),Acell,1); 
 
 for k = 1:n
    Acell{k} = -spconvert([k, n+1, gam(k); n+1, k, gam(k); n+1, n+1, -2*q(k)]); 
 end
 Acell{n+1} = -spconvert([n+1,n+1,1]); 
 At(2,1) = svec(blk(2,:),Acell,1);  
%%***********************************************************************
