%%*************************************************************************
%% ToeplitzApproxSQQ: find the nearest symmetric positive definite Toeplitz 
%%                    matrix to a given symmetric matrix F. 
%%
%%  max -y0
%%  s.t.   y0*B + T(y)                (S>=) 0 
%%        [y0; gam.*y] + [0; q./gam]  (Q>=) 0 
%%
%%  where B = diag([zeros(n,1); 1])
%%        q(1) = - Tr(F); q(k+1) = -sum of upper and lower kth diagonals of F
%%        gam(1) = sqrt(n); gam(k) = sqrt(2*(n-k+1)) for k=2:n 
%%*************************************************************************

 function [blk,At,C,b] =  ToeplitzApproxSQQ(F)

 n = length(F); 
 gam = sqrt([n, 2*(n-1:-1:1)]');
 q = zeros(n,1); 
 q(1) = -sum(diag(F));
 for k=1:n-1
    q(k+1) = -2*sum(diag(F,k)); 
 end
 beta = norm(F,'fro')^2;
 
 blk{1,1} = 's'; blk{1,2} = n+1;
 blk{2,1} = 'q'; blk{2,2} = n+1; 
 
 b = [-1; zeros(n,1)]; 
 C{1,1} = sparse(n+1,n+1); 
 C{2,1} = [0; q./gam]; 
 
 Acell = cell(1,n+1);
 Acell{1} = -spconvert([n+1,n+1,1]); 
 Acell{2} = -spdiags([ones(n,1); 0],0,n+1,n+1); 
 for k = 1:n-1
    tmp = -spdiags([ones(n,1); 0],k,n+1,n+1);  
    Acell{k+2} = tmp + tmp';
 end
 At(1,1) = svec(blk(1,:),Acell,1); 
 
 At{2,1} = -spdiags([1; gam],0,n+1,n+1);
%%***********************************************************************