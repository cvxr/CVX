%%**********************************************************
%% chebyinf.
%%
%%    minimize || p(d) ||_infty 
%%    p = monic polynomial of degree m.
%% 
%%    Here d = n-vector of complex numbers 
%%---------------------------------------------------------
%% [blk,Avec,C,b,X0,y0,Z0,objval,p] = chebyinf(d,m,solve);
%%
%% d = a vector. 
%% m = degree of polynomial. 
%% feas  = 1 if want feasible starting point
%%       = 0 if otherwise.
%% solve = 0 if just want initialization
%%       = 1 if want to solve using sdp.m
%%
%% SDPT3: version 3.0 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last modified: 2 Feb 01
%%**********************************************************
 
  function [blk,Avec,C,b,X0,y0,Z0,objval,p] = chebyinf(d,m,solve);

  if nargin <= 2; solve = 0; end;
  if (size(d,1) < size(d,2)); d = d.'; end;
  cmp = 1-isreal(d); 

  tstart=cputime; 
  n = length(d);
  e = ones(n,1);
  V(1:n,1) = e/norm(e); R(1,1) = 1/norm(e);
  for i =1:m    
      v = d.*V(:,i);  
      for j = 1:i                 %% Arnoldi iterations:
          H(j,i) = (V(:,j))'*v;   %% constructing upper-Hessenberg matrix.
          v = v - H(j,i)*V(:,j);  %% orthonormaliztion of Krylov basis.
      end;
      H(i+1,i) = norm(v);
      V(:,i+1) = v/H(i+1,i);
      R(1:i+1,i+1) = (1/H(i+1,i))*([0; R(1:i,i)] - [R(1:i,1:i)*H(1:i,i); 0]);
  end 
  if (cmp)
     blk{1,1} = 'q'; blk{1,2} = 3*ones(1,n);
     C = zeros(3*n,1); C(2:3:3*n) = real(V(:,m+1)); C(3:3:3*n) = imag(V(:,m+1));
     b = [zeros(2*m,1); -1];
     Atmp = [];
     II = [0:3:3*n-3]'; ee = ones(n,1); 
     for k=1:m
         Vk = V(:,k);
         Atmp = [Atmp; [2+II, k*ee, real(Vk)]; [3+II, k*ee, imag(Vk)]]; 
         Atmp = [Atmp; [2+II, (m+k)*ee, -imag(Vk)]; [3+II, (m+k)*ee, real(Vk)]]; 
     end
     Atmp = [Atmp;  [1+II, (2*m+1)*ee, -ones(n,1)]]; 
  else
     blk{1,1} = 'l'; blk{1,2} = 2*ones(1,n);
     b = [zeros(m,1); -1];
     C = [V(:,m+1); -V(:,m+1)];
     Atmp = [];
     II = [1:n]'; ee = ones(n,1); 
     for k=1:m
         Vk = V(:,k);
         Atmp = [Atmp; [II, k*ee, Vk]; [n+II, k*ee, -Vk]]; 
     end 
     Atmp = [Atmp; [II, (m+1)*ee, -ee]; [n+II, (m+1)*ee, -ee]];
  end
  Avec = spconvert(Atmp);
  [X0,y0,Z0] = infeaspt(blk,Avec,C,b); 
%%
  if (solve)
     [obj,X,y,Z] = sqlp(blk,Avec,C,b,[],X0,y0,Z0);
     if (cmp)
        y  = y(1:m) + sqrt(-1)*y(m+1:2*m); 
     else
        y = y(1:m); 
     end
     h  = R(1:m+1,1:m+1) \ [zeros(m,1); 1]; 
     x1 = R(1:m,1:m)*(h(m+1)*y(1:m) + h(1:m));
     p = [1; -x1(m:-1:1)];  
     objval  = max(abs(polyval(p,d)));
  else
     objval = []; p = [];
  end 
%%**********************************************************
