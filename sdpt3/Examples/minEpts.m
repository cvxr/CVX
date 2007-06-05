%%********************************************************************
%% minEpts: minimum volume ellipsoid containing given 
%%          points, V(:,1),...,V(:,m). 
%%
%% max log(det(B))
%% s.t. || Bx + d || <= 1, for all x = V(:,1),...,V(:,m)
%%
%% [blk,At,C,b,OPTIONS,B,d] = minEpts(V,solve)
%% 
%% For problem formulation, see 
%% Vandenberghe, Boyd, Wu, Determinant maximization with linear 
%% matrix inequalities, SIAM J. Matrix Analysis and Applications,
%% 19 (1998), pp. 499--533.
%%
%%********************************************************************

  function  [blk,At,C,b,OPTIONS,B,d] = minEpts(V,solve)

  if (nargin < 2); solve = 0; end; 
%%
%% form data matrices 
%%
  [p,m] =  size(V); N = (p+1)*m;   
  blk{1,1} = 's'; blk{1,2} = (p+1)*ones(1,m);
  blk{2,1} = 's'; blk{2,2} = p; 
%%
  count = 0; 
  for j = 1:p
     s1 = V(j,:)'; i1 = [j:(p+1):(p+1)*(m-1)+j]; j1 = [p+1:p+1:(p+1)*m];
     tmp = sparse(i1,j1,s1,N,N);
     tmp = tmp + tmp';   
     count = count + 1; 
     F{1,count} = -tmp;
     F{2,count} = sparse(j,j,-1,p,p); 
  end
  for j = 2:p
     for k = 1:j-1    
        s1 = V(k,:)'; i1 = [j:(p+1):(p+1)*(m-1)+j]; j1 = [p+1:p+1:(p+1)*m];  
        s2 = V(j,:)'; i2 = [k:(p+1):(p+1)*(m-1)+k]; j2 = [p+1:p+1:(p+1)*m]; 
        tmp  = sparse(i1,j1,s1,N,N);   
        tmp  = tmp + sparse(i2,j2,s2,N,N);
        tmp = tmp + tmp'; 
        count = count + 1; 
        F{1,count} = -tmp;         
        F{2,count} = sparse([j k],[k j],[-1 -1],p,p);     
     end
  end
  for j = 1:p
     s1 = ones(m,1); 
     i1 = [j:(p+1):(p+1)*(m-1)+j]; j1 = [p+1:p+1:(p+1)*m];
     tmp = sparse(i1,j1,s1,N,N);  
     tmp = tmp + tmp'; 
     count = count + 1; 
     F{1,count} = -tmp;
     F{2,count} = sparse(p,p); 
  end        
  At = svec(blk,F,ones(2,1));
  C{1,1} = speye(N); 
  C{2,1} = zeros(p);     
  b = zeros(p*(p+1)/2+p,1);       
  parbarrier{1,1} = 0;
  parbarrier{2,1} = 1; 
  OPTIONS.parbarrier = parbarrier; 
%%
%% || Bx + d || <= 1. 
%% x'*(B'*B)*x + 2(B'*d)'*x  + d'*d <= 1.
%%
  if (solve) 
     [obj,X,y,Z,info] = sqlp(blk,At,C,b,OPTIONS);
     if (length(y) ~= p*(p+3)/2)
        error('length of y not compatible with p'); 
     end
     B = diag(y(1:p));
     tmp = p; 
     for k = 1:p-1
         B(k+1:p,k) = y(tmp + [1:p-k]);
         B(k,k+1:p) = B(k+1:p,k)';
         tmp = tmp + p-k;
     end; 
     d = y(tmp+1:length(y)); 
  else 
     B = []; d = []; 
  end 
%%********************************************************************
