%%******************************************************
%% control: a SDP problem from control theory. 
%%
%%    max  t
%%   s.t   Bk'*P + P*Bk <= 0, k = 1:L
%%            P <= I,   
%%         P >= t*I,  P = P'.
%%
%%------------------------------------------------------
%%
%%  [blk,Avec,C,b,X0,y0,Z0,objval,P] = control(B,solve), 
%%
%%  where B{k} = Bk,  k = 1:L. 
%%  
%%  For example,  B1 = [-0.8    1.2   -0.5;
%%                      -1.1   -1.0   -2.5;
%%                       2.0    0.2   -1.0];
%%
%%                B2 = [-1.5    0.5   -2.0; 
%%                       1.1   -2.0    0.2;
%%                      -1.4    1.1   -1.5];
%%
%% solve = 0 just to initialize
%%       = 1 if want to solve the problem
%%
%% SDPT3: version 3.0 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last modified: 2 Feb 01
%%******************************************************

  function [blk,Avec,C,b,X0,y0,Z0,objval,P] = control(B,solve); 

  if nargin < 2; solve = 0; end; 
  if ~iscell(B); error('B must be a cell array such that B{k} = Bk'); end; 

  L = size(B,2); 
  N = length(B{1});

  m = N*(N+1)/2 + 1; 
  n = N*(L+2); 

  blk{1,1} = 's'; blk{1,2} = N*ones(1,L+2); 

  C = sparse(n,n); 
  C(L*N+1:(L+1)*N,L*N+1:(L+1)*N) = speye(N);   
  b = [zeros(m-1,1); 1]; 

  A = cell(1,N*(N+1)/2+1);
  count = 1; 
  for k = 1:N  
      for j = k:N
          ak = []; aj = []; row1 = []; row2 = []; col1 = []; col2 = []; 
          for l = 1:L 
              G = B{l}; 
              ak = [ak; G(k,:)'];
              aj = [aj; G(j,:)'];
              col1 = [col1; ((l-1)*N+j)*ones(N,1)];
              col2 = [col2; ((l-1)*N+k)*ones(N,1)];
              row1 = [row1; [(l-1)*N+1:l*N]']; 
              row2 = [row2; [(l-1)*N+1:l*N]']; 
          end;    
          ek = [zeros(k-1,1); 0.5 ;zeros(N-k,1)];        
          ej = [zeros(j-1,1); 0.5 ;zeros(N-j,1)];        
          ak = [ak; ek; -ek];            
          aj = [aj; ej; -ej]; 
          col1 = [col1; (L*N+j)*ones(N,1); (L*N+N+j)*ones(N,1)];
          col2 = [col2; (L*N+k)*ones(N,1); (L*N+N+k)*ones(N,1)];
          row1 = [row1; [L*N+1:(L+1)*N]'; [L*N+N+1:(L+2)*N]' ]; 
          row2 = [row2; [L*N+1:(L+1)*N]'; [L*N+N+1:(L+2)*N]' ]; 
          if j ~= k; 
             tmp = sparse([row1; row2],[col1; col2],[ak; aj],n,n); 
          elseif (j == k);                           
             tmp = sparse(row1,col1,ak,n,n); 
          end; 
          A{count} = tmp + tmp'; 
          count = count+1; 
      end;
  end;
  tmp = sparse(n,n); 
  tmp((L+1)*N+1:(L+2)*N,(L+1)*N+1:(L+2)*N) = speye(N); 
  A{N*(N+1)/2+1} = tmp;
%%
%% initial iterate
%% 

  Avec = svec(blk,A,ones(size(blk,1),1)); 
  [X0,y0,Z0] = infeaspt(blk,Avec,C,b); 
%%
  if (solve)
     [obj,X,y,Z] = sqlp(blk,Avec,C,b,[],X0,y0,Z0);
     objval = mean(obj); 
     idx1 = 1; 
     for k = 1:N
         idx2 = idx1 + N-k; 
         P(k:N,k) = y(idx1 : idx2);  
         idx1 = idx2+1; 
     end;
     P = P + tril(P,-1)';
  else
     objval = []; P = [];
  end
%%======================================================
