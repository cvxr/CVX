%%*******************************************************************
%% randsdp.m : creates random feasible SDP problems with various block
%%              diagonal structures.
%%
%% [blk,Avec,C,b,X0,y0,Z0] = randsdp(dense_blk,sparse_blk,diag_blk,m,feas,solve);
%%
%% E.g.
%%      randsdp([32 20],[10 5],100,10,feas,solve);
%%
%%  dense_blk : for generating dense blocks, where
%%              dense_blk(i) is the dimension of the ith block
%%  sparse_blk: for generating a sparse block of small subblocks, where 
%%              sparse_blk(i) is the size of the ith subblock.
%%  diag_blk :  for generating a column vector of length specified by 
%%              diag_blk (this corresponds to a diagonal block). 
%%
%%  feas = 2 if want centered feasible starting point
%%       = 1 if want feasible starting point
%%       = 0 if otherwise.  (default)
%%
%%  solve = 0 just to initialize  (default)
%%        = 1 if want to solve the problem. 
%%
%% SDPT3: version 3.0 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last modified: 2 Feb 01
%%******************************************************************

  function  [blk,Avec,C,b,X0,y0,Z0] = ...
                randsdp(dense_blk,sparse_blk,diag_blk,m,feas,solve);

   if nargin < 6; solve = 0; end;
   if nargin < 5; feas  = 0; end;
   
   blk = [];
   if ~isempty(dense_blk);  
      for k = 1:length(dense_blk); 
          blk{k,1} = 's';  blk{k,2} = dense_blk(k);
      end; 
   end;
   if ~isempty(sparse_blk); 
      if size(sparse_blk,1) > size(sparse_blk,2); sparse_blk = sparse_blk'; end; 
      tmp = size(blk,1); 
      blk{tmp+1,1} = 's';  
      blk{tmp+1,2} = sparse_blk; 
   end;
   if ~isempty(diag_blk); 
      tmp = size(blk,1); 
      blk{tmp+1,1} = 'l';        
      blk{tmp+1,2} = diag_blk;  
   end; 
   N = size(blk,1); 
%%
%% generate X0, Z0 
%%
   tmp_sp = sparse(sum(sparse_blk));  
   for L = 1:2
       T = [];
       if ~isempty(dense_blk); 
          for k = 1:length(dense_blk);
              n = dense_blk(k); 
              tmp = randn(n); 
              tmp = tmp*tmp';
              tmp = 0.5*(tmp + tmp');
              %mineig = min(real(eig(tmp)));
              %if (mineig < 0); tmp = tmp - 1.1*mineig*eye(n); end;
              if (feas == 2);  tmp = eye(n); end; 
              T{k,1} = tmp; 
          end;
       end;
       if ~isempty(sparse_blk); 
          tmp_sp = sparse(sum(sparse_blk),sum(sparse_blk)); 
          for k = 1:length(sparse_blk);
              n = sparse_blk(k);
              pos = [sum(sparse_blk(1:k-1))+1 : sum(sparse_blk(1:k))]; 
              tmp = randn(n);  
              tmp = tmp*tmp';
              tmp = 0.5*(tmp + tmp');
              %mineig = min(real(eig(tmp)));
              %if (mineig < 0); tmp = tmp - 1.1*mineig*eye(n); end;
              if (feas == 2);  tmp = eye(n); end; 
              tmp_sp(pos,pos) = sparse(tmp); 
          end; 
          T{size(T,1)+1,1} = tmp_sp; 
       end;
       if ~isempty(diag_blk);           
          if (feas == 2);  
             T{size(T,1)+1,1} = ones(diag_blk,1); 
          else; 
             tmp = randn(diag_blk,1);
             tmp = tmp.*tmp;
             %mineig = min(tmp);
             %if (mineig < 0); tmp = tmp - 1.1*mineig*ones(diag_blk,1); end;
             T{size(T,1)+1,1} = tmp; 
          end;
       end; 
       if (L == 1);  X0 = T; else; Z0 = T; end;
   end;
%%
%% set up the matrices Ak and b
%%
   b  = zeros(m,1);
   y0 = randn(m,1); 
   A  = cell(N,m);  
   Ak_sp = sparse(sum(sparse_blk),sum(sparse_blk)); 
   for k = 1:m;
       Ak = []; 
       if ~isempty(dense_blk); 
          for j = 1:length(dense_blk); 
              tmp = randn(dense_blk(j)); tmp = 0.5*(tmp+tmp');
              Ak{j,1} = tmp;
          end;
       end;
       if ~isempty(sparse_blk);   
          for j = 1:length(sparse_blk);
              n = sparse_blk(j);
              tmp = randn(n); tmp = 0.5*(tmp+tmp');
              pos = [sum(sparse_blk(1:j-1))+1 : sum(sparse_blk(1:j))]; 
              Ak_sp(pos,pos) = tmp; 
          end;
          Ak{size(Ak,1)+1,1} = Ak_sp; 
       end;
       if ~isempty(diag_blk);  
          Ak{size(Ak,1)+1,1} = randn(diag_blk,1); 
       end;
       A(:,k) = Ak;
       b(k) = blktrace(blk,Ak,X0); 
   end;
   C = ops(Z0,'+',Asum(blk,A,y0)); 
%%
%% infeasible initial iterate
%%
   Avec = svec(blk,A,ones(size(blk,1),1)); 
   if (feas == 0); 
      [X0,y0,Z0] = infeaspt(blk,Avec,C,b);    
   end;
%%
   if (solve); 
      [obj,X,y,Z] = sqlp(blk,Avec,C,b,[],X0,y0,Z0);
   end; 
%%=================================================






