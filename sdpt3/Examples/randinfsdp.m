%%*******************************************************************
%% randinfsdp.m : creates random infeasible SDP problems with various block
%%                diagonal structures. 
%%
%% [blk,Avec,C,b,X0,y0,Z0] = 
%%      randinfsdp(dense_blk,sparse_blk,diag_blk,m,infeas,solve);
%%
%% E.g.
%%      randinfsdp([32 20],[10 5],100,10,infeas,solve);
%%
%%  dense_blk : for generating dense blocks, where
%%              dense_blk(i) is the dimension of the ith block
%%  sparse_blk: for generating a sparse block of small subblocks, where 
%%              sparse_blk(i) is the size of the ith subblock.
%%  diag_blk: for generating a column vector of length specified by 
%%            diag_blk (this corresponds to a diagonal block). 
%%
%%  infeas = 1 if want primal infeasible pair of problems
%%         = 2 if want dual infeasible pair of problems
%%
%%  solve = 0 just to initialize
%%        = 1 if want to solve the problem. 
%%
%% SDPT3: version 3.0 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last modified: 2 Feb 01
%%******************************************************************

   function  [blk,Avec,C,b,X0,y0,Z0] = ...
                randinfsdp(dense_blk,sparse_blk,diag_blk,m,infeas,solve);

   if nargin < 6; solve = 0; end;
   if nargin < 5; error(' insufficient number of inputs '); end;   
   if all(infeas-[1 2]); 
      error(' infeas must be 1 or 2 ');
   end

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

   if (infeas == 1),
%%
%% primal infeasible
%%
%% generate infeasibility certificate Zi and Z0
%%
   tmp_sp = sparse(sum(sparse_blk));  
   for L = 1:2,
       T = [];
       if ~isempty(dense_blk); 
          for k = 1:length(dense_blk);
              n = dense_blk(k); 
              tmp = randn(n); 
              tmp = tmp*tmp'; 
              tmp = 0.5*(tmp + tmp');
              %mineig = min(real(eig(tmp)));
              %if (mineig < 0); tmp = tmp - 1.1*mineig*eye(n); end;
              T{k,1} = tmp; 
          end;
       end;
       if ~isempty(sparse_blk); 
          for k = 1:length(sparse_blk);
              n = sparse_blk(k);
              pos = [sum(sparse_blk(1:k-1))+1 : sum(sparse_blk(1:k))]; 
              tmp = randn(n);  
              tmp = tmp*tmp'; 
              tmp = 0.5*(tmp + tmp');
              %mineig = min(real(eig(tmp)));
              %if (mineig < 0); tmp = tmp - 1.1*mineig*eye(n); end;
              tmp_sp(pos,pos) = sparse(tmp); 
          end; 
          T{size(T,1)+1,1} = tmp_sp; 
       end;
       if ~isempty(diag_blk);           
          tmp = randn(diag_blk,1); 
          tmp = tmp.*tmp;
          %mineig = min(tmp);
          %if (mineig < 0); tmp = tmp - 1.1*mineig*ones(diag_blk,1); end;
          T{size(T,1)+1,1} = tmp; 
       end; 
       if (L == 1);  Zi = T; else; Z0 = T; end;
   end;
%%
%% set up the matrices Ak and b
%%
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
   end;

   y = ones(m,1);
   SAZm = ops(ops(Zi,'+',Asum(blk,A,y)),'/',m);
   SAZm = ops(ops(SAZm,'+',ops(SAZm,'transpose')),'*',0.5);
   yi  = randn(m,1);
   for k = 1:m,
      Ak = A(:,k);
      Ak = ops(ops(Ak,'-',SAZm),'/',yi(k));
      Ak = ops(ops(Ak,'+',ops(Ak,'transpose')),'*',0.5);
      A(:,k) = Ak;
   end;
   b = randn(m,1);
   if (b'*yi < 0), b = -b; end;
   y0  = randn(m,1);
   C = ops(Z0,'+',Asum(blk,A,y0));
   C = ops(ops(C,'+',ops(C,'transpose')),'*',0.5);
%%
%% (yi,Zi) is a primal infeasibility certificate
%%
   elseif (infeas == 2),
%%
%% dual infeasible
%%
%% generate infeasibility certificate Xi, positive definite X0, and C
%%
   tmp_sp_Xi  = sparse(sum(sparse_blk));
   tmp_sp_X0 = sparse(sum(sparse_blk));
   tmp_sp_C   = sparse(sum(sparse_blk));
   Xi = [];
   X0 = [];
   C  = [];
       if ~isempty(dense_blk);
          for k = 1:length(dense_blk);
              n = dense_blk(k);
              tmp = randn(n); 
              tmp = tmp*tmp'; 
              tmp = 0.5*(tmp + tmp');
              %mineig = min(real(eig(tmp)));
              %if (mineig < 0); tmp = tmp - 1.1*mineig*eye(n); end;
              Xi{k,1} = tmp;
              tmp = randn(n); 
              tmp = tmp*tmp'; 
              tmp = 0.5*(tmp + tmp');
              %mineig = min(real(eig(tmp)));
              %if (mineig < 0); tmp = tmp - 1.1*mineig*eye(n); end;
              X0{k,1} = tmp;
              tmp = randn(n); tmp = 0.5*(tmp + tmp');
              C{k,1} = tmp;
          end;
       end;
       if ~isempty(sparse_blk);
          for k = 1:length(sparse_blk);
              n = sparse_blk(k);
              pos = [sum(sparse_blk(1:k-1))+1 : sum(sparse_blk(1:k))];
              tmp = randn(n);  
              tmp = tmp*tmp'; 
              tmp = 0.5*(tmp + tmp');
              %mineig = min(real(eig(tmp)));
              %if (mineig < 0); tmp = tmp - 1.1*mineig*eye(n); end;
              tmp_sp_Xi(pos,pos) = sparse(tmp);
              tmp = randn(n);  
              tmp = tmp*tmp'; 
              tmp = 0.5*(tmp + tmp');
              %mineig = min(real(eig(tmp)));
              %if (mineig < 0); tmp = tmp - 1.1*mineig*eye(n); end;
              tmp_sp_X0(pos,pos) = sparse(tmp);
              tmp = randn(n);  tmp = 0.5*(tmp + tmp');
              tmp_sp_C(pos,pos) = sparse(tmp);
          end;
          Xi{size(Xi,1)+1,1} = tmp_sp_Xi;
          X0{size(X0,1)+1,1} = tmp_sp_X0;
          C{size(C,1)+1,1}   = tmp_sp_C;
       end;
       if ~isempty(diag_blk);
          n = diag_blk;
          tmp = randn(n,1); 
          tmp = tmp.*tmp;
          %mineig = min(tmp);
          %if (mineig < 0); tmp = tmp - 1.1*mineig*ones(n,1); end;
          Xi{size(X0,1)+1,1} = tmp;
          tmp = randn(n,1); 
          tmp = tmp.*tmp;
          %mineig = min(tmp);
          %if (mineig < 0); tmp = tmp - 1.1*mineig*ones(n,1); end;
          X0{size(X0,1)+1,1} = tmp;
          tmp = randn(n,1);
          C{size(C,1)+1,1} = tmp;
       end;
%%
%% set up the matrices Ak and b
%%
   trXX = blktrace(blk,Xi,Xi);
   A  = cell(N,m);
   Ak_sp = sparse(sum(sparse_blk),sum(sparse_blk));
   b = ones(m,1);
   for k = 1:m;
       Ak = [];
       if ~isempty(dense_blk);
          for j = 1:length(dense_blk);
              n = dense_blk(j);
              Aj  = randn(n);
              Ak{j,1} = (Aj + Aj')/2;
          end;
       end;
       if ~isempty(sparse_blk);
          for j = 1:length(sparse_blk);
              pos = [sum(sparse_blk(1:j-1))+1 : sum(sparse_blk(1:j))];
              n = sparse_blk(j);
              Aj  = randn(n);
              Ak_sp(pos,pos) = (Aj + Aj')/2;
          end;
          Ak{size(Ak,1)+1,1} = Ak_sp;
       end;
       if ~isempty(diag_blk);
          n = diag_blk;
          Ak{size(Ak,1)+1,1} = randn(n,1);
       end;
       trAX = blktrace(blk,Ak,Xi);
       Ak = ops(Ak,'-',ops(Xi,'*',(trAX/trXX)));
       Ak = ops(ops(Ak,'+',ops(Ak,'transpose')),'*',0.5);
       A(:,k) = Ak;
       b(k) = blktrace(blk,Ak,X0);
   end;
   if (blktrace(blk,C,Xi) > 0), C = ops(C,'*',-1); end;
%%
%% Xi is a dual infeasibility certificate
%%
   end;

%%
%% infeasible initial iterate
%%
   Avec = svec(blk,A,ones(size(blk,1),1));    
   [X0,y0,Z0] = infeaspt(blk,Avec,C,b);    
%%
   if solve; 
      [obj,X,y,Z] = sqlp(blk,Avec,C,b,[],X0,y0,Z0);
   end; 

%%=================================================
