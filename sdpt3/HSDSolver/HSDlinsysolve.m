%%***************************************************************
%% linsysolve: solve linear system to get dy, and direction
%%             corresponding to unrestricted variables. 
%%
%% [xx,coeff,L,resnrm] = linsysolve(schur,UU,EE,Bmat,rhs); 
%%
%% child functions: mybicgstable.m
%%
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%***************************************************************
   
   function [xx,coeff,L,resnrm] = HSDlinsysolve(par,schur,UU,EE,Bmat,rhs); 
   
    global solve_ok msg
    global nnzmat nnzmatold matfct_options matfct_options_old use_LU
    global Lsymb 

    spdensity  = par.spdensity;
    printlevel = par.printlevel;
    iter       = par.iter; 

    m = length(schur); 
    if (iter==1); use_LU = 0; matfct_options_old = ''; end
    if isempty(nnzmatold); nnzmatold = 0; end
%%
%% diagonal perturbation
%% old: pertdiag = 1e-15*max(1,diagschur);
%% 
    diagschur = abs(full(diag(schur)));
    const = 1e-2/max(1,norm(par.dy2)); 
    alpha = max(1e-14,min(1e-10,const*norm(par.rp))/(1+norm(diagschur.*par.dy2)));
    pertdiag = alpha*max(1e-8,diagschur); %% Note: alpha is close to 1e-15.
    mexschurfun(schur,pertdiag); 
    %%if (printlevel); fprintf(' %3.1e ',alpha); end
    if (par.depconstr) | (min(diagschur) < min([1e-20*max(diagschur), 1e-4])) 
       lambda = 0.1*min(1e-14,const*norm(par.rp)/(1+norm(par.diagAAt.*par.dy2)));
       mexschurfun(schur,lambda*par.diagAAt);  
       if (printlevel); fprintf('#'); end
    end
    if (max(diagschur)/min(diagschur) > 1e14) & (par.blkdim(2) == 0) ...
       & (iter > 10)
       tol = 1e-6; 
       idx = find(diagschur < tol); len = length(idx);
       pertdiagschur = zeros(m,1);  
       if (len > 0 & len < 5) & (norm(rhs(idx)) < tol) 
          pertdiagschur(idx) = 1*ones(length(idx),1); 
          mexschurfun(schur,pertdiagschur); 
          if (printlevel); fprintf('$'); end
       end
    end
%%
%%
%%
    UU = [UU, Bmat]; 
    if ~isempty(EE)
       len = max(max(EE(:,1)),max(EE(:,2))); 
    else
       len = 0;
    end
    tmp = [len+1,len+3,-1; len+2,len+4,1; len+3,len+1,1; len+4,len+2,-1; 
           len+2,len+2,par.addschur];  %% this is the -inverse
    EE = [EE; tmp]; 
    ncolU = size(UU,2);
%%
%% assemble coefficient matrix
%% 
    if isempty(EE)
       coeff.mat22 = []; 
    else
       coeff.mat22 = spconvert(EE);
    end
    coeff.mat12 = UU; 
    coeff.mat11 = schur; %% important to use perturbed schur matrix
%%
%% pad rhs with zero vector
%% decide which solution methods to use
%%
    rhs = [rhs; zeros(m+ncolU-length(rhs),1)]; 
    if (ncolU > 300); use_LU = 1; end
%%
%% Cholesky factorization
%%
    L = []; resnrm = []; xx = inf*ones(m,1);
    if (~use_LU)
       nnzmat = mexnnz(coeff.mat11);
       nnzmatdiff = (nnzmat ~= nnzmatold);     
       solve_ok = 1;  solvesys = 1;    
       if (nnzmat > spdensity*m^2) | (m < 500) 
          matfct_options = 'chol';
       else
          matfct_options = 'spchol'; 
       end
       if (printlevel > 2); fprintf(' %s',matfct_options); end 
       if strcmp(matfct_options,'chol')
          if issparse(schur); schur = full(schur); end;
          if (iter<=5); %%--- to fix strange anonmaly in Matlab
             mexschurfun(schur,1e-20,2); 
          end 
          L.matfct_options = 'chol';    
          [L.R,indef] = chol(schur); 
          L.perm = [1:m];
      elseif strcmp(matfct_options,'spchol')
          if ~issparse(schur); schur = sparse(schur); end;
          L.matfct_options = 'spchol'; 
          [L.R,indef,L.perm] = chol(schur,'vector'); 
          L.Rt = L.R'; 
       end    
       if (indef)
          solve_ok = -2; solvesys = 0; 
          msg = 'HSDlinsysolve: Schur complement matrix not positive definite'; 
          if (printlevel); fprintf('\n  %s',msg); end
       end
       if (solvesys)
          if (ncolU)
             tmp = coeff.mat12'*linsysolvefun(L,coeff.mat12)-coeff.mat22;
	     if issparse(tmp); tmp = full(tmp); end
             [L.Ml,L.Mu,L.Mp] = lu(tmp);
             tol = 1e-16;
             condest = max(abs(diag(L.Mu)))/min(abs(diag(L.Mu))); 
             idx = find(abs(diag(L.Mu)) < tol);
             if ~isempty(idx) | (condest > 1e30*sqrt(norm(par.diagAAt))); %% default=1e25 
                solvesys = 0; solve_ok = -4; 
                use_LU = 1; 
                msg = 'SMW too ill-conditioned, switch to LDL factor'; 
                if (printlevel); fprintf('\n  %s, %2.1e.',msg,condest); end
             end         
          end
          if (solvesys)
             [xx,resnrm,solve_ok] = HSDbicgstab(coeff,rhs,L,[],[],printlevel);
             if (solve_ok<=0) & (printlevel)
                fprintf('\n  warning: HSDbicgstab fails: %3.1f.',solve_ok); 
             end
          end
       end
       if (solve_ok < 0) 
          if (m < 6000 & strcmp(matfct_options,'chol')) | ...
             (m < 1e5 & strcmp(matfct_options,'spchol')) 
             use_LU = 1;
             if (printlevel); fprintf('\n  switch to LDL factor'); end
          end
       end
    end
%%
%% symmetric indefinite LDL factorization
%%
    if (use_LU)
       nnzmat = mexnnz(coeff.mat11)+mexnnz(coeff.mat12); 
       nnzmatdiff = (nnzmat ~= nnzmatold);  
       solve_ok = 1; 
       if ~isempty(coeff.mat22)
          raugmat = [coeff.mat11, coeff.mat12; coeff.mat12', coeff.mat22]; 
       else
          raugmat = coeff.mat11; 
       end
       if (nnzmat > spdensity*m^2) | (m+ncolU < 500) 
          matfct_options = 'ldl'; 
       else
          matfct_options = 'spldl';
       end
       if (printlevel > 2); fprintf(' %s ',matfct_options); end 
       if strcmp(matfct_options,'ldl') 
          if issparse(raugmat); raugmat = full(raugmat); end
          L.matfct_options = 'ldl'; 
          [L.L,L.D] = ldl(raugmat); 
          L.perm = [1:length(raugmat)]';        
       elseif strcmp(matfct_options,'spldl') 
          if ~issparse(raugmat); raugmat = sparse(raugmat); end  
          L.matfct_options = 'spldl';  
          [L.L,L.D,L.perm] = ldl(raugmat,'vector');
          L.Lt = L.L';
       end
       [xx,resnrm,solve_ok] = HSDbicgstab(coeff,rhs,L,[],[],printlevel);
       if (solve_ok<=0) & (printlevel)
          fprintf('\n  warning: HSDbicgstab fails: %3.1f,',solve_ok); 
       end
    end
    if (printlevel>=3); fprintf('%2.0f ',length(resnrm)-1); end
%%
    nnzmatold = nnzmat; matfct_options_old = matfct_options; 
%%***************************************************************
