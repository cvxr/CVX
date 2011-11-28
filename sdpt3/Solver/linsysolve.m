%%***************************************************************
%% linsysolve: solve linear system to get dy, and direction
%%             corresponding to unrestricted variables. 
%%
%% [xx,coeff,L,resnrm] = linsysolve(schur,UU,Afree,EE,rhs); 
%%
%% child functions: symqmr.m, mybicgstable.m, linsysolvefun.m
%%
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%***************************************************************
   
   function [xx,coeff,L,resnrm] = linsysolve(par,schur,UU,Afree,EE,rhs); 
   
    global solve_ok exist_analytic_term
    global nnzmat nnzmatold matfct_options matfct_options_old use_LU 
    global msg diagR diagRold numpertdiagschur

    spdensity  = par.spdensity;
    printlevel = par.printlevel; 
    iter       = par.iter; 
    if isfield(par,'relgap') & isfield(par,'pinfeas') & isfield(par,'dinfeas')
       err = max([par.relgap,par.pinfeas,par.dinfeas]);  
    else
       err = inf;
    end
%%
    m = length(schur); 
    if (iter==1); 
       use_LU = 0; matfct_options_old = ''; 
       diagR = ones(m,1); 
       numpertdiagschur = 0; 
    end
    if isempty(nnzmatold); nnzmatold = 0; end
    diagRold = diagR; 
%%
%% schur = schur + rho*diagschur + lam*AAt
%%
    diagschur = abs(full(diag(schur))); 
    if (par.ublksize) 
       minrho(1) = 1e-15; 
    else
       minrho(1) = 1e-17;
    end
    minrho(1) = max(minrho(1), 1e-6/3.0^iter); %% old: 1e-6/3.0^iter
    minrho(2) = max(1e-04, 0.7^iter);
    minlam = max(1e-10, 1e-4/2.0^iter); 
    rho = min(minrho(1), minrho(2)*(1+norm(rhs))/(1+norm(diagschur.*par.y)));
    lam = min(minlam, 0.1*rho*norm(diagschur)/par.normAAt);
    if (exist_analytic_term); rho = 0; end; %% important
    ratio = max(diagR)/min(diagR); 
    if (par.depconstr) | (ratio > 1e10) | (iter < 5)
       %% important: do not perturb beyond certain threshold 
       %% since it will adversely affect prim_infeas of fp43
       %%
       pertdiagschur = min(rho*diagschur,1e-4./max(1,abs(par.dy)));
       mexschurfun(schur,full(pertdiagschur));
       %%if (printlevel>2); fprintf(' %2.1e',rho); end
    end
    if (par.depconstr) | (par.ZpATynorm > 1e10) | (par.ublksize) | (iter < 10)
       %% Note: do not add this perturbation even if ratio is large.
       %% It adversely affects hinf15.
       %%       
       lam = min(lam,1e-4/max(1,norm(par.AAt*par.dy)));
       if (exist_analytic_term); lam = 0; end
       mexschurfun(schur,lam*par.AAt); 
       %%if (printlevel>2); fprintf('*'); end
    end
    if (max(diagschur)/min(diagschur) > 1e14) & (par.blkdim(2) == 0) ...
       & (iter > 10)
       tol = 1e-8; 
       idx = find(diagschur < tol); len = length(idx);
       pertdiagschur = zeros(m,1);  
       if (len > 0 & len < 5) & (norm(rhs(idx)) < tol) 
          pertdiagschur(idx) = 1*ones(length(idx),1); 
          mexschurfun(schur,pertdiagschur); 
          numpertdiagschur = numpertdiagschur + 1; 
          if (printlevel>2); fprintf('#'); end   
       end
    end
%% 
%% assemble coefficient matrix
%% 
    len = size(Afree,2);
    if ~isempty(EE)
       EE(:,[1 2]) = len + EE(:,[1 2]); %% adjust for ublk
    end
    EE = [(1:len)' (1:len)' zeros(len,1); EE]; 
    if isempty(EE)
       coeff.mat22 = []; 
    else
       coeff.mat22 = spconvert(EE);
    end
    if (size(Afree,2) | size(UU,2))
       coeff.mat12 = [Afree, UU]; 
    else
       coeff.mat12 = []; 
    end
    coeff.mat11 = schur; %% important to use perturbed schur matrix
    ncolU = size(coeff.mat12,2); 
%%
%% pad rhs with zero vector
%% decide which solution methods to use
%%
    rhs = [rhs; zeros(m+ncolU-length(rhs),1)]; 
    if (ncolU > 300); use_LU = 1; end
%%
%% Cholesky factorization
%%
    L = []; resnrm = norm(rhs); xx = inf*ones(m,1);
    if (~use_LU)
       solve_ok = 1;  solvesys = 1;    
       nnzmat = mexnnz(coeff.mat11);
       nnzmatdiff = (nnzmat ~= nnzmatold);     
       if (nnzmat > spdensity*m^2) | (m < 500) 
          matfct_options = 'chol';
       else
          matfct_options = 'spchol'; 
       end
       if (printlevel>2); fprintf(' %s ',matfct_options); end 
       L.matdim = length(schur); 
       if strcmp(matfct_options,'chol')
          if issparse(schur); schur = full(schur); end;
          if (iter<=5); %%--- to fix strange anonmaly in Matlab
             mexschurfun(schur,1e-20,2); 
          end 
          L.matfct_options = 'chol';
          [L.R,indef] = chol(schur);
          L.perm = [1:m]; 
          diagR  = diag(L.R).^2; 
       elseif strcmp(matfct_options,'spchol')
          if ~issparse(schur); schur = sparse(schur); end;
          L.matfct_options = 'spchol'; 
          [L.R,indef,L.perm] = chol(schur,'vector'); 
          L.Rt  = L.R';
          diagR = full(diag(L.R)).^2; 
       end    
       if (indef)
          diagR = diagRold; 
          solve_ok = -2; solvesys = 0;
          msg = 'linsysolve: Schur complement matrix not positive definite'; 
          if (printlevel); fprintf('\n  %s',msg); end
       end
       if (solvesys)
          if (ncolU)
             tmp = coeff.mat12'*linsysolvefun(L,coeff.mat12)-coeff.mat22; 
	     if issparse(tmp); tmp = full(tmp); end
             tmp = 0.5*(tmp + tmp'); 
             [L.Ml,L.Mu,L.Mp] = lu(tmp);
             tol = 1e-16; 
             condest = max(abs(diag(L.Mu)))/min(abs(diag(L.Mu))); 
             idx = find(abs(diag(L.Mu)) < tol);
             if ~isempty(idx) | (condest > 1e30);  %%old: 1e18 
                solvesys = 0; solve_ok = -4;  
                use_LU = 1; 
                msg = 'SMW too ill-conditioned, switch to LU factor'; 
                if (printlevel); fprintf('\n  %s, %2.1e.',msg,condest); end
             end         
          end
          [xx,resnrm,solve_ok] = symqmr(coeff,rhs,L,[],[],printlevel);
          if (solve_ok <= 0.3) & (printlevel)
             fprintf('\n  warning: symqmr failed: %3.1f ',solve_ok); 
          end
       end
       if (solve_ok <= 0.3) 
          tol = 1e-10; 
          if (m < 1e4 & strcmp(matfct_options,'chol') & (err > tol)) ...
             | (m < 2e5 & strcmp(matfct_options,'spchol') & (err > tol)) 
             use_LU = 1;
             if (printlevel); fprintf('\n  switch to LU factor.'); end
          end
       end
    end
%%
%% LU factorization
%%
    if (use_LU)
       nnzmat = mexnnz(coeff.mat11)+mexnnz(coeff.mat12); 
       nnzmatdiff = (nnzmat ~= nnzmatold);  
       solve_ok = 1; solvesys = 1; 
       if ~isempty(coeff.mat22)
          raugmat = [coeff.mat11, coeff.mat12; coeff.mat12', coeff.mat22]; 
       else
          raugmat = coeff.mat11; 
       end
       if (nnzmat > spdensity*m^2) | (m+ncolU < 500) 
          matfct_options = 'lu';  %% lu is better than ldl
       else
          matfct_options = 'splu'; %% faster than spldl
       end
       if (printlevel>2); fprintf(' %s ',matfct_options); end 
       L.matdim = length(raugmat); 
       if strcmp(matfct_options,'lu') 
          if issparse(raugmat); raugmat = full(raugmat); end
          L.matfct_options = 'lu';         
          [L.L,L.U,L.p] = lu(raugmat,'vector'); 
       elseif strcmp(matfct_options,'splu')
          if ~issparse(raugmat); raugmat = sparse(raugmat); end  
          L.matfct_options = 'splu';
          [L.L,L.U,L.p,L.q,L.s] = lu(raugmat,'vector'); 
          L.s = full(diag(L.s)); 
       elseif strcmp(matfct_options,'ldl') 
          if issparse(raugmat); raugmat = full(raugmat); end
          L.matfct_options = 'ldl'; 
          [L.L,L.D,L.p] = ldl(raugmat,'vector');   
          L.D = sparse(L.D); 
       elseif strcmp(matfct_options,'spldl') 
          if ~issparse(raugmat); raugmat = sparse(raugmat); end  
          L.matfct_options = 'spldl'; 
          [L.L,L.D,L.p,L.s] = ldl(raugmat,'vector');
          L.s  = full(diag(L.s));           
          L.Lt = L.L'; 
       end
       if (solvesys)
          [xx,resnrm,solve_ok] = symqmr(coeff,rhs,L,[],[],printlevel);
          %%[xx,resnrm,solve_ok] = mybicgstab(coeff,rhs,L,[],[],printlevel);
          if (solve_ok<=0) & (printlevel)
             fprintf('\n  warning: bicgstab fails: %3.1f,',solve_ok); 
          end
       end
    end
    if (printlevel>2); fprintf('%2.0d ',length(resnrm)-1); end
%%
    nnzmatold = nnzmat; matfct_options_old = matfct_options; 
%%***************************************************************
