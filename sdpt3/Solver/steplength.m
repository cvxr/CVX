%%***************************************************************************
%% steplength: compute xstep such that  X + xstep*dX >= 0.
%%
%% [xstep] = steplength(blk,X,dX,Xchol,invXchol);
%%
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%***************************************************************************

    function [xstep,invXchol] = steplength(blk,X,dX,Xchol,invXchol);

%%
    for p = 1:size(blk,1)
       pblk = blk(p,:); 
       numblk = length(pblk{2}); 
       pblksize = sum(pblk{2});
       if (any(isnan(dX{p})) | any(isinf(dX{p}))); xstep = 0; break; end; 
       if strcmp(pblk{1},'s') 
          if (max(pblk{2}) >= 200) 
             use_lanczos = 1; 
          else
             use_lanczos = 0; 
          end 
          if (use_lanczos)
             tol = 1e-3;     
             maxit = max(min(pblksize,30),round(sqrt(pblksize))); 
             [lam,delta,res] = lanczosfun(Xchol{p},-dX{p},maxit,tol);
             %%
             %% Note: lam <= actual largest eigenvalue <= lam + delta.
             %%   
	     d = lam+delta; 
          else
             if isempty(invXchol{p}); 
                invXchol{p} = inv(Xchol{p}); 
             end
	     tmp = Prod2(pblk,dX{p},invXchol{p},0); 
             M = Prod2(pblk,invXchol{p}',tmp,1); 
             if (exist('mexblkeig')==3)
                d = mexblkeig(pblk,-M); 
             else
                d = blkeig(pblk,-M); 
             end
          end
          tmp = max(d) + 1e-15*max(abs(d)); 
          if (tmp > 0);  
             xstep(p) = 1/max(tmp);
          else          
             xstep(p) = 1e12; 
          end
       elseif strcmp(pblk{1},'q')
          aa = qops(pblk,dX{p},dX{p},2); 
          bb = qops(pblk,dX{p},X{p},2); 
          cc = qops(pblk,X{p},X{p},2);
          dd = bb.*bb - aa.*cc; 
          tmp = min(aa,bb); 
          idx = find(dd > 0 & tmp < 0); 
          steptmp = 1e12*ones(numblk,1); 
          if ~isempty(idx)
             steptmp(idx) = -(bb(idx)+sqrt(dd(idx)))./aa(idx);       
          end
          idx = find(abs(aa) < eps & bb < 0); 
          if ~isempty(idx)
             steptmp(idx) = -cc(idx)./(2*bb(idx)); 
          end
          %%
          %% also need first component to be non-negative
          %%
          ss = 1 + [0, cumsum(pblk{2})];
          ss = ss(1:length(pblk{2})); 
          dX0 = dX{p}(ss); 
          X0 = X{p}(ss); 
          idx = find(dX0 < 0 & X0 > 0); 
          if ~isempty(idx)
             steptmp(idx) = min(steptmp(idx),-X0(idx)./dX0(idx)); 
          end
          xstep(p) = min(steptmp); 
       elseif strcmp(pblk{1},'l')
          idx = find(dX{p} < 0); 
          if ~isempty(idx)
             xstep(p) = min(-X{p}(idx)./dX{p}(idx));  
          else 
             xstep(p) = 1e12;
          end
       elseif strcmp(pblk{1},'u')
          xstep(p) = 1e12; 
       end
    end
    xstep = min(xstep); 
%%***************************************************************************
%%***************************************************************************
%% lanczos: find the largest eigenvalue of 
%%          invXchol'*dX*invXchol via the lanczos iteration.
%%
%% [lam,delta] = lanczosfun(Xchol,dX,maxit,tol,v)
%%
%% lam:  an estimate of the largest eigenvalue.
%% lam2: an estimate of the second largest eigenvalue.
%% res:  residual norm of the largest eigen-pair.
%% res2: residual norm of the second largest eigen-pair.
%%***************************************************************************

   function [lam,delta,res] = lanczosfun(Xchol,dX,maxit,tol,v) 
                                            
   if (norm(dX,'fro') < 1e-13) 
      lam = 0; delta = 0; res = 0; 
      return;
   end
   n = length(dX); 
   if (nargin < 5); 
      state = randn('state'); 
      randn('state',0); 
      v = randn(n,1); 
      randn('state',state); 
   end
   if (nargin < 4); tol = 1e-3; end
   if (nargin < 3); maxit = 30; end
   V = zeros(n,maxit+1); H = zeros(maxit+1,maxit); 
   v = v/norm(v); 
   V(:,1) = v; 
   if issparse(Xchol); Xcholtransp = Xchol'; end
%%
%% lanczos iteration. 
%%
   for k = 1:maxit
      if issparse(Xchol)
         w = dX*mextriangsp(Xcholtransp,v,1); 
         w = mextriangsp(Xchol,w,2); 
      else         
         w = dX*mextriang(Xchol,v,1); 
         w = mextriang(Xchol,w,2); 
      end      
      wold = w;
      if (k > 1); 
         w = w - H(k,k-1)*V(:,k-1); 
      end;
      alp = w'*V(:,k); 
      w   = w - alp*V(:,k); 
      H(k,k) = alp; 
      %%
      %% one step of iterative refinement if necessary. 
      %%
      if (norm(w) <= 0.8*norm(wold));
         s = (w'*V(:,1:k))'; 
         w = w - V(:,1:k)*s;
         H(1:k,k) = H(1:k,k) + s;
      end; 
      nrm = norm(w); 
      v = w/nrm; 
      V(:,k+1) = v; 
      H(k+1,k) = nrm;  H(k,k+1) = nrm; 
      %%
      %% compute ritz pairs and test for convergence
      %%
      if (rem(k,5) == 0) | (k == maxit); 
         Hk = H(1:k,1:k); Hk = 0.5*(Hk+Hk'); 
         [Y,D] = eig(Hk); 
         eigH  = real(diag(D)); 
         [dummy,idx] = sort(eigH);
         res_est = abs(H(k+1,k)*Y(k,idx(k)));
         if (res_est <= 0.1*tol) | (k == maxit);
            lam = eigH(idx(k));  
            lam2 = eigH(idx(k-1)); 
            z  = V(:,1:k)*Y(:,idx(k));
            z2 = V(:,1:k)*Y(:,idx(k-1));
            if issparse(Xchol) 
               tmp = dX*mextriangsp(Xcholtransp,z,1); 
               res = norm(mextriangsp(Xchol,tmp,2) -lam*z); 
               tmp = dX*mextriangsp(Xcholtransp,z2,1); 
               res2 = norm(mextriangsp(Xchol,tmp,2) -lam*z2);   
            else
               tmp = dX*mextriang(Xchol,z,1); 
               res = norm(mextriang(Xchol,tmp,2) -lam*z); 
               tmp = dX*mextriang(Xchol,z2,1); 
               res2 = norm(mextriang(Xchol,tmp,2) -lam*z2);   
            end
            tmp = lam-lam2 -res2; 
            if (tmp > 0); beta = tmp; else; beta = eps; end;  
            delta = min(res,res^2/beta); 
            if (delta <= tol); break; end;
         end 
      end 
   end
%%***************************************************************************
