%%***************************************************************
%% degeneracy: determine if an SDP problem is non-degenerate. 
%%
%% [ddx,ddz,B1,B2,sig1,sig12] = degeneracy(blk,At,X,y,Z);
%%
%% Assume that strict complementarity holds:
%% for primal non-degeneracy, we need rank([B1 B2]) = m
%% for dual non-degeneracy,   we need B1 to have full column rank. 
%%
%%***************************************************************

    function [ddx,ddz,XB1,XB2,ZB1,sig1,sig12] = degeneracy(blk,At,X,y,Z); 

    if ~iscell(X); tmp = X; clear X; X{1} = tmp; end;
    if ~iscell(Z); tmp = Z; clear Z; Z{1} = tmp; end;

    m = length(y); XB1 = []; XB2 = []; ZB1 = []; 
    numblk = size(blk,1); 
%%
%%
%%
    for p = 1:size(blk,1)
       n = length(X{p}); 
       if strcmp(blk{p,1},'s')
          mu = sum(sum(X{p}.*Z{p}))/n; 
          [Qx,Dx] = eig(full(X{p})); [dx,idx] = sort(diag(Dx)); 
          Qx = Qx(:,idx(n:-1:1)); dx = dx(n:-1:1); dx = dx + max(dx)*eps;  
          Dx = diag(dx); 
          [Qz,Dz] = eig(full(Z{p})); [dz,idx] = sort(diag(Dz)); 
          Qz = Qz(:,idx); dz = dz + max(dz)*eps;  
          Dz = diag(dz);   
          sep_option = 1;
          if (sep_option==1)
             tolx(p) = mean(sqrt(dx.*dz));  
             tolz(p) = tolx(p); 
	  elseif (sep_option==2)
             ddtmp = dx./dz;  
             idxtmp = find(ddtmp<1e12); 
             len = max(3,min(idxtmp)); 
             ddratio = ddtmp(len:n-1)./ddtmp(len+1:n); 
             [dummy,idxmax] = max(ddratio); 
             idxmax2 = idxtmp(find(idxtmp==len+idxmax-1)); 
             tmptolx = mean(dx(idxmax2+[0:1])); 
             tmptolz = mean(dz(idxmax2+[0:1])); 
             tolx(p) = exp(mean(log([tmptolx tmptolz])));
             tolz(p) = tolx(p); 
          end
          idxr = find(dx > tolx(p));  rp = length(idxr); 
          idxs = find(dz > tolz(p));  sp = length(idxs); 
          r(p) = rp; s(p) = sp; 
          prim_rank(p)    = n*(n+1)/2 - (n-rp)*(n-rp+1)/2;
          dual_rank(p)   = (n-sp)*(n-sp+1)/2;
          strict_comp(p) = (r(p)+s(p) == n); 
          if (nargout > 2)    
             Q1 = Qx(:,idxr); Q2 = Qx(:,setdiff([1:n],idxr)); 
             B11 = zeros(m,rp*(rp+1)/2); 
             B22 = zeros(m,rp*(n-rp));                  
             for k = 1:m
		Ak = smat(blk(p,:),At{p}(:,k)); 
                tmp = Q1'*Ak*Q2;
                B22(k,:) = tmp(:)';
                tmp = Q1'*Ak*Q1;
                B11(k,:) = svec(blk(p,:),tmp)';
             end
	     XB1 = [XB1, B11]; 
	     XB2 = [XB2, sqrt(2)*B22];  
             Qz1 = Qz(:,setdiff([1:n],idxs));
             ZB11 = zeros(m,(n-sp)*(n-sp+1)/2);
             for k = 1:m
		Ak = smat(blk(p,:),At{p}(:,k)); 
                tmp = Qz1'*Ak*Qz1;
                ZB11(k,:) = svec(blk(p,:),tmp)';
             end 
	     ZB1 = [ZB1, ZB11];
          end                    
       elseif strcmp(blk{p,1},'q')
          error('qblk is not allowed at the moment.'); 
       elseif strcmp(blk{p,1},'l')
          mu = X{p}'*Z{p}/n; 
          dx = sort(X{p}); dx = dx(n:-1:1);
          dz = sort(Z{p});
          tolx(p) = mean(sqrt(dx.*dz));
          tolz(p) = tolx(p);   
          idxr = find(dx > tolx(p));  rp = length(idxr); 
          idxs = find(dz > tolz(p));  sp = length(idxs); 
          r(p) = rp; s(p) = sp; 
          prim_rank(p)    = rp; 
          dual_rank(p)    = n-sp; 
          strict_comp(p) = (r(p)+s(p) == n);          
          if (nargout > 2)
             idx = find(X{p} > tolx(p)); 
	     XB1 = [XB1, full(At{p}(idx,:))']; 
             zidx = find(Z{p} < tolz(p)); 
             ZB1 = [ZB1, full(At{p}(zidx,:))'];  
          end
       end
       ddx{p} = dx; ddz{p} = dz; 
       fprintf('\n    blkno = %1.0d, tol = %3.1e,%3.1e,  m = %2.0d',...
               p,tolx(p),tolz(p),m); 
       fprintf('\n    n= %2.0d,  r= %2.0d, s= %2.0d',n,r(p),s(p)); 
       fprintf('\n    n2-(n-r)2 = %2.0d',prim_rank(p));
       fprintf('\n    (n-s)2    = %2.0d',dual_rank(p));
       fprintf('\n    complemen = %2.1e  %2.1e\n',max(dx+dz),min(dx+dz));
       subplot(121) 
       color = (1-p/numblk)*[0 0 1] + (p/numblk)*[1 0 0];
       semilogy(dx,'+','color',color); hold on; 
       semilogy(dz,'o','color',color); 
       semilogy([1 n],tolx(p)*[1 1],'color',color);
       semilogy([1 n],tolz(p)*[1 1],'--','color',color);
       title('eig(X) and eig(Z)'); 
       subplot(122)
       semilogy(dx+dz,'*','color',color); hold on; 
       title('dx+dz')
    end
%%
%%
%%
    subplot(121); hold off; subplot(122); hold off; 
    prim_non_degen = (sum(prim_rank)  >= m);
    dual_non_degen = (sum(dual_rank) <= m);
    fprintf('\n    sum(n2-(n-r)2)            = %2.0d (>=m)',sum(prim_rank)); 
    fprintf('\n    sum((n-s)2)               = %2.0d (<=m)',sum(dual_rank)); 
    fprintf('\n    nec. cond. prim_non_degen = %1d',prim_non_degen); 
    fprintf('\n    nec. cond. dual_non_degen = %1d',dual_non_degen); 
    fprintf('\n    strict_comp               = %1d\n',all(strict_comp)); 
    sig1 = svd(ZB1); 
    sig12 = svd([XB1, XB2]'); 
    fprintf('\n    svd(ZB1):          max, min = %2.1e %2.1e, cond = %2.1e\n',...
            max(sig1),min(sig1),max(sig1)/min(sig1)); 
    fprintf('    svd([XB1, XB2]^T): max, min = %2.1e %2.1e, cond = %2.1e\n',...
            max(sig12),min(sig12),max(sig12)/min(sig12)); 
%%***************************************************************

