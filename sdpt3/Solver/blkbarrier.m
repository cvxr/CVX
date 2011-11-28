%%********************************************************************
%% blkbarrier: calculate 
%% [-v(p)*logdet(X{p}),   v(p)*logdet(Z{p}) + n*v(p)*(1-log(v(p)))]
%% [-v(p)*log(gam(X{p})), v(p)*log(gam(Z{p})) + v(p)]
%% [-v(p)*log(X{p}),      v(p)*log(Z{p}) + n*v(p)*(1-log(v(p)))]
%%********************************************************************

  function objadd = blkbarrier(blk,X,Z,Xchol,Zchol,v); 

   objadd = zeros(1,2); tmp = zeros(1,2); 
   for p = 1:size(blk,1)
      pblk = blk(p,:);
      vp = v{p};
      idx = find(vp > 0);
      if ~isempty(idx) 
         vpsub = vp(idx); 
         if size(vpsub,1) < size(vpsub,2); vpsub = vpsub'; end
         if strcmp(pblk{1},'s')
            ss = [0, cumsum(pblk{2})]; 
            logdetX = 2*log(diag(Xchol{p})); 
            logdetZ = 2*log(diag(Zchol{p})); 
            logdetXsub = zeros(length(idx),1); 
            logdetZsub = zeros(length(idx),1); 
            for k = 1:length(idx)
               idxtmp = [ss(idx(k))+1:ss(idx(k)+1)]; 
               logdetXsub(k) = sum(logdetX(idxtmp)); 
               logdetZsub(k) = sum(logdetZ(idxtmp)); 
            end
            tmp(1) = -sum(vpsub.*logdetXsub); 
            tmp(2) = sum(vpsub.*logdetZsub + (pblk{2}(idx)').*vpsub.*(1-log(vpsub))); 
         elseif strcmp(pblk{1},'q')
            gamX = sqrt(qops(pblk,X{p},X{p},2)); 
            gamZ = sqrt(qops(pblk,Z{p},Z{p},2)); 
            tmp(1) = -sum(vpsub.*log(gamX(idx))); 
            tmp(2) = sum(vpsub.*log(gamZ(idx)) + vpsub); 
         elseif strcmp(pblk{1},'l')
            logX = log(X{p}); logZ = log(Z{p});
            tmp(1) = -sum(vpsub.*logX(idx)); 
            tmp(2) = sum(vpsub.*logZ(idx) + vpsub.*(1-log(vpsub))); 
         end 
         objadd = objadd + tmp; 
      end
   end 
%%********************************************************************
