%%**********************************************************************
%% NTscaling: Compute NT scaling matrix
%%                       
%% compute SVD of Xchol*Zchol via eigenvalue decompostion of
%%     Zchol * X * Zchol' = V * diag(sv2) * V'. 
%% compute W satisfying W*Z*W = X. 
%%     W = G'*G,  where G = diag(sqrt(sv)) * (invZchol*V)'
%%     important to keep W symmertic.
%%
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%**********************************************************************

 function [W,G,sv,gamx,gamz,dd,ee,ff] = ...
          NTscaling(blk,X,Z,Zchol,invZchol);

    numblk = size(blk,1);
    W = cell(numblk,1); G = cell(numblk,1); sv = cell(numblk,1); 
    gamx = cell(numblk,1); gamz = cell(numblk,1); 
    dd = cell(numblk,1); ee = cell(numblk,1); ff = cell(numblk,1);   
%%
    for p = 1:size(blk,1)
        pblk = blk(p,:); 
        numblk = length(pblk{2});  
        n = sum(pblk{2});  
        if strcmp(pblk{1},'l')
           dd{p} = X{p}./Z{p}; %% do not add perturbation, it badly affects cre-a   
        elseif strcmp(pblk{1},'q');  
           gamx{p} = sqrt(qops(pblk,X{p},X{p},2)); 
           gamz{p} = sqrt(qops(pblk,Z{p},Z{p},2)); 
           w2 = gamz{p}./gamx{p};  w = sqrt(w2); 
           dd{p} = qops(pblk,1./w2,ones(n,1),4);
           tt = qops(pblk,1./w,Z{p},3) - qops(pblk,w,X{p},4);
           gamtt = sqrt(qops(pblk,tt,tt,2)); 
           ff{p} = qops(pblk,1./gamtt,tt,3); 
           ee{p} = qops(pblk,sqrt(2)./w,ff{p},4); 
        elseif strcmp(pblk{1},'s')   
           tmp = Prod2(pblk,Zchol{p},X{p},0); 
           tmp = Prod2(pblk,tmp,Zchol{p}',1); 
           [sv2,V] = blkeig(pblk,tmp); 
           sv2 = max(1e-20,sv2); 
           sv{p} = sqrt(sv2); 
           tmp  = Prod2(pblk,invZchol{p},V); 
           G{p} = Prod2(pblk,spdiags(sqrt(sv{p}),0,n,n),tmp'); 
           W{p} = Prod2(pblk,G{p}',G{p},1);                    
        end
    end
%%**********************************************************************



