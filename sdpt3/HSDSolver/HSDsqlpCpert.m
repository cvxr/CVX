%%*****************************************************************************
%% HSDsqlpCpert: perturb C. 
%%
%%
%%*****************************************************************************
  
    function [At,Cpert] = HSDsqlpCpert(blk,At,par,C,X,Cpert,runhist); 

    iter = length(runhist.pinfeas); 
    prim_infeas = runhist.pinfeas(iter); 
    dual_infeas = runhist.dinfeas(iter); 
    relgap      = runhist.relgap(iter); 
    infeas      = runhist.infeas(iter); 
    theta       = runhist.theta(iter); 
%%
    Cpertold = Cpert; 
    err = max(relgap,infeas); 
    for p = 1:size(blk,1) 
       pblk = blk(p,:); 
       n = sum(pblk{2}); 
       tmp = max(1,norm(C{p},'fro'))/sqrt(n);
       if (err < 1e-6)
          if (norm(X{p},'fro') < 1e2); const=0.2; else; const=0.3; end 
          Cpert(p) = max(const*Cpert(p),1e-10*tmp);  
       elseif (err < 1e-2) 
          if (norm(X{p},'fro') < 1e2); const=0.4; else; const=0.5; end 
          Cpert(p) = max(const*Cpert(p),1e-8*tmp); 
       else
	  Cpert(p) = max(0.9*Cpert(p),1e-6*tmp); 
       end
       Cpert = min(Cpert,Cpertold); 
       if (prim_infeas < min([0.1*dual_infeas, 1e-7*runhist.pinfeas(1)])) ...
      	  & (iter > 1 & dual_infeas > 0.8*runhist.dinfeas(iter-1) & relgap < 1e-4)
          Cpert(p) = 0.5*Cpert(p);        
       elseif (dual_infeas < min([0.1*prim_infeas, 1e-7*runhist.dinfeas(1)])) ...
          & (iter > 1 & prim_infeas > 0.8*runhist.pinfeas(iter-1) & relgap < 1e-4)
          Cpert(p) = 0.5*Cpert(p);       
       elseif (max(relgap,1e-2*infeas) < 1e-6 & relgap < 0.1*infeas) 
          Cpert(p) = 0.5*Cpert(p);
       end
       if (prim_infeas < min([1e-4*dual_infeas,1e-7]) & theta < 1e-6) ...
	  | (prim_infeas < 1e-4 & theta < 1e-10) 
          Cpert(p) = 0.1*Cpert(p);  
       elseif (dual_infeas < min([1e-4*prim_infeas,1e-7]) & theta < 1e-6) ...
	  | (dual_infeas < 1e-4 & theta < 1e-10) 
          Cpert(p) = 0.1*Cpert(p);         
       elseif (iter > 1 & theta > 0.9*runhist.theta(iter-1) & infeas < 1e-3)
          Cpert(p) = 0.1*Cpert(p); 
       end
       if strcmp(pblk{1},'s')
          Cnew = C{p} + Cpert(p)*speye(n); 
          At{p}(:,par.invpermA(p,end-1)) = -svec(pblk,Cnew,1); 
       else
          Cnew = C{p} + Cpert(p)*ones(n,1); 
          At{p}(:,par.invpermA(p,end-1)) = -Cnew; 
       end
    end
%%*****************************************************************************
