%%*****************************************************************************
%% HSDsqlpcheckconvg: check convergence.
%%
%% ZpATynorm, AX, normX, normZ are with respect to the 
%% original variables, not the HSD variables. 
%%
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%*****************************************************************************

  function [param,breakyes,use_olditer,msg] = HSDsqlpcheckconvg(param,runhist);

      termcode    = param.termcode; 
      iter        = param.iter; 
      obj         = param.obj;
      relgap      = param.relgap; 
      gap         = param.gap; 
      prim_infeas = param.prim_infeas;
      dual_infeas = param.dual_infeas; 
      mu          = param.mu;  
      prim_infeas_bad = param.prim_infeas_bad; 
      dual_infeas_bad = param.dual_infeas_bad; 
      printlevel  = param.printlevel;
      stoplevel   = param.stoplevel; 
      inftol      = param.inftol; 
      gaptol      = param.gaptol;
      kap         = param.kap; 
      tau         = param.tau;       
      theta       = param.theta; 
      breakyes    = 0; 
      use_olditer = 0; 
      msg         = [];
      infeas      = max(prim_infeas,dual_infeas); 
      prim_infeas_min = min(param.prim_infeas_min, max(prim_infeas,1e-10));  
      dual_infeas_min = min(param.dual_infeas_min, max(dual_infeas,1e-10)); 
%%
      err = max(infeas,relgap); 
      if (obj(2) > 0); homRd = param.ZpATynorm/obj(2); else; homRd = inf; end
      if (obj(1) < 0); homrp = norm(param.AX)/(-obj(1)); else; homrp = inf; end
      if (param.normX > 1e15*param.normX0 | param.normZ > 1e15*param.normZ0)
         termcode = 3;
         breakyes = 1; 
      end
      if (homRd < min(1e-6,1e-2*sqrt(err*inftol)) & tau < 1e-4 ...
          & prim_infeas > 0.5*runhist.pinfeas(iter)) ...
         | (homRd < 10*tau & tau < 1e-7)
         termcode = 1;
         breakyes = 1;
      end
      if (homrp < min(1e-6,1e-2*sqrt(err*inftol)) & tau < 1e-4 ...
         & dual_infeas > 0.5*runhist.dinfeas(iter)) ...
         | (homrp < 10*tau & tau < 1e-7)
         termcode = 2;
         breakyes = 1;
      end
      if (err < gaptol)
         msg = sprintf('Stop: max(relative gap,infeasibilities) < %3.2e',gaptol);
         if (printlevel); fprintf('\n  %s',msg); end
         termcode = 0;
         breakyes = 1;
      end
      min_prim_infeas = min(runhist.pinfeas(1:iter)); 
      prim_infeas_bad = prim_infeas_bad + (prim_infeas > ...
           max(1e-10,5*min_prim_infeas) & (min_prim_infeas < 1e-2));
      if (mu < 1e-6)
         idx = [max(1,iter-1): iter];
      elseif (mu < 1e-3);
         idx = [max(1,iter-2): iter]; 
      else
         idx = [max(1,iter-3): iter];
      end
      idx2 = [max(1,iter-4): iter]; 
      gap_ratio2 = runhist.gap(idx2+1)./runhist.gap(idx2);
      gap_slowrate = min(0.8,max(0.6,2*mean(gap_ratio2)));
      gap_ratio = runhist.gap(idx+1)./runhist.gap(idx);
      pstep = runhist.step(iter+1);  
      if (infeas < 1e-4 | prim_infeas_bad) & (relgap < 1e-3) ...
         & (iter > 5) & (prim_infeas > (1-pstep/2)*runhist.pinfeas(iter)) 
         gap_slow = all(gap_ratio > gap_slowrate) & (relgap < 1e-3);
         min_pinfeas = min(runhist.pinfeas); 
         if (relgap < 0.1*infeas) ...
	    & ((runhist.step(iter+1) < 0.5) | (min_pinfeas < min(1e-6,0.1*prim_infeas))) ...
            & (dual_infeas > 0.9*runhist.dinfeas(iter) | (dual_infeas < 1e-2*gaptol))
            msg = 'Stop: relative gap < infeasibility'; 
            if (printlevel); fprintf('\n  %s',msg); end
            termcode = -1;
            breakyes = 1;           
         elseif (gap_slow) & (infeas > 0.9*runhist.infeas(iter)) ...
            & (theta < 1e-8)
            msg = 'Stop: progress is too slow'; 
            if (printlevel); fprintf('\n  %s',msg); end
            termcode = -5; 
            breakyes = 1;
         end  
      elseif (prim_infeas_bad) & (iter >50) & all(gap_ratio > gap_slowrate)
         msg = 'Stop: progress is bad';
         if (printlevel); fprintf('\n  %s',msg); end
         termcode = -5;
         breakyes = 1; 
      elseif (infeas < 1e-8) & (gap > 1.2*mean(runhist.gap(idx)))
         msg = 'Stop: progress is bad*'; 
         if (printlevel); fprintf('\n  %s',msg); end
         termcode = -5;
         breakyes = 1;  
      end
      if (err < 1e-3) & (iter > 10) ...
         & (runhist.pinfeas(iter+1) > 0.9*runhist.pinfeas(max(1,iter-5))) ...
         & (runhist.dinfeas(iter+1) > 0.9*runhist.dinfeas(max(1,iter-5))) ...
         & (runhist.relgap(iter+1)  > 0.1*runhist.relgap(max(1,iter-5))); 
         msg = 'Stop: progress is bad**';
         if (printlevel); fprintf('\n  %s',msg); end
         termcode = -5;
         breakyes = 1;  
      end
      if (infeas > 100*max(1e-12,min(runhist.infeas)) & relgap < 1e-4)
         msg = 'Stop: infeas has deteriorated too much'; 
         if (printlevel); fprintf('\n  %s, %3.1e',msg,infeas); end
         use_olditer = 1; 
         termcode = -7; 
         breakyes = 1; 
      end
      if (min(runhist.infeas) < 1e-4 | prim_infeas_bad) ...
         & (max(runhist.infeas) > 1e-4) & (iter > 5)
         relgap2 = abs(diff(obj))/(1+mean(abs(obj))); 
         if (relgap2 < 1e-3); 
            step_short = all(runhist.step([iter:iter+1]) < 0.1) ;
         elseif (relgap2 < 1) 
            idx = [max(1,iter-3): iter+1];
            step_short = all(runhist.step(idx) < 0.05); 
         else
            step_short = 0; 
         end
         if (step_short) 
            msg = 'Stop: steps too short consecutively'; 
            if (printlevel); fprintf('\n  %s',msg); end
            termcode = -5; 
            breakyes = 1;      
         end
      end
      if (iter > 3 & iter < 20) & (max(runhist.step(max(1,iter-3):iter+1)) < 1e-3) ...
         & (infeas > 1) & (min(homrp,homRd) > 1000*inftol) 
         if (stoplevel >= 2)
            msg = 'Stop: steps too short consecutively'; 
            if (printlevel); fprintf('\n  %s',msg); end
            termcode = -5;
            breakyes = 1;              
         end
      end
      if (pstep < 1e-4) & (err > 1.1*max(runhist.relgap(iter),runhist.infeas(iter)))
         msg = 'Stop: steps are too short';
         if (printlevel); fprintf('\n  %s',msg); end
         use_olditer = 1; 
         termcode = -5; 
         breakyes = 1; 
      end
      if (iter == param.maxit) 
         termcode = -6; 
         msg = 'Stop: maximum number of iterations reached'; 
         if (printlevel); fprintf('\n  %s',msg); end
      end
      if (infeas < 1e-8 & relgap < 1e-10 & kap < 1e-13 & theta < 1e-15)
         msg = 'Stop: obtained accurate solution';
         if (printlevel); fprintf('\n  %s',msg); end
         termcode = 0; 
         breakyes = 1;          
      end
      param.prim_infeas_bad = prim_infeas_bad;
      param.prim_infeas_min = prim_infeas_min;
      param.dual_infeas_bad = dual_infeas_bad;
      param.dual_infeas_min = dual_infeas_min;
      param.termcode = termcode;
%%*****************************************************************************
