%%**************************************************************************************
%% sqlpcheckconvg: check convergence.
%%
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%**************************************************************************************

   function [param,breakyes,restart,msg] = sqlpcheckconvg(param,runhist);

      termcode    = param.termcode; 
      iter        = param.iter; 
      obj         = param.obj;
      rel_gap     = param.rel_gap; 
      gap         = param.gap; 
      mu          = param.mu; 
      prim_infeas = param.prim_infeas;
      dual_infeas = param.dual_infeas;
      homRd       = param.homRd;
      homrp       = param.homrp;
      normX       = param.normX; 
      normZ       = param.normZ; 
      prim_infeas_bad = param.prim_infeas_bad; 
      prim_infeas_min = min(param.prim_infeas_min, max(prim_infeas,1e-9)); 
      normX0      = param.normX0; 
      normZ0      = param.normZ0;
      gaptol      = param.gaptol;
      inftol      = param.inftol;  
      printlevel  = param.printlevel;
      stoplevel   = param.stoplevel; 
      ublksize    = param.ublksize; 
      infeas_meas = max(prim_infeas,dual_infeas); 
      restart     = 0; 
      breakyes    = 0; 
      msg         = []; 
%%
      if (normX > 1e15*normX0 | normZ > 1e15*normZ0)
         termcode = 3;
         breakyes = 1; 
      end
      if (homRd < max(inftol,1e-13))
         termcode = 1;
         breakyes = 1;
      end
      if (homrp < max(inftol,1e-13))
         termcode = 2;
         breakyes = 1;
      end
      if (stoplevel) & (iter > 2)
         prim_infeas_bad = ... %%*** prim_infeas_bad (modified 25 Aug 2006)  
            + (prim_infeas > max(1e-10,1e2*prim_infeas_min) & (prim_infeas_min < 1e-2)) ...
            + (prim_infeas > prod(1.5-runhist.step(iter+1:iter-1))*runhist.pinfeas(iter-2));
         if (mu < 1e-8)
            idx = [max(1,iter-1): iter];
         elseif (mu < 1e-4);
            idx = [max(1,iter-2): iter]; 
         else
            idx = [max(1,iter-3): iter];
         end
         gap_ratio = runhist.gap(idx+1)./runhist.gap(idx); 
         idx2 = [max(1,iter-4): iter]; 
         gap_ratio2 = runhist.gap(idx2+1)./runhist.gap(idx2);
         gap_slowrate = min(0.8,max(0.6,2*mean(gap_ratio2)));
         idx3 = [max(1,iter-10): iter]; 
         gap_ratio3 = runhist.gap(idx3+1)./runhist.gap(idx3); 
         if (infeas_meas < 1e-4 | prim_infeas_bad) & (max(infeas_meas,rel_gap) < 1) ...
   	    & (iter > 20) & ~(min(runhist.step(idx)) > 0.2 & ublksize)
            gap_slow  = all(gap_ratio > gap_slowrate);  
            gap_slow3 = all(gap_ratio3 > gap_slowrate); 
            if ((gap_slow & (rel_gap < 5e-3)) | (gap_slow3 & prim_infeas_bad))
               msg = 'sqlp stop: progress is too slow'; 
               if (printlevel); fprintf('\n  %s',msg); end
               termcode = -5; 
               breakyes = 1;
            elseif (max(infeas_meas,rel_gap) < 1e-2) & (prim_infeas_bad) %% (modified 25 Aug 2006)
               if (rel_gap < max(0.2*prim_infeas,1e-2*dual_infeas)) 
                  msg = 'sqlp stop: relative gap < infeasibility'; 
                  if (printlevel); fprintf('\n  %s',msg); end
                  termcode = -1;
                  breakyes = 1; 
               end
            end
         end  
         if (infeas_meas < 1e-8) & (gap > 1.2*mean(runhist.gap(idx))) & (rel_gap < 5e-3) ...
	    & (prim_infeas_bad) %% (added 22 May 2007)
            msg = 'sqlp stop: progress is bad'; 
            if (printlevel); fprintf('\n  %s',msg); end
            termcode = -5;
            breakyes = 1;  
         end
	 if (prim_infeas_bad) & (iter > 30) & all(gap_ratio > gap_slowrate) ...
            & (dual_infeas < 1e-3) & ~(min(runhist.step(idx)) > 0.2 & ublksize)
            msg = 'sqlp stop: progress is bad*'; 
            if (printlevel); fprintf('\n  %s',msg); end
            termcode = -5;
            breakyes = 1; 
         end
         if (min(runhist.infeas) < 1e-4 | (prim_infeas_bad & iter > 10)) ...
            & (max(runhist.infeas) > 1e-5) | (iter > 20) 
            rel_gap2 = abs(diff(obj))/(1+sum(abs(obj))); 
            if (rel_gap2 < 1e-3); 
               step_short = all(runhist.step([iter:iter+1]) < 0.05) ;
            elseif (rel_gap2 < 1) 
               idx = [max(1,iter-3): iter+1];
               step_short = all(runhist.step(idx) < 0.03); 
	    else
               step_short = 0; 
            end
            if (step_short) 
               msg = 'sqlp stop: steps too short consecutively'; 
               if (printlevel); fprintf('\n  %s',msg); end
               termcode = -5; 
               breakyes = 1;      
            end
         end
         if (iter > 3 & iter < 20) & (max(runhist.step(max(1,iter-3):iter+1)) < 1e-3) ...
            & (infeas_meas > 1) & (min(homrp,homRd) > min(1e-8,inftol)) 
            if (stoplevel == 2) 
               msg = 'sqlp stop: steps too short consecutively'; 
               if (printlevel)
                  fprintf('\n *** Too many tiny steps, advisable to restart sqlp'); 
                  fprintf(' with the following iterate.')
                  fprintf('\n *** Suggestion: [X0,y0,Z0] = infeaspt(blk,At,C,b,2,1e5);'); 
                  fprintf('\n  %s',msg); 
               end
               termcode = -5; 
               breakyes = 1;             
            elseif (stoplevel == 3)
               msg = 'sqlp stop: steps too short consecutively'; 
               if (printlevel)
                  fprintf('\n *** Too many tiny steps even')
                  fprintf(' after restarting sqlp');                
                  fprintf('\n  %s',msg);
               end
               termcode = -5;
               breakyes = 1;              
            else 
               if (printlevel)
                  fprintf('\n *** Too many tiny steps:')
                  fprintf('  restarting sqlp with the following iterate.')
                  fprintf('\n *** [X,y,Z] = infeaspt(blk,At,C,b,2,1e5);'); 
               end
	       prim_infeas_min = 1e20; 
               prim_infeas_bad = 0; 
   	       restart = 1; 
            end
         end
      end
      if (max(rel_gap,infeas_meas) < gaptol)
         msg = sprintf('sqlp stop: max(relative gap, infeasibilities) < %3.2e',gaptol); 
         if (printlevel); fprintf('\n  %s',msg); end
         termcode = 0;
         breakyes = 1;
      end
%%
      param.prim_infeas_bad = prim_infeas_bad;
      param.prim_infeas_min = prim_infeas_min;
      param.termcode = termcode;
%%**************************************************************************************
