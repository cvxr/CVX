%%*****************************************************************************
%% sqlpsummary: print summary
%%
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%*****************************************************************************

  function  sqlpsummary(info,ttime,infeas_org,printlevel) 

   iter        = info.iter; 
   obj         = info.obj; 
   gap         = info.gap; 
   relgap      = info.relgap;  
   prim_infeas = info.pinfeas;
   dual_infeas = info.dinfeas;  
   termcode    = info.termcode; 
   reldist     = info.reldist; 
   resid       = info.resid;
   dimacs      = info.dimacs; 
   totaltime   = info.cputime;
%%
   preproctime = ttime.preproc; pcholtime = ttime.pchol; dcholtime = ttime.dchol;
   predtime  = ttime.pred; pstep_predtime = ttime.pred_pstep; dstep_predtime = ttime.pred_pstep;
   corrtime  = ttime.corr; pstep_corrtime = ttime.corr_pstep; dstep_corrtime = ttime.corr_dstep;
   misctime  = ttime.misc; 
%%
   if (printlevel >= 2) 
      fprintf('\n------------------------------------------------');
      fprintf('-------------------\n');
      fprintf(' number of iterations   = %2.0f\n',iter);
   end
   if (termcode <= 0)
      if (printlevel >=2)
         fprintf(' primal objective value = %- 9.8e\n',obj(1));
         fprintf(' dual   objective value = %- 9.8e\n',obj(2));
         fprintf(' gap := trace(XZ)       = %3.2e\n',gap);
         fprintf(' relative gap           = %3.2e\n',relgap);
         fprintf(' actual relative gap    = %3.2e\n',-diff(obj)/(1+sum(abs(obj))));
         if norm(infeas_org)
            fprintf(' rel. primal infeas (scaled problem)   = %3.2e\n',prim_infeas);
            fprintf(' rel. dual     "        "       "      = %3.2e\n',dual_infeas);
            fprintf(' rel. primal infeas (unscaled problem) = %3.2e\n',infeas_org(1));
            fprintf(' rel. dual     "        "       "      = %3.2e\n',infeas_org(2));
         else
            fprintf(' rel. primal infeas     = %3.2e\n',prim_infeas);
            fprintf(' rel. dual   infeas     = %3.2e\n',dual_infeas);
         end
         fprintf(' norm(X), norm(y), norm(Z) = %3.1e, %3.1e, %3.1e\n',...
                   info.normX,info.normy,info.normZ);
         fprintf(' norm(A), norm(b), norm(C) = %3.1e, %3.1e, %3.1e\n',...
                   info.normA,info.normb,info.normC);
      end
   elseif (termcode == 1)
      if (printlevel >=2)
         fprintf(' residual of primal infeasibility      \n')
         fprintf(' certificate (y,Z)      = %3.2e\n',resid);
         fprintf(' reldist to infeas.    <= %3.2e\n',reldist);
      end
   elseif (termcode == 2)
      if (printlevel >=2)
         fprintf(' residual of dual infeasibility        \n')
         fprintf(' certificate X          = %3.2e\n',resid);
         fprintf(' reldist to infeas.    <= %3.2e\n',reldist);
      end
   end
   if (printlevel >=2)
      fprintf(' Total CPU time (secs)  = %3.2f  \n',totaltime);
      fprintf(' CPU time per iteration = %3.2f  \n',totaltime/iter);
      fprintf(' termination code       = %2.0f\n',termcode);
      fprintf(' DIMACS: %.1e  %.1e  %.1e  %.1e  %.1e  %.1e\n',dimacs); 
      fprintf('------------------------------------------------');
      fprintf('-------------------\n');
      if (printlevel > 3) 
         fprintf(' Percentage of CPU time spent in various parts \n'); 
         fprintf('------------------------------------------------');
         fprintf('-------------------\n');
         fprintf(' preproc Xchol Zchol pred pred_steplen  corr corr_steplen  misc\n')
         tt = [preproctime, pcholtime, dcholtime, predtime, pstep_predtime, dstep_predtime];
         tt = [tt, corrtime, pstep_corrtime, dstep_corrtime, misctime];
         tt = tt/sum(tt)*100; 
         fprintf('   %3.1f   %3.1f   %3.1f   %3.1f  %3.1f  %3.1f ',...
                 tt(1),tt(2),tt(3),tt(4),tt(5),tt(6));
         fprintf('    %3.1f  %3.1f  %3.1f     %3.1f\n',tt(7),tt(8),tt(9),tt(10)); 
         fprintf('------------------------------------------------');
         fprintf('-------------------\n');
      end
   end
%%*****************************************************************************
