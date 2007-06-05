%%*************************************************************************
%% sqlp: main solver 
%%
%%*************************************************************************
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%*************************************************************************

  function [obj,X,y,Z,info,runhist] = sqlpmain(blk,At,C,b,par,parbarrier,X0,y0,Z0);

   global spdensity printlevel msg
   global solve_ok  exist_analytic_term
   global schurfun  schurfun_par 
 %%
   matlabversion = par.matlabversion;
   vers          = par.vers;
   predcorr      = par.predcorr;
   gam           = par.gam; 
   expon         = par.expon;
   gaptol        = par.gaptol;
   inftol        = par.inftol;
   steptol       = par.steptol;
   maxit         = par.maxit;
   printlevel    = par.printlevel;
   stoplevel     = par.stoplevel;
   scale_data    = par.scale_data;
   spdensity     = par.spdensity;
   rmdepconstr   = par.rmdepconstr;
   cachesize     = par.cachesize; 
   smallblkdim   = par.smallblkdim;
   schurfun      = par.schurfun;
   schurfun_par  = par.schurfun_par;
   ublksize      = par.ublksize; 
%%
   tstart = cputime; 
   X = X0; y = y0; Z = Z0; 
   for p = 1:size(blk,1)
      if strcmp(blk{p,1},'u'); Z{p} = zeros(blk{p,2},1); end
   end
%%
%%-----------------------------------------
%% convert unrestricted blk to linear blk. 
%%-----------------------------------------
%%
   randstate = rand('state');
   rand('state',0);
%%
   ublkidx = zeros(size(blk,1),1); 
   for p = 1:size(blk,1) 
      pblk = blk(p,:); 
      if strcmp(pblk{1},'u') & (pblk{2} > 0) 
         convert2lblk = 1; 
         if (sum(pblk{2}) < min(length(b),20))
            AAt = At{p}*At{p}';       
            if (matlabversion < 7.3)
               [Lsymb,flag] = symbcholfun(AAt,cachesize);
               if (flag==0); 
                  L = sparcholfun(Lsymb,AAt);  
                  if ~any(L.skip) & (max(L.d)/min(L.d) < 1e6)
                     convert2lblk = 0; 
                     if (printlevel); fprintf(' *** no conversion for ublk'); end
                  end
               end
            else
              [L.R,L.p,L.perm] = chol(AAt,'vector'); 
              L.d = full(diag(L.R)).^2; 
              L.skip = []; 
              if ~any(L.skip) & (max(L.d)/min(L.d) < 1e6)
                 convert2lblk = 0; 
                 if (printlevel); fprintf(' *** no conversion for ublk'); end
              end
            end
         end
         if (convert2lblk) 
            ublkidx(p) = 1; 
            n = 2*blk{p,2}; 
            blk{p,1} = 'l'; 
            blk{p,2} = n;
            parbarrier{p} = [parbarrier{p}, parbarrier{p}];
            At{p} = [At{p}; -At{p}];       
            C{p} = [C{p}; -C{p}]; 
            b2 = 1 + abs(b');  
            normCtmp = 1+norm(C{p});
            normAtmp = 1+sqrt(sum(At{p}.*At{p}));
            if (n > 1000)
               const = sqrt(n); 
            else
	       const = n; 
            end
            if (par.startpoint == 1)
               X{p} = const* max([1,b2./normAtmp]) *ones(n,1); 
               Z{p} = const* max([1,normAtmp/sqrt(n),normCtmp/sqrt(n)]) *ones(n,1);
               X{p} = X{p}.*(1+1e-10*rand(n,1)); 
               Z{p} = Z{p}.*(1+1e-10*rand(n,1)); 
	    else
               const = max(abs(X{p})) + 100; 
               X{p} = [X{p}+const; const*ones(n/2,1)]; 
               const = 100; 
               Z{p} = [const*ones(n/2,1); const*ones(n/2,1)]; 
            end
         end
      end
   end 
   rand('state',randstate);
%%
%%-----------------------------------------
%% check whether {A1,...,Am} is 
%% linearly independent. 
%%-----------------------------------------
%%
   m0 = length(b); 
   [At,b,y,indeprows,par.depconstr,feasible,par.AAt] = ...
    checkdepconstr(blk,At,b,y,rmdepconstr);
   if (~feasible)
      obj = []; X = cell(size(blk,1),1); y = []; Z = cell(size(blk,1),1); 
      runhist = [];      
      msg = 'SQLP is not feasible'; 
      if (printlevel); fprintf('\n %s \n',msg); end
      return;
   end
   par.normAAt = norm(par.AAt,'fro'); 
%%
%%-----------------------------------------
%% scale SQLP data. Note: must be done only 
%% after checkdepconstr
%%-----------------------------------------
%%
   normb2 = 1+norm(b); 
   normC2 = 1+ops(C,'norm');
   normA2 = 1+ops(At,'norm'); 
   normX0 = 1+ops(X0,'norm'); 
   normZ0 = 1+ops(Z0,'norm'); 
   if (scale_data)
      [At,C,b,normA,normC,normb,X,y,Z] = scaling(blk,At,C,b,X,y,Z);
   else
      normA = 1; normC = 1; normb = 1; 
   end 
%%
%%-----------------------------------------
%% find the combined list of non-zero 
%% elements of Aj, j = 1:k, for each k. 
%% IMPORTANT NOTE: Ak, C are permuted.
%%-----------------------------------------
%% 
   par.numcolAt = length(b); 
   [At,C,X,Z,par.permA,par.permZ] = sortA(blk,At,C,b,X,Z);
   [par.isspA,par.nzlistA,par.nzlistAsum,par.isspAy,par.nzlistAy] = nzlist(blk,At,par);
%%
%%-----------------------------------------
%% create an artifical non-negative block 
%% for a purely log-barrier problem
%%-----------------------------------------
%%
   numblkold = size(blk,1);  
   nn = 0; 
   for p = 1:size(blk,1);
      pblk = blk(p,:);  
      idx = find(parbarrier{p}==0); 
      if ~isempty(idx); 
         if strcmp(pblk{1},'l') 
            nn = nn + length(idx); 
         elseif strcmp(pblk{1},'s') | strcmp(pblk{1},'q')   
            nn = nn + sum(pblk{2}(idx)); 
         end
      end
   end
   if (nn==0)
      analytic_prob = 1; 
      numblk = size(blk,1)+1; 
      blk{numblk,1} = 'l'; blk{numblk,2} = 1; 
      At{numblk,1} = sparse(1,length(b)); 
      C{numblk,1} = 1; 
      X{numblk,1} = 1e3; 
      Z{numblk,1} = 1e3;
      parbarrier{numblk,1} = 0; 
      ublkidx(numblk,1) = 0;
      nn = nn + 1; 
   else
      analytic_prob = 0;       
   end
%%
   exist_analytic_term = 0; 
   for p = 1:size(blk,1);
      idx = find(parbarrier{p} > 0); 
      if ~isempty(idx); 
         exist_analytic_term = 1; 
      end
   end
%%
%%-----------------------------------------
%% initialization
%%-----------------------------------------
%%
%%
   [Xchol,indef(1)] = blkcholfun(blk,X); 
   [Zchol,indef(2)] = blkcholfun(blk,Z); 
   if any(indef)
      msg = 'sqlp stop: X, Z are not both positive definite';
      if (printlevel); fprintf('\n  %s\n',msg); end
      info.termcode = -3;
      info.msg1 = msg;
      obj = []; X = cell(size(blk,1),1); y = []; Z = cell(size(blk,1),1); 
      runhist = [];      
      return;
   end 
   AX = AXfun(blk,At,par.permA,X); 
   rp = b-AX;
   ZpATy = ops(Z,'+',Atyfun(blk,At,par.permA,par.isspAy,y));
   ZpATynorm = ops(ZpATy,'norm');
   Rd = ops(C,'-',ZpATy);
   objadd0 = 0; 
   if (scale_data)
      for p = 1:size(blk,1)
         pblk = blk(p,:); 
         objadd0 = objadd0 + sum(parbarrier{p}.*pblk{2})*log(normA{p}); 
      end
   end
   objadd = blkbarrier(blk,X,Z,Xchol,Zchol,parbarrier) + objadd0;
   obj = (normb*normC)*[blktrace(blk,C,X), b'*y] + objadd;      
   gap = (normb*normC)*blktrace(blk,X,Z) - diff(objadd); 
   rel_gap = gap/(1+sum(abs(obj)));
   prim_infeas = norm(rp)/normb2;
   dual_infeas = ops(Rd,'norm')/normC2;
   infeas_meas = max(prim_infeas,dual_infeas); 
   if (scale_data)
      infeas_org(1) = prim_infeas*normb;
      infeas_org(2) = dual_infeas*normC;
   else
      infeas_org = [0,0]; 
   end
   trXZ = blktrace(blk,X,Z,parbarrier); 
   if (nn > 0); mu  = trXZ/nn; else; mu = gap/ops(X,'getM'); end
%%   
   termcode = 0; 
   pstep = 0; dstep = 0; pred_convg_rate = 1; corr_convg_rate = 1;
   prim_infeas_bad = 0;  prim_infeas_min = prim_infeas; homRd = inf; homrp = inf; 
   msg = []; msg2 = []; msg3 = [];
   runhist.pobj = obj(1);
   runhist.dobj = obj(2); 
   runhist.gap  = gap;
   runhist.relgap  = rel_gap;
   runhist.pinfeas = prim_infeas;
   runhist.dinfeas = dual_infeas;
   runhist.infeas  = infeas_meas;  
   runhist.step    = 0; 
   runhist.cputime = cputime-tstart; 
   ttime.preproc   = runhist.cputime; 
   ttime.pred = 0; ttime.pred_pstep = 0; ttime.pred_dstep = 0; 
   ttime.corr = 0; ttime.corr_pstep = 0; ttime.corr_dstep = 0; 
   ttime.pchol = 0; ttime.dchol = 0; ttime.misc = 0; 
%%
%%-----------------------------------------
%% display parameters and initial info
%%-----------------------------------------
%%
   if (printlevel >= 2)
      fprintf('\n********************************************');
      fprintf('***********************\n');
      fprintf('   SDPT3: Infeasible path-following algorithms'); 
      fprintf('\n********************************************');
      fprintf('***********************\n');
      [hh,mm,ss] = mytime(ttime.preproc); 
      if (printlevel>=3)       
         fprintf(' version  predcorr  gam  expon  scale_data\n');
         if (vers == 1); fprintf('   HKM '); elseif (vers == 2); fprintf('    NT '); end
         fprintf('     %1.0f      %4.3f',predcorr,gam);
         fprintf('   %1.0f        %1.0f    %1.0f\n',expon,scale_data); 
         fprintf('\nit  pstep dstep p_infeas d_infeas  gap')
         fprintf('     mean(obj)    cputime\n');
         fprintf('------------------------------------------------');
         fprintf('-------------------\n');
         fprintf('%2.0f  %4.3f %4.3f %2.1e %2.1e',0,0,0,prim_infeas,dual_infeas);
         fprintf('  %2.1e %- 7.6e  %s:%s:%s',gap,mean(obj),hh,mm,ss);
      end
   end
%%
%%---------------------------------------------------------------
%% start main loop
%%---------------------------------------------------------------
%%
   param.termcode    = termcode; 
   param.iter        = 0; 
   param.normA       = normA; 
   param.normb       = normb;
   param.normC       = normC;
   param.normX0      = normX0; 
   param.normZ0      = normZ0; 
   param.m0          = m0;
   param.indeprows   = indeprows;
   param.prim_infeas_bad = prim_infeas_bad; 
   param.prim_infeas_min = prim_infeas_min; 
   param.gaptol      = gaptol;
   param.inftol      = inftol; 
   param.maxit       = maxit;
   param.scale_data  = scale_data;
   param.printlevel  = printlevel; 
   param.ublksize    = ublksize; 
%%
   for iter = 1:maxit;  

       tstart  = cputime;  
       timeold = cputime;
       update_iter = 0; breakyes = 0; pred_slow = 0; corr_slow = 0; step_short = 0; 
       par.parbarrier = parbarrier; 
       par.iter = iter; 
       par.obj  = obj; 
       par.y    = y; 
%%
%%---------------------------------------------------------------
%% predictor step.
%%---------------------------------------------------------------
%%
       if (predcorr)
          sigma = 0.05; 
       else 
          sigma = 1-0.9*min(pstep,dstep); 
          if (iter == 1); sigma = 0.5; end; 
       end
       sigmu = cell(size(blk,1),1);
       for p = 1:size(blk,1)
          sigmu{p} = max(sigma*mu, parbarrier{p}');  
       end
       invXchol = cell(size(blk,1),1); 
       invZchol = ops(Zchol,'inv'); 
       if (vers == 1);
          [par,dX,dy,dZ,coeff,L,hRd] = ...
           HKMpred(blk,At,par,rp,Rd,sigmu,X,Z,invZchol);
       elseif (vers == 2);
          [par,dX,dy,dZ,coeff,L,hRd] = ...
           NTpred(blk,At,par,rp,Rd,sigmu,X,Z,Zchol,invZchol);
       end
       if (solve_ok <= 0)
          msg = 'sqlp stop: difficulty in computing predictor directions'; 
          if (printlevel); fprintf('\n  %s',msg); end
          runhist.cputime(iter+1) = cputime-tstart; 
          termcode = -4;
          break; %% do not ues breakyes = 1
       end
       timenew = cputime;
       ttime.pred = ttime.pred + timenew-timeold; timeold = timenew; 
%%
%%-----------------------------------------
%% step-lengths for predictor step
%%-----------------------------------------
%%
      if (gam == 0) 
         gamused = 0.9 + 0.09*min(pstep,dstep); 
      else
         gamused = gam;
      end 
      [Xstep,invXchol] = steplength(blk,X,dX,Xchol,invXchol); 
      pstep = min(1,gamused*full(Xstep));
      if (Xstep > .99e12) & (blktrace(blk,C,dX) < -1e-3) & (prim_infeas < 1e-3)
         pstep = Xstep; 
         msg = 'Predictor: dual seems infeasible'; 
         if (printlevel); fprintf('\n %s',msg); end
      end
      timenew = cputime; 
      ttime.pred_pstep = ttime.pred_pstep + timenew-timeold; timeold = timenew;
      Zstep = steplength(blk,Z,dZ,Zchol,invZchol); 
      dstep = min(1,gamused*full(Zstep));
      if (Zstep > .99e12) & (b'*dy > 1e-3) & (dual_infeas < 1e-3)
         dstep = Zstep; 
         msg = 'Predictor: primal seems infeasible'; 
         if (printlevel); fprintf('\n %s',msg); end
      end
      trXZpred = trXZ + pstep*blktrace(blk,dX,Z,parbarrier) ...
                 + dstep*blktrace(blk,X,dZ,parbarrier) ...
                 + pstep*dstep*blktrace(blk,dX,dZ,parbarrier);
      if (nn > 0); mupred  = trXZpred/nn; else; mupred = 1e-16; end
      mupredhist(iter) = mupred; 
      timenew = cputime;        
      ttime.pred_dstep = ttime.pred_dstep + timenew-timeold; timeold = timenew;   
%%
%%-----------------------------------------
%%  stopping criteria for predictor step.
%%-----------------------------------------
%%
      if (min(pstep,dstep) < steptol) & (stoplevel) & (iter > 10)
         msg = 'sqlp stop: steps in predictor too short';
         if (printlevel) 
            fprintf('\n  %s',msg);
            fprintf(': pstep = %3.2e,  dstep = %3.2e\n',pstep,dstep);
         end
         runhist.cputime(iter+1) = cputime-tstart; 
         termcode = -2; 
         breakyes = 1; 
      end
      if (iter >= 2) 
         idx = [max(2,iter-2) : iter];
         pred_slow = all(mupredhist(idx)./mupredhist(idx-1) > 0.4);
         idx = [max(2,iter-5) : iter];
         pred_convg_rate = mean(mupredhist(idx)./mupredhist(idx-1));
         pred_slow = pred_slow + (mupred/mu > 5*pred_convg_rate);
      end 
      if (~predcorr)
         if (max(mu,infeas_meas) < 1e-6) & (pred_slow) & (stoplevel)
            msg = 'sqlp stop: lack of progress in predictor'; 
            if (printlevel) 
               fprintf('\n  %s',msg);
               fprintf(': mupred/mu = %3.2f, pred_convg_rate = %3.2f.',...
               mupred/mu,pred_convg_rate);
            end
            runhist.cputime(iter+1) = cputime-tstart; 
            termcode = -2; 
            breakyes = 1;
         else 
            update_iter = 1; 
         end
      end
%%
%%---------------------------------------------------------------
%% corrector step.
%%---------------------------------------------------------------
%%
      if (predcorr) & (~breakyes)
         step_pred = min(pstep,dstep);
         if (mu > 1e-6)
            if (step_pred < 1/sqrt(3)); 
               expon_used = 1; 
            else
               expon_used = max(expon,3*step_pred^2); 
            end
         else 
            expon_used = max(1,min(expon,3*step_pred^2)); 
         end 
         if (nn > 0); 
            sigma = min( 1, max(0.05,(mupred/mu)^expon_used) );
         else
            sigma = 0.2; 
         end
         sigmu = cell(size(blk,1),1); 
         for p = 1:size(blk,1)
            sigmu{p} = max(sigma*mu, parbarrier{p}'); 
         end
%%
         if (vers == 1)
            [dX,dy,dZ] = HKMcorr(blk,At,par,rp,Rd,sigmu,hRd,...
             dX,dZ,coeff,L,X,Z);
         elseif (vers == 2)
            [dX,dy,dZ] = NTcorr(blk,At,par,rp,Rd,sigmu,hRd,...
             dX,dZ,coeff,L,X,Z); 
         end
         if (solve_ok <= 0)
            msg = 'sqlp stop: difficulty in computing corrector directions'; 
            if (printlevel); fprintf('\n  %s',msg); end
            runhist.cputime(iter+1) = cputime-tstart; 
            termcode = -4;
            break; %% do not ues breakyes = 1
         end
         timenew = cputime;
         ttime.corr = ttime.corr + timenew-timeold; timeold = timenew; 
%%
%%-----------------------------------
%% step-lengths for corrector step
%%-----------------------------------
%%
         if (gam == 0) 
            gamused = 0.9 + 0.09*min(pstep,dstep); 
         else
            gamused = gam;
         end            
         Xstep = steplength(blk,X,dX,Xchol,invXchol);
         pstep = min(1,gamused*full(Xstep));
         if (Xstep > .99e12) & (blktrace(blk,C,dX) < -1e-3) & (prim_infeas < 1e-3)
            pstep = Xstep;
            msg = 'Corrector: dual seems infeasible'; 
            if (printlevel); fprintf('\n %s',msg); end
         end
         timenew = cputime;
         ttime.corr_pstep = ttime.corr_pstep + timenew-timeold; timeold = timenew;
         Zstep = steplength(blk,Z,dZ,Zchol,invZchol);
         dstep = min(1,gamused*full(Zstep));
         if (Zstep > .99e12) & (b'*dy > 1e-3) & (dual_infeas < 1e-3)
            dstep = Zstep;
            msg = 'Corrector: primal seems infeasible'; 
            if (printlevel); fprintf('\n %s',msg); end
         end     
         trXZcorr = trXZ + pstep*blktrace(blk,dX,Z,parbarrier) ...
                    + dstep*blktrace(blk,X,dZ,parbarrier)...
                    + pstep*dstep*blktrace(blk,dX,dZ,parbarrier); 
         if (nn > 0); mucorr  = trXZcorr/nn; else; mucorr = 1e-16; end
         timenew = cputime;
         ttime.corr_dstep = ttime.corr_dstep + timenew-timeold; timeold = timenew; 
%%
%%-----------------------------------------
%%  stopping criteria for corrector step
%%-----------------------------------------
%%
         if (iter >= 2) 
            idx = [max(2,iter-2) : iter];
            corr_slow = all(runhist.gap(idx)./runhist.gap(idx-1) > 0.8); 
            idx = [max(2,iter-5) : iter];
            corr_convg_rate = mean(runhist.gap(idx)./runhist.gap(idx-1));
            corr_slow = corr_slow + (mucorr/mu > max(min(1,5*corr_convg_rate),0.8));
         end 
	 if (max(mu,infeas_meas) < 1e-6) & (iter > 20) & (corr_slow) & (stoplevel)
            msg = 'sqlp stop: lack of progress in corrector'; 
   	    if (printlevel) 
               fprintf('\n  %s',msg);
               fprintf(': mucorr/mu = %3.2f, corr_convg_rate = %3.2f',...
               mucorr/mu,corr_convg_rate); 
            end
            runhist.cputime(iter+1) = cputime-tstart; 
            termcode = -1; 
            breakyes = 1;
         else
            update_iter = 1;
         end
      end 
%%
%%---------------------------------------------------------------
%% udpate iterate
%%---------------------------------------------------------------
%%
      indef = [1,1]; 
      if (update_iter)
         for t = 1:5
            [Xchol,indef(1)] = blkcholfun(blk,ops(X,'+',dX,pstep)); 
            timenew = cputime;
            ttime.pchol = ttime.pchol + timenew-timeold; timeold = timenew;
            if (indef(1)); pstep = 0.8*pstep; else; break; end            
         end
	 if (t > 1); pstep = gamused*pstep; end
	 for t = 1:5
            [Zchol,indef(2)] = blkcholfun(blk,ops(Z,'+',dZ,dstep)); 
            timenew = cputime;
            ttime.dchol = ttime.dchol + timenew-timeold; timeold = timenew; 
            if (indef(2)); dstep = 0.8*dstep; else; break; end             
         end
	 if (t > 1); dstep = gamused*dstep; end
         AXtmp = AX + pstep*AXfun(blk,At,par.permA,dX);
         prim_infeasnew = norm(b-AXtmp)/normb2;
         if (rel_gap < 5*infeas_meas); alpha = 1e2; else; alpha = 1e3; end
         if any(indef)
            msg = 'sqlp stop: X, Z not both positive definite'; 
            if (printlevel); fprintf('\n  %s',msg); end
            termcode = -3;
            breakyes = 1;         
         elseif (prim_infeasnew > max([1e-8,rel_gap,20*prim_infeas])) ...
            | (prim_infeasnew > 1e3*max([1e-10,prim_infeas]) & rel_gap < 1e-2) ...
            | (prim_infeasnew > alpha*max([1e-9,param.prim_infeas_min]) ...
               & (prim_infeasnew > 3*prim_infeas) ...
               & (iter > 25) & (dual_infeas < 1e-6) & (rel_gap < 0.1))
            if (stoplevel) & (iter > 1)
               msg = 'sqlp stop: primal infeas has deteriorated too much'; 
               if (printlevel); fprintf('\n  %s, %2.1e',msg,prim_infeasnew); end
               termcode = -7; 
               breakyes = 1; 
            end
         else
            X = ops(X,'+',dX,pstep);  
            y = y + dstep*dy;           
            Z = ops(Z,'+',dZ,dstep);
         end
      end
%%---------------------------------------------------------------
%% adjust linear blk arising from unrestricted blk
%%---------------------------------------------------------------
%%
      for p = 1:size(blk,1)
         if (ublkidx(p) == 1)
            len = blk{p,2}/2;
            alpha = 0.8; 
            xtmp = min(X{p}([1:len]),X{p}(len+[1:len])); 
            X{p}([1:len]) = X{p}([1:len]) - alpha*xtmp;
            X{p}(len+[1:len]) = X{p}(len+[1:len]) - alpha*xtmp;
            if (mu < 1e-8)
               Z{p} = 0.5*mu./max(1,X{p});
	    else
               ztmp = min(1,max(Z{p}([1:len]),Z{p}(len+[1:len]))); 
               beta1 = xtmp'*(Z{p}([1:len])+Z{p}(len+[1:len]));
               beta2 = (X{p}([1:len])+X{p}(len+[1:len])-2*xtmp)'*ztmp;
               beta = max(0.1,min(beta1/beta2,0.5));
               Z{p}([1:len]) = Z{p}([1:len]) + beta*ztmp;
               Z{p}(len+[1:len]) = Z{p}(len+[1:len]) + beta*ztmp;
            end
         end
      end
%%
%%---------------------------------------------------------------
%% compute rp, Rd, infeasibities, etc.
%%---------------------------------------------------------------
%%
      AX  = AXfun(blk,At,par.permA,X); 
      rp  = b-AX;
      ZpATy = ops(Z,'+',Atyfun(blk,At,par.permA,par.isspAy,y));
      ZpATynorm = ops(ZpATy,'norm');
      Rd   = ops(C,'-',ZpATy);
      objadd = blkbarrier(blk,X,Z,Xchol,Zchol,parbarrier) + objadd0; 
      obj = (normb*normC)*[blktrace(blk,C,X), b'*y] + objadd;  
      gap = (normb*normC)*blktrace(blk,X,Z) - diff(objadd);
      rel_gap = gap/(1+sum(abs(obj))); 
      prim_infeas = norm(rp)/normb2;
      dual_infeas = ops(Rd,'norm')/normC2;
      infeas_meas = max(prim_infeas,dual_infeas); 
      if (scale_data)
         infeas_org(1) = prim_infeas*normb;
         infeas_org(2) = dual_infeas*normC;
      end
      homRd = inf; homrp = inf; 
      if (ops(parbarrier,'norm') == 0)
         if (obj(2) > 0); homRd = ZpATynorm/(obj(2)); end
         if (obj(1) < 0); homrp = norm(AX)/(-obj(1))/(normC); end
      end
      trXZ = blktrace(blk,X,Z,parbarrier); 
      if (nn > 0); mu = trXZ/nn; else; mu = gap/ops(X,'getM'); end
%%
      runhist.pobj(iter+1)  = obj(1); 
      runhist.dobj(iter+1)  = obj(2); 
      runhist.gap(iter+1)   = gap;
      runhist.relgap(iter+1)  = rel_gap;
      runhist.pinfeas(iter+1) = prim_infeas;
      runhist.dinfeas(iter+1) = dual_infeas;
      runhist.infeas(iter+1)  = infeas_meas;
      runhist.step(iter+1)    = min(pstep,dstep); 
      runhist.cputime(iter+1) = cputime-tstart; 
      timenew = cputime;
      ttime.misc = ttime.misc + timenew-timeold; timeold = timenew;  
      [hh,mm,ss] = mytime(sum(runhist.cputime)); 
      if (printlevel>=3)
         fprintf('\n%2.0f  %4.3f %4.3f',iter,pstep,dstep);
         fprintf(' %2.1e %2.1e  %2.1e',prim_infeas,dual_infeas,gap);
         fprintf(' %- 7.6e  %s:%s:%s',mean(obj),hh,mm,ss);
      end
%%
%%--------------------------------------------------
%% check convergence.
%%--------------------------------------------------
%%
      param.iter        = iter; 
      param.obj         = obj;
      param.rel_gap     = rel_gap; 
      param.gap         = gap; 
      param.mu          = mu; 
      param.prim_infeas = prim_infeas;
      param.dual_infeas = dual_infeas;
      param.homRd       = homRd; 
      param.homrp       = homrp; 
      param.AX          = AX; 
      param.ZpATynorm   = ZpATynorm;
      param.normX       = ops(X,'norm'); 
      param.normZ       = ops(Z,'norm'); 
      param.stoplevel   = stoplevel; 
      param.termcode    = termcode; 
%%
      if (~breakyes)
         [param,breakyes,restart,msg2] = sqlpcheckconvg(param,runhist); 
      end
      if (breakyes); break; end
      if (restart)
         [X,y,Z] = infeaspt(blk,At,C,b,2,1e5); 
         rp  = b-AXfun(blk,At,par.permA,X); 
         ZpATy = ops(Z,'+',Atyfun(blk,At,par.permA,par.isspAy,y));
         Rd  = ops(C,'-',ZpATy); 
         trXZ = blktrace(blk,X,Z,parbarrier); 
         mu   = trXZ/nn;
         gap  =  (normb*normC)*blktrace(blk,X,Z) - diff(objadd);
         prim_infeas = norm(rp)/normb2;
         dual_infeas = ops(Rd,'norm')/normC2;
         infeas_meas = max(prim_infeas,dual_infeas); 
         [Xchol,indef(1)] = blkcholfun(blk,X); 
         [Zchol,indef(2)] = blkcholfun(blk,Z); 
         stoplevel = 3;
      end
   end
%%
%%---------------------------------------------------------------
%% end of main loop
%%---------------------------------------------------------------
%%
%%---------------------------------------------------------------
%% unscale and produce infeasibility certificates if appropriate
%%---------------------------------------------------------------
%%
   if (iter >= 1)
      [X,y,Z,termcode,resid,reldist,msg3] = ...
      sqlpmisc(blk,At,C,b,X,y,Z,par.permZ,param); 
   end
%%
%%---------------------------------------------------------------
%% recover unrestricted blk from linear blk
%%---------------------------------------------------------------
%% 
   for p = 1:size(blk,1)
      if (ublkidx(p) == 1)
         n = blk{p,2}/2; 
         X{p} = X{p}(1:n)-X{p}(n+[1:n]); 
         Z{p} = Z{p}(1:n); 
      end
   end
   if (analytic_prob)
      X = X(1:numblkold); Z = Z(1:numblkold); 
   end
%%
%%---------------------------------------------------------------
%% print summary
%%---------------------------------------------------------------
%%
   maxC = 1+ops(ops(C,'abs'),'max'); 
   maxb = 1+max(abs(b)); 
   if (scale_data)
      dimacs = [infeas_org(1)*normb2/maxb; 0; infeas_org(2)*normC2/maxC; 0]; 
   else
      dimacs = [prim_infeas*normb2/maxb; 0; dual_infeas*normC2/maxC; 0];
   end
   dimacs = [dimacs; [-diff(obj); gap]/(1+sum(abs(obj)))];
   info.dimacs   = dimacs; 
   info.termcode = termcode;
   info.iter     = iter; 
   info.obj      = obj; 
   info.gap      = gap; 
   info.relgap   = rel_gap;
   info.pinfeas  = prim_infeas;
   info.dinfeas  = dual_infeas;
   info.cputime  = sum(runhist.cputime); 
   info.ttime    = ttime; 
   info.resid    = resid;
   info.reldist  = reldist; 
   info.normX    = ops(X,'norm'); 
   info.normy    = norm(y); 
   info.normZ    = ops(Z,'norm'); 
   info.normb    = normb2; 
   info.normC    = normC2; 
   info.normA    = normA2;
   info.msg1     = msg; 
   info.msg2     = msg2;
   info.msg3     = msg3;
%%
   sqlpsummary(info,ttime,infeas_org,printlevel);
%%*****************************************************************************
