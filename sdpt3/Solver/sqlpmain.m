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

   global spdensity smallblkdim  printlevel msg
   global solve_ok  use_LU  exist_analytic_term numpertdiagschur  
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
   smallblkdim   = par.smallblkdim;
   schurfun      = par.schurfun;
   schurfun_par  = par.schurfun_par;
   ublksize      = par.ublksize; 
%%
   tstart = clock; 
   X = X0; y = y0; Z = Z0; 
   for p = 1:size(blk,1)
      if strcmp(blk{p,1},'u'); Z{p} = zeros(blk{p,2},1); end
   end
%%
%%-----------------------------------------
%% convert unrestricted blk to linear blk. 
%%-----------------------------------------
%%
   convertlen = 0; 
   [blk,At,C,X,Z,u2lblk,ublkidx] = sqlpu2lblk(blk,At,C,X,Z,par,convertlen);
   for p = 1:size(blk,1) 
      pblk = blk(p,:); 
      if (u2lblk(p) == 1) 
         n = 2*blk{p,2}; 
         blk{p,1} = 'l';  blk{p,2} = n;
         parbarrier{p} = zeros(1,n);
         At{p} = [At{p}; -At{p}];  
         tau = max(1,norm(C{p})); 
         C{p} = [C{p}; -C{p}]; 
         msg = 'convert ublk to lblk'; 
         if (printlevel); fprintf(' *** %s',msg); end
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
            X{p} = X{p}.*(1+1e-10*randmat(n,1,0,'u')); 
            Z{p} = Z{p}.*(1+1e-10*randmat(n,1,0,'u')); 
	 else
            const = max(abs(X{p})) + 100; 
            X{p} = [X{p}+const; const*ones(n/2,1)]; 
            %%old: const = 100; Z{p} = [const*ones(n/2,1); const*ones(n/2,1)];
            Z{p} = [abs(Z0{p}); abs(Z0{p})] + 1e-4; 
         end
      end
   end
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
   normA2 = 1+ops(At,'norm'); 
   normb2 = 1+norm(b); 
   normC2 = 1+ops(C,'norm'); 
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
         elseif strcmp(pblk{1},'q') 
            nn = nn + sum(pblk{2}(idx)); 
         elseif strcmp(pblk{1},'s')   
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
      u2lblk(numblk,1) = 0;
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
%%-----------------------------------------
%% initialization
%%-----------------------------------------
%%
   EE = ops(blk,'identity');
   normE2 = ops(EE,'norm'); Zpertold = 1; 
   for p = 1:size(blk,1) 
      normCC(p) = 1+ops(C(p),'norm');
      normEE(p) = 1+ops(EE(p),'norm'); 
   end
   [Xchol,indef(1)] = blkcholfun(blk,X); 
   [Zchol,indef(2)] = blkcholfun(blk,Z); 
   if any(indef)
      msg = 'stop: X or Z not positive definite'; 
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
   relgap = gap/(1+sum(abs(obj)));
   prim_infeas = norm(rp)/normb2;
   dual_infeas = ops(Rd,'norm')/normC2;
   infeas = max(prim_infeas,dual_infeas); 
   if (scale_data)
      infeas_org(1) = prim_infeas*normb;
      infeas_org(2) = dual_infeas*normC;
   else
      infeas_org = [0,0]; 
   end
   trXZ = blktrace(blk,X,Z,parbarrier); 
   if (nn > 0); mu  = trXZ/nn; else; mu = gap/ops(X,'getM'); end
   normX = ops(X,'norm'); 
%%   
   termcode = 0; restart = 0; 
   pstep = 1; dstep = 1; pred_convg_rate = 1; corr_convg_rate = 1;
   prim_infeas_min  = prim_infeas; 
   dual_infeas_min  = dual_infeas; 
   prim_infeas_best = prim_infeas; 
   dual_infeas_best = dual_infeas; 
   infeas_best = infeas; 
   relgap_best = relgap; 
   homRd = inf; homrp = inf; dy = zeros(length(b),1);    
   msg = []; msg2 = []; msg3 = [];
   runhist.pobj    = obj(1);
   runhist.dobj    = obj(2); 
   runhist.gap     = gap;
   runhist.relgap  = relgap;
   runhist.pinfeas = prim_infeas;
   runhist.dinfeas = dual_infeas;
   runhist.infeas  = infeas;  
   runhist.pstep   = 0; 
   runhist.dstep   = 0; 
   runhist.step    = 0; 
   runhist.normX   = normX; 
   runhist.cputime = etime(clock,tstart); 
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
         fprintf('\nit pstep dstep pinfeas dinfeas  gap')
         fprintf('      prim-obj      dual-obj    cputime\n');
         fprintf('------------------------------------------------');
         fprintf('-------------------\n');
         fprintf('%2.0f|%4.3f|%4.3f|%2.1e|%2.1e|',0,0,0,prim_infeas,dual_infeas);
         fprintf('%2.1e|%- 7.6e %- 7.6e| %s:%s:%s|',gap,obj(1),obj(2),hh,mm,ss);
      end
   end
%%
%%---------------------------------------------------------------
%% start main loop
%%---------------------------------------------------------------
%%
   param.termcode    = termcode; 
   param.iter        = 0; 
   param.obj         = obj;
   param.relgap      = relgap; 
   param.prim_infeas = prim_infeas;   param.dual_infeas = dual_infeas;    
   param.homRd       = homRd;         param.homrp       = homrp; 
   param.AX          = AX;            param.ZpATynorm   = ZpATynorm;
   param.normA       = normA;  
   param.normb       = normb;         param.normC       = normC;
   param.normX0      = normX0;        param.normZ0      = normZ0; 
   param.m0          = m0;            param.indeprows   = indeprows;
   param.prim_infeas_bad = 0;         
   param.dual_infeas_bad = 0; 
   param.prim_infeas_min = prim_infeas; 
   param.dual_infeas_min = dual_infeas; 
   param.gaptol      = gaptol;
   param.inftol      = inftol; 
   param.maxit       = maxit;
   param.scale_data  = scale_data;
   param.printlevel  = printlevel; 
   param.ublksize    = ublksize; 
   Xbest = X; ybest = y; Zbest = Z; 
%%
   for iter = 1:maxit;
      tstart  = clock;  
      timeold = tstart;
      update_iter = 0; breakyes = 0; 
      pred_slow = 0; corr_slow = 0; step_short = 0; 
      par.parbarrier = parbarrier; 
      par.iter    = iter; 
      par.obj     = obj; 
      par.relgap  = relgap; 
      par.pinfeas = prim_infeas; 
      par.dinfeas = dual_infeas;
      par.rp      = rp; 
      par.y       = y; 
      par.dy      = dy; 
      par.normX   = normX; 
      par.ZpATynorm = ZpATynorm; 
      %%if (printlevel > 2); fprintf(' %2.1e',par.normX); end
      if (iter == 1 | restart); Cpert = min(1,normC2/ops(EE,'norm')); end
      if (runhist.dinfeas(1) > 1e-3) & (~exist_analytic_term) ...
         & (relgap > 1e-4) 
         if (par.normX > 5e3 & iter < 20)
            Cpert = Cpert*0.5; 
         elseif (par.normX > 5e2 & iter < 20); 
            Cpert = Cpert*0.3; 
         else; 
            Cpert = Cpert*0.1; 
         end
         Rd = ops(Rd,'+',EE,Cpert); 
         %%if (printlevel > 2); fprintf('|%2.1e',Cpert); end
      end
%%---------------------------------------------------------------
%% predictor step.
%%---------------------------------------------------------------
%%
      if (predcorr)
         sigma = 0; 
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
         msg = 'stop: difficulty in computing predictor directions'; 
         if (printlevel); fprintf('\n  %s',msg); end
         runhist.pinfeas(iter+1) = runhist.pinfeas(iter); 
         runhist.dinfeas(iter+1) = runhist.dinfeas(iter); 
         runhist.relgap(iter+1)  = runhist.relgap(iter); 
         runhist.cputime(iter+1) = etime(clock,tstart); 
         termcode = -4;
         break; %% do not ues breakyes = 1
      end
      timenew = clock;
      ttime.pred = ttime.pred + etime(timenew,timeold); timeold = timenew; 
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
      timenew = clock; 
      ttime.pred_pstep = ttime.pred_pstep + etime(timenew,timeold); timeold = timenew;
      Zstep = steplength(blk,Z,dZ,Zchol,invZchol); 
      dstep = min(1,gamused*full(Zstep));
      trXZnew = trXZ + pstep*blktrace(blk,dX,Z,parbarrier) ...
                 + dstep*blktrace(blk,X,dZ,parbarrier) ...
                 + pstep*dstep*blktrace(blk,dX,dZ,parbarrier);
      if (nn > 0); mupred  = trXZnew/nn; else; mupred = 1e-16; end
      mupredhist(iter) = mupred;
      timenew = clock;        
      ttime.pred_dstep = ttime.pred_dstep + etime(timenew,timeold); timeold = timenew;
%%
%%-----------------------------------------
%%  stopping criteria for predictor step.
%%-----------------------------------------
%%
      if (min(pstep,dstep) < steptol) & (stoplevel) & (iter > 10)
         msg = 'stop: steps in predictor too short';
         if (printlevel) 
            fprintf('\n  %s',msg);
            fprintf(': pstep = %3.2e,  dstep = %3.2e\n',pstep,dstep);
         end
         runhist.cputime(iter+1) = etime(clock,tstart); 
         termcode = -2; 
         breakyes = 1; 
      end
      if (~predcorr)
         if (iter >= 2) 
            idx = [max(2,iter-2) : iter];
            pred_slow = all(mupredhist(idx)./mupredhist(idx-1) > 0.4);
            idx = [max(2,iter-5) : iter];
            pred_convg_rate = mean(mupredhist(idx)./mupredhist(idx-1));
            pred_slow = pred_slow + (mupred/mu > 5*pred_convg_rate);
         end 
         if (max(mu,infeas) < 1e-6) & (pred_slow) & (stoplevel)
            msg = 'stop: lack of progress in predictor'; 
            if (printlevel) 
               fprintf('\n  %s',msg);
               fprintf(': mupred/mu = %3.2f, pred_convg_rate = %3.2f.',...
               mupred/mu,pred_convg_rate);
            end
            runhist.cputime(iter+1) = etime(clock,tstart); 
            termcode = -2; 
            breakyes = 1;
         else 
            update_iter = 1; 
         end
      end
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
         if (nn==0)
             sigma = 0.2; 
         elseif (mupred < 0) 
             sigma = 0.8; 
         else
            sigma = min(1, (mupred/mu)^expon_used);
         end
         sigmu = cell(size(blk,1),1); 
         for p = 1:size(blk,1)
            sigmu{p} = max(sigma*mu, parbarrier{p}'); 
         end	 
         if (vers == 1)
            [dX,dy,dZ] = HKMcorr(blk,At,par,rp,Rd,sigmu,hRd,...
             dX,dZ,coeff,L,X,Z);
         elseif (vers == 2)
            [dX,dy,dZ] = NTcorr(blk,At,par,rp,Rd,sigmu,hRd,...
             dX,dZ,coeff,L,X,Z); 
         end
         if (solve_ok <= 0)
            msg = 'stop: difficulty in computing corrector directions'; 
            if (printlevel); fprintf('\n  %s',msg); end
            runhist.pinfeas(iter+1) = runhist.pinfeas(iter); 
            runhist.dinfeas(iter+1) = runhist.dinfeas(iter); 
            runhist.relgap(iter+1)  = runhist.relgap(iter); 
            runhist.cputime(iter+1) = etime(clock,tstart); 
            termcode = -4;
            break; %% do not ues breakyes = 1
         end
         timenew = clock;
         ttime.corr = ttime.corr + etime(timenew,timeold); timeold = timenew; 
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
         timenew = clock;
         ttime.corr_pstep = ttime.corr_pstep+etime(timenew,timeold); timeold = timenew;
         Zstep = steplength(blk,Z,dZ,Zchol,invZchol);
         dstep = min(1,gamused*full(Zstep));
         trXZnew = trXZ + pstep*blktrace(blk,dX,Z,parbarrier) ...
                    + dstep*blktrace(blk,X,dZ,parbarrier)...
                    + pstep*dstep*blktrace(blk,dX,dZ,parbarrier); 
         if (nn > 0); mucorr  = trXZnew/nn; else; mucorr = 1e-16; end
         timenew = clock;
         ttime.corr_dstep = ttime.corr_dstep+etime(timenew,timeold); timeold = timenew;
%%
%%-----------------------------------------
%%  stopping criteria for corrector step
%%-----------------------------------------
         if (iter >= 2) 
            idx = [max(2,iter-2) : iter];
            corr_slow = all(runhist.gap(idx)./runhist.gap(idx-1) > 0.8); 
            idx = [max(2,iter-5) : iter];
            corr_convg_rate = mean(runhist.gap(idx)./runhist.gap(idx-1));
            corr_slow = corr_slow + (mucorr/mu > max(min(1,5*corr_convg_rate),0.8));
         end 
	 if (max(relgap,infeas) < 1e-6) & (iter > 20) ...
            & (corr_slow > 1) & (stoplevel)
            msg = 'stop: lack of progress in corrector'; 
   	    if (printlevel) 
               fprintf('\n  %s',msg);
               fprintf(': mucorr/mu = %3.2f, corr_convg_rate = %3.2f',...
               mucorr/mu,corr_convg_rate); 
            end
            runhist.cputime(iter+1) = etime(clock,tstart); 
            termcode = -2; 
            breakyes = 1;
         else
            update_iter = 1;
         end
      end 
%%---------------------------------------------------------------
%% udpate iterate
%%---------------------------------------------------------------
      indef = [1,1]; 
      if (update_iter)
         for t = 1:5
            [Xchol,indef(1)] = blkcholfun(blk,ops(X,'+',dX,pstep)); 
            timenew = clock;
            ttime.pchol = ttime.pchol + etime(timenew,timeold); timeold = timenew;
            if (indef(1)); pstep = 0.8*pstep; else; break; end            
         end
	 if (t > 1); pstep = gamused*pstep; end
	 for t = 1:5
            [Zchol,indef(2)] = blkcholfun(blk,ops(Z,'+',dZ,dstep)); 
            timenew = clock;
            ttime.dchol = ttime.dchol + etime(timenew,timeold); timeold = timenew; 
            if (indef(2)); dstep = 0.8*dstep; else; break; end             
         end
	 if (t > 1); dstep = gamused*dstep; end
         %%-------------------------------------------
         AXtmp = AX + pstep*AXfun(blk,At,par.permA,dX);
         prim_infeasnew = norm(b-AXtmp)/normb2;
         if (relgap < 5*infeas); alpha = 1e2; else; alpha = 1e3; end
         if any(indef)
            if indef(1); msg = 'stop: X not positive definite'; end
            if indef(2); msg = 'stop: Z not positive definite'; end
            if (printlevel); fprintf('\n  %s',msg); end
            termcode = -3;
            breakyes = 1;         
         elseif (prim_infeasnew > max([1e-8,relgap,20*prim_infeas]) & iter > 10) ...
            | (prim_infeasnew > max([1e-7,1e3*prim_infeas,0.1*relgap]) & relgap < 1e-2) ...
            | (prim_infeasnew > alpha*max([1e-9,param.prim_infeas_min]) ...
               & (prim_infeasnew > max([3*prim_infeas,0.1*relgap])) ...
               & (iter > 25) & (dual_infeas < 1e-6) & (relgap < 0.1)) ...
            | ((prim_infeasnew > 1e3*prim_infeas & prim_infeasnew > 1e-12) ...
               & (max(relgap,dual_infeas) < 1e-8))
            if (stoplevel) 
               msg = 'stop: primal infeas has deteriorated too much'; 
               if (printlevel); fprintf('\n  %s, %2.1e',msg,prim_infeasnew); end
               termcode = -7; 
               breakyes = 1; 
            end
         elseif (trXZnew > 1.05*runhist.gap(iter)) & (~exist_analytic_term) ...
	    & ((infeas < 1e-5) & (relgap < 1e-4) & (iter > 20) ...
	       | (max(infeas,relgap) < 1e-7) & (iter > 10)) 
            if (stoplevel) 
               msg = 'stop: progress in duality gap has deteriorated'; 
               if (printlevel); fprintf('\n  %s, %2.1e',msg,trXZnew); end
               termcode = -8; 
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
      if (~breakyes)
         for p = 1:size(blk,1)
            if (u2lblk(p) == 1)
               len = blk{p,2}/2;              
               xtmp = min(X{p}([1:len]),X{p}(len+[1:len])); 
               alpha = 0.8; 
               X{p}([1:len])     = X{p}([1:len]) - alpha*xtmp;
               X{p}(len+[1:len]) = X{p}(len+[1:len]) - alpha*xtmp;
               if (mu < 1e-4) %% old: (mu < 1e-7)
                  Z{p} = 0.5*mu./max(1,X{p}); %% good to keep this step
               else
                  ztmp = min(1,max(Z{p}([1:len]),Z{p}(len+[1:len])));
                  if (dual_infeas > 1e-4 & dstep < 0.2)
                     beta = 0.3; 
                  else  
                     beta = 0.0; 
                  end
                  %% important to set beta = 0 at later stage. 
                  Z{p}([1:len])     = Z{p}([1:len]) + beta*ztmp;
                  Z{p}(len+[1:len]) = Z{p}(len+[1:len]) + beta*ztmp;
               end
            end
         end
      end
%%--------------------------------------------------
%% perturb Z: do this step before checking for break
%%--------------------------------------------------
      if (~breakyes) & (~exist_analytic_term)
         trXZtmp = blktrace(blk,X,Z);
         trXE  = blktrace(blk,X,EE);
         Zpert = max(1e-12,0.2*min(relgap,prim_infeas)).*normC2./normE2;
         Zpert = min(Zpert,0.1*trXZtmp./trXE);
         Zpert = min([1,Zpert,1.5*Zpertold]); 
         if (infeas < 0.1) 
            Z = ops(Z,'+',EE,Zpert); 
            [Zchol,indef(2)] = blkcholfun(blk,Z);
            if any(indef(2))
               msg = 'stop: Z not positive definite';      
               if (printlevel); fprintf('\n  %s',msg); end
               termcode = -3;
               breakyes = 1; 
            end
            %%if (printlevel > 2); fprintf(' %2.1e',Zpert); end
         end
         Zpertold = Zpert; 
      end
%%---------------------------------------------------------------
%% compute rp, Rd, infeasibities, etc
%%---------------------------------------------------------------
%%
      AX  = AXfun(blk,At,par.permA,X); 
      rp  = b-AX;
      ZpATy = ops(Z,'+',Atyfun(blk,At,par.permA,par.isspAy,y));
      ZpATynorm = ops(ZpATy,'norm');
      Rd  = ops(C,'-',ZpATy);
      objadd = blkbarrier(blk,X,Z,Xchol,Zchol,parbarrier) + objadd0; 
      obj = (normb*normC)*[blktrace(blk,C,X), b'*y] + objadd;  
      gap = (normb*normC)*blktrace(blk,X,Z) - diff(objadd);
      relgap = gap/(1+sum(abs(obj))); 
      prim_infeas = norm(rp)/normb2;
      dual_infeas = ops(Rd,'norm')/normC2;
      infeas = max(prim_infeas,dual_infeas); 
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
      normX = ops(X,'norm');
%%
      runhist.pobj(iter+1)  = obj(1); 
      runhist.dobj(iter+1)  = obj(2); 
      runhist.gap(iter+1)   = gap;
      runhist.relgap(iter+1)  = relgap;
      runhist.pinfeas(iter+1) = prim_infeas;
      runhist.dinfeas(iter+1) = dual_infeas;
      runhist.infeas(iter+1)  = infeas;
      runhist.pstep(iter+1)   = pstep; 
      runhist.dstep(iter+1)   = dstep; 
      runhist.step(iter+1)    = min(pstep,dstep); 
      runhist.normX(iter+1)   = normX; 
      runhist.cputime(iter+1) = etime(clock,tstart); 
      timenew = clock;
      ttime.misc = ttime.misc + etime(timenew,timeold); timeold = timenew;  
      [hh,mm,ss] = mytime(sum(runhist.cputime)); 
      if (printlevel>=3)
         fprintf('\n%2.0f|%4.3f|%4.3f',iter,pstep,dstep);
         fprintf('|%2.1e|%2.1e|%2.1e|',prim_infeas,dual_infeas,gap);
         fprintf('%- 7.6e %- 7.6e| %s:%s:%s|',obj(1),obj(2),hh,mm,ss);
      end
%%--------------------------------------------------
%% check convergence
%%--------------------------------------------------
      param.use_LU      = use_LU; 
      param.stoplevel   = stoplevel; 
      param.termcode    = termcode; 
      param.iter        = iter; 
      param.obj         = obj;
      param.gap         = gap; 
      param.relgap      = relgap; 
      param.prim_infeas = prim_infeas;
      param.dual_infeas = dual_infeas;
      param.mu        = mu; 
      param.homRd     = homRd; 
      param.homrp     = homrp; 
      param.AX        = AX; 
      param.ZpATynorm = ZpATynorm;
      param.normX     = ops(X,'norm'); 
      param.normZ     = ops(Z,'norm'); 
      param.numpertdiagschur = numpertdiagschur; 
      if (~breakyes)
         [param,breakyes,restart,msg2] = sqlpcheckconvg(param,runhist); 
      end
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
         infeas = max(prim_infeas,dual_infeas); 
         [Xchol,indef(1)] = blkcholfun(blk,X); 
         [Zchol,indef(2)] = blkcholfun(blk,Z); 
         stoplevel = 3;
      end
%%--------------------------------------------------
%% check for break
%%--------------------------------------------------
      if ((prim_infeas < 1.5*prim_infeas_best) ...                
         | (max(relgap,infeas) < 0.8*max(relgap_best,infeas_best))) ...
         & (max(relgap,dual_infeas) < 0.8*max(relgap_best,dual_infeas_best)) 
         Xbest = X; ybest = y; Zbest = Z; 
         prim_infeas_best = prim_infeas; 
         dual_infeas_best = dual_infeas; 
         relgap_best = relgap; infeas_best = infeas; 
         update_best(iter+1) = 1; 
         %%fprintf('#')
      else
         update_best(iter+1) = 0; 
      end   
      if (max(relgap_best,infeas_best) < 1e-4 ...
          & norm(update_best(max(1,iter-1):iter+1)) == 0)
         msg = 'lack of progress in infeas'; 
         if (printlevel); fprintf('\n  %s',msg); end
         termcode = -9; 
         breakyes = 1; 
      end
      if (breakyes); break; end
   end
%%---------------------------------------------------------------
%% end of main loop
%%---------------------------------------------------------------
%%
   use_bestiter = 1; 
   if (use_bestiter) & (param.termcode <= 0)
      X = Xbest; y = ybest; Z = Zbest; 
      Xchol = blkcholfun(blk,X); 
      Zchol = blkcholfun(blk,Z);      
      AX = AXfun(blk,At,par.permA,X); 
      rp = b-AX;
      ZpATy = ops(Z,'+',Atyfun(blk,At,par.permA,par.isspAy,y));
      Rd = ops(C,'-',ZpATy);
      objadd = blkbarrier(blk,X,Z,Xchol,Zchol,parbarrier) + objadd0; 
      obj = (normb*normC)*[blktrace(blk,C,X), b'*y] + objadd;  
      gap = (normb*normC)*blktrace(blk,X,Z) - diff(objadd);
      relgap = gap/(1+sum(abs(obj)));
      prim_infeas = norm(rp)/normb2; 
      dual_infeas = ops(Rd,'norm')/normC2; 
      infeas = max(prim_infeas,dual_infeas); 
      runhist.pobj(iter+1)  = obj(1); 
      runhist.dobj(iter+1)  = obj(2); 
      runhist.gap(iter+1)   = gap;
      runhist.relgap(iter+1)  = relgap;
      runhist.pinfeas(iter+1) = prim_infeas;
      runhist.dinfeas(iter+1) = dual_infeas;
      runhist.infeas(iter+1)  = infeas; 
   end
%%---------------------------------------------------------------
%% unscale and produce infeasibility certificates if appropriate
%%---------------------------------------------------------------
   if (iter >= 1)
      [X,y,Z,termcode,resid,reldist,msg3] = ...
      sqlpmisc(blk,At,C,b,X,y,Z,par.permZ,param); 
   end
%%---------------------------------------------------------------
%% recover unrestricted blk from linear blk
%%---------------------------------------------------------------
%% 
   for p = 1:size(blk,1)
      if (u2lblk(p) == 1)
         n = blk{p,2}/2; 
         X{p} = X{p}(1:n)-X{p}(n+[1:n]); 
         Z{p} = Z{p}(1:n); 
      end
   end
   for p = 1:size(ublkidx,1) 
      if ~isempty(ublkidx{p,2})
         n0 = ublkidx{p,1}; idxB = setdiff([1:n0]',ublkidx{p,2});
         tmp = zeros(n0,1); tmp(idxB) = X{p}; X{p} = tmp; 
         tmp = zeros(n0,1); tmp(idxB) = Z{p}; Z{p} = tmp; 
      end
   end
   if (analytic_prob)
      X = X(1:numblkold); Z = Z(1:numblkold); 
   end
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
   info.relgap   = relgap;
   info.pinfeas  = prim_infeas;
   info.dinfeas  = dual_infeas;
   info.cputime  = sum(runhist.cputime); 
   info.time     = ttime; 
   info.resid    = resid;
   info.reldist  = reldist; 
   info.normX    = ops(X,'norm'); 
   info.normy    = norm(y); 
   info.normZ    = ops(Z,'norm'); 
   info.normb    = normb2; info.maxb = maxb; 
   info.normC    = normC2; info.maxC = maxC; 
   info.normA    = normA2;
   info.msg1     = msg; 
   info.msg2     = msg2;
   info.msg3     = msg3;
   sqlpsummary(info,ttime,infeas_org,printlevel);
%%*****************************************************************************
