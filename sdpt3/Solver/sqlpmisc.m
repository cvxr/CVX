%%*****************************************************************************
%% sqlpmisc: 
%% unscale and produce infeasibility certificates if appropriate
%%
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004.
%%*****************************************************************************

   function [X,y,Z,termcode,resid,reldist,msg] = sqlpmisc(blk,At,C,b,X,y,Z,permZ,param);

   termcode    = param.termcode;
   iter        = param.iter; 
   obj         = param.obj;
   relgap      = param.relgap; 
   prim_infeas = param.prim_infeas;
   dual_infeas = param.dual_infeas;
   homRd       = param.homRd; 
   homrp       = param.homrp; 
   AX          = param.AX;
   ZpATynorm   = param.ZpATynorm; 
   m0          = param.m0;    
   indeprows   = param.indeprows;
   normX0      = param.normX0;
   normZ0      = param.normZ0;
   inftol      = param.inftol;
   maxit       = param.maxit; 
   scale_data  = param.scale_data; 
   printlevel  = param.printlevel;
%%
   resid = []; reldist = []; msg = []; 
   if (scale_data)
      normA = param.normA; normC = param.normC; normb = param.normb;
   else
      normA = 1; normC = 1; normb = 1;
   end
   Anorm = ops(At,'norm'); xnorm = ops(X,'norm'); ynorm = norm(y);
   infeas = max(prim_infeas,dual_infeas); 
%%
   if (iter >= maxit) 
      termcode = -6;
      msg = 'sqlp stop: maximum number of iterations reached'; 
      if (printlevel); fprintf('\n  %s',msg); end
   end
   if (termcode <= 0)
      %%
      %% To detect near-infeasibility when the algorithm provides 
      %% a "better" certificate of infeasibility than of optimality.
      %%
      err = max(infeas,relgap);
      iflag = 0;
      if (obj(2) > 0)
         if (homRd < 0.1*sqrt(err*inftol))  
            iflag = 1;
            msg = sprintf('prim_inf,dual_inf,relgap = %3.2e, %3.2e, %3.2e',...
		  prim_infeas,dual_infeas,relgap); 
            if (printlevel); fprintf('\n  %s',msg); end
            termcode = 1;
         end
      end	 
      if (obj(1) < 0)
         if (homrp < 0.1*sqrt(err*inftol)) 
            iflag = 1; 
            msg = sprintf('prim_inf,dual_inf,relgap = %3.2e, %3.2e, %3.2e',...
                  prim_infeas,dual_infeas,relgap); 
            if (printlevel); fprintf('\n  %s',msg); end
            termcode = 2;
         end
      end
      if (iflag == 0)
         if (scale_data == 1)
            X = ops(ops(X,'./',normA),'*',normb);
            y = y*normC;  
            Z = ops(ops(Z,'.*',normA),'*',normC);
         end
      end
   end
   if (termcode == 1) & (iter > 3)
      msg = 'sqlp stop: primal problem is suspected of being infeasible'; 
      if (printlevel); fprintf('\n  %s',msg); end
      if (scale_data == 1)
         X = ops(X,'./',normA); b = b*normb; 
      end
      rby = 1/(b'*y); y = rby*y; Z = ops(Z,'*',rby);
      resid = ZpATynorm * rby;
      reldist = ZpATynorm/(Anorm*ynorm);
   end  
   if (termcode == 2) & (iter > 3)
      msg = 'sqlp stop: dual problem is suspected of being infeasible'; 
      if (printlevel); fprintf('\n  %s',msg); end
      if (scale_data == 1)
         C = ops(C,'.*',normC); 
         Z = ops(Z,'.*',normA);
      end
      tCX = blktrace(blk,C,X);
      X = ops(X,'*',1/(-tCX));
      resid = norm(AX)/(-tCX);
      reldist = norm(AX)/(Anorm*xnorm);
   end
   if (termcode == 3)
      maxblowup = max(ops(X,'norm')/normX0,ops(Z,'norm')/normZ0);
      msg = sprintf('sqlp stop: primal or dual is diverging, %3.1e',maxblowup); 
      if (printlevel); fprintf('\n  %s',msg); end
   end
   [X,Z] = unperm(blk,permZ,X,Z);
   if ~isempty(indeprows)
      ytmp = zeros(m0,1); 
      ytmp(indeprows) = y;
      y = ytmp; 
   end
%%*****************************************************************************
%% unperm: undo the permutations applied in validate.
%%
%% [X,Z] = unperm(blk,permZ,X,Z);
%%
%% undoes the permutation introduced in validate.
%%*****************************************************************************
 
  function [X,Z] = unperm(blk,permZ,X,Z);
%%
  for p = 1:size(blk,1)
     if (strcmp(blk{p,1},'s') & ~isempty(permZ{p}))
        per = permZ{p};
        X{p} = X{p}(per,per);
        Z{p} = Z{p}(per,per);
     end
  end
%%*****************************************************************************
