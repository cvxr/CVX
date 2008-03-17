%%*****************************************************************************
%% HSDsqlpmisc: 
%% produce infeasibility certificates if appropriate
%%
%% Input: X,y,Z are the original variables, not the HSD variables. 
%%
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004.
%%*****************************************************************************

   function [X,y,Z,resid,reldist,param,msg] = HSDsqlpmisc(blk,At,C,b,X,y,Z,permZ,param);

   obj         = param.obj;
   relgap      = param.relgap; 
   prim_infeas = param.prim_infeas;
   dual_infeas = param.dual_infeas;
   ZpATynorm   = param.ZpATynorm; 
   inftol      = param.inftol;
   m0          = param.m0;    
   indeprows   = param.indeprows;
   termcode    = param.termcode;
   AX          = param.AX;
   normX0      = param.normX0;
   normZ0      = param.normZ0;
   printlevel  = param.printlevel; 
%%
   resid = []; reldist = []; msg = []; 
   Anorm = ops(At,'norm'); xnorm = ops(X,'norm'); ynorm = norm(y);
%% 
   if (termcode <= 0)
      %%
      %% To detect near-infeasibility when the algorithm provides 
      %% a "better" certificate of infeasibility than of optimality.
      %%
      err = max([prim_infeas,dual_infeas,relgap]);
      iflag = 0;
      if (obj(2) > 0)
         homRd = ZpATynorm/obj(2);
         if (homRd < 1e-2*sqrt(err*inftol))
            iflag = 1;
            termcode = 1;
            param.termcode = 1;
         end
      elseif (obj(1) < 0)
         homrp = norm(AX)/(-obj(1)); 
         if (homrp < 1e-2*sqrt(err*inftol)) 
            iflag = 1; 
            termcode = 2;
            param.termcode = 2;
         end
      end
   end
   if (termcode == 1)
      rby = 1/(b'*y); y = rby*y; Z = ops(Z,'*',rby);
      resid = ZpATynorm * rby;
      reldist = ZpATynorm/(Anorm*ynorm);
      msg = 'Stop: primal problem is suspected of being infeasible'; 
      if (printlevel); fprintf('\n  %s',msg); end
   end  
   if (termcode == 2)
      tCX = blktrace(blk,C,X);
      X = ops(X,'*',1/(-tCX));
      resid = norm(AX) /(-tCX);
      reldist = norm(AX)/(Anorm*xnorm);
      msg = 'Stop: dual problem is suspected of being infeasible'; 
      if (printlevel); fprintf('\n  %s',msg); end
   end
   if (termcode == 3)
      maxblowup = max(ops(X,'norm')/normX0,ops(Z,'norm')/normZ0);
      msg = sprintf('Stop: primal or dual is diverging, %3.1e',maxblowup); 
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
%% [X,Z,Xiter,Ziter] = unperm(blk,permZ,X,Z,Xiter,Ziter);
%%
%% undoes the permutation introduced in validate.
%% can also be called if Xiter and Ziter have not been set as
%%
%% [X,Z] = unperm(blk,permZ,X,Z);
%%
%% SDPT3: version 3.0
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last modified: 2 Feb 01
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
