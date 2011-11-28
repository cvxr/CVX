%%*****************************************************************************
%% sdpt3: solve an semidefinite-quadratic-linear program 
%%        by infeasible path-following method on the original problem 
%%        or the 3-parameter homogeneous self-dual model. 
%%
%%  [obj,X,y,Z,info,runhist] = sdpt3(blk,At,C,b,OPTIONS,X0,y0,Z0);
%%
%%  Input: blk: a cell array describing the block diagonal structure of SQL data.
%%          At: a cell array with At{p} = [svec(Ap1) ... svec(Apm)] 
%%         b,C: data for the SQL instance.
%%  (X0,y0,Z0): an initial iterate (if it is not given, the default is used).
%%     OPTIONS: a structure that specifies parameters required in sdpt3.m,
%%              (if it is not given, the default in sqlparameters.m is used). 
%%
%%  Output: obj = [<C,X>, <b,y>].
%%          (X,y,Z): an approximately optimal solution or a primal or dual
%%                   infeasibility certificate. 
%%          info.termcode = termination-code  
%%          info.iter     = number of iterations
%%          info.obj      = [primal-obj, dual-obj]
%%          info.cputime  = total-time
%%          info.gap      = gap
%%          info.pinfeas  = primal_infeas
%%          info.dinfeas  = dual_infeas  
%%          runhist.pobj    = history of primal objective value. 
%%          runhist.dobj    = history of dual   objective value.
%%          runhist.gap     = history of <X,Z>. 
%%          runhist.pinfeas = history of primal infeasibility. 
%%          runhist.dinfeas = history of dual   infeasibility. 
%%          runhist.cputime = history of cputime spent.
%%----------------------------------------------------------------------------
%%  The OPTIONS structure specifies the required parameters: 
%%      vers  gam  predcorr  expon  gaptol  inftol  steptol  
%%      maxit  printlevel  scale_data ...
%%      (all have default values set in sqlparameters.m).
%%*************************************************************************
%% SDPT3: version 4.0
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 20 Dec 2007
%%*************************************************************************

  function [obj,X,y,Z,info,runhist] = sdpt3(blk,At,C,b,OPTIONS,X0,y0,Z0);

%%                                      
%%-----------------------------------------
%% get parameters from the OPTIONS structure 
%%-----------------------------------------
%%
   warning off; 

   matlabversion = sscanf(version,'%f');
   if strcmp(computer,'PCWIN64') | strcmp(computer,'GLNXA64')
      par.computer = 64; 
   else
      par.computer = 32; 
   end
   model = 0; %% automatically decide between sqlp.m and HSDsqlp.m
   par.matlabversion = matlabversion(1); 
   par.vers        = 0; 
   par.predcorr    = 1; 
   par.gam         = 0; 
   par.expon       = 1; 
   par.gaptol      = 1e-8;
   par.inftol      = 1e-8;
   par.steptol     = 1e-6;
   par.maxit       = 100;
   par.printlevel  = 3;
   par.stoplevel   = 1; 
   par.scale_data  = 0;
   par.spdensity   = 0.4; 
   par.rmdepconstr = 0; 
   par.smallblkdim = 40; 
   par.schurfun     = cell(size(blk,1),1);
   par.schurfun_par = cell(size(blk,1),1); 
%%
   parbarrier  = cell(size(blk,1),1); 
   for p = 1:size(blk,1)
      pblk = blk(p,:); 
      if strcmp(pblk{1},'s') | strcmp(pblk{1},'q')
         parbarrier{p} = zeros(1,length(pblk{2}));
      elseif strcmp(pblk{1},'l') | strcmp(pblk{1},'u' )
         parbarrier{p} = zeros(1,sum(pblk{2}));
      end
   end
   parbarrier_0 = parbarrier;   
%%
   if exist('OPTIONS')
      if isfield(OPTIONS,'model')  
         model = OPTIONS.model; 
         if all(model-[0,1,2]); error(' model must be 0, 1 or 2'); end
      end
      if isfield(OPTIONS,'vers');        par.vers     = OPTIONS.vers; end
      if isfield(OPTIONS,'predcorr');    par.predcorr = OPTIONS.predcorr; end 
      if isfield(OPTIONS,'gam');         par.gam      = OPTIONS.gam; end
      if isfield(OPTIONS,'expon');       par.expon    = OPTIONS.expon; end
      if isfield(OPTIONS,'gaptol');      par.gaptol   = OPTIONS.gaptol; end
      if isfield(OPTIONS,'inftol');      par.inftol   = OPTIONS.inftol; end
      if isfield(OPTIONS,'steptol');     par.steptol  = OPTIONS.steptol; end
      if isfield(OPTIONS,'maxit');       par.maxit    = OPTIONS.maxit; end
      if isfield(OPTIONS,'printlevel');  par.printlevel  = OPTIONS.printlevel; end 
      if isfield(OPTIONS,'stoplevel');   par.stoplevel   = OPTIONS.stoplevel; end 
      if isfield(OPTIONS,'scale_data');  par.scale_data  = OPTIONS.scale_data; end
      if isfield(OPTIONS,'spdensity');   par.spdensity   = OPTIONS.spdensity; end
      if isfield(OPTIONS,'rmdepconstr'); par.rmdepconstr = OPTIONS.rmdepconstr; end
      if isfield(OPTIONS,'smallblkdim'); par.smallblkdim = OPTIONS.smallblkdim; end
      if isfield(OPTIONS,'parbarrier');     
         parbarrier = OPTIONS.parbarrier;
         if isempty(parbarrier); parbarrier = parbarrier_0; end
         if ~iscell(parbarrier); 
            tmp = parbarrier; clear parbarrier; parbarrier{1} = tmp; 
         end
         if (length(parbarrier) < size(blk,1))
            len = length(parbarrier); 
            parbarrier(len+1:size(blk,1)) = parbarrier_0(len+1:size(blk,1)); 
         end
      end      
      if isfield(OPTIONS,'schurfun');    
         par.schurfun = OPTIONS.schurfun; 
         if ~isempty(par.schurfun); par.scale_data = 0; end
      end
      if isfield(OPTIONS,'schurfun_par'); par.schurfun_par = OPTIONS.schurfun_par; end
      if isempty(par.schurfun);     par.schurfun     = cell(size(blk,1),1); end
      if isempty(par.schurfun_par); par.schurfun_par = cell(size(blk,1),1); end
   end
   if (size(blk,2) > 2); par.smallblkdim = 0; end
%%
%%-----------------------------------------
%% Add a redundant constraint if there is
%% no linear constraint
%%-----------------------------------------
%%
    isemptyAtb = 0; 
    if isempty(At) & isempty(b);
       %% Add redundant constraint: <-I,X> <= 0
       b = 0; 
       At = ops(ops(blk,'identity'),'*',-1);  
       numblk = size(blk,1); 
       blk{numblk+1,1} = 'l'; blk{numblk+1,2} = 1;  
       At{numblk+1,1} = 1; C{numblk+1,1} = 0; 
       parbarrier{numblk+1} = 0; 
       isemptyAtb = 1; 
    end
%%
%%-----------------------------------------
%% convert matrices to cell arrays 
%%-----------------------------------------
%%
   if ~iscell(At); At = {At}; end;
   if ~iscell(C);  C = {C}; end;
   if all(size(At) == [size(blk,1), length(b)]); 
      convertyes = zeros(size(blk,1),1); 
      for p = 1:size(blk,1)
         if strcmp(blk{p,1},'s') & all(size(At{p,1}) == sum(blk{p,2}))
            convertyes(p) = 1;    
         end
      end
      if any(convertyes)
         if (par.printlevel); 
            fprintf('\n sdpt3: converting At into required format'); 
         end
         At = svec(blk,At,ones(size(blk,1),1));
      end
   end
%%
%%-----------------------------------------
%% validate SQLP data 
%%-----------------------------------------
%%
   tstart = cputime; 
   [blk,At,C,b,blkdim,numblk,parbarrier] = validate(blk,At,C,b,par,parbarrier);
   [blk,At,C,b,iscmp] = convertcmpsdp(blk,At,C,b);
   if (iscmp) & (par.printlevel>=2); 
      fprintf('\n SQLP has complex data'); 
   end
   exist_analytic_term = 0; 
   for p = 1:size(blk,1);
      idx = find(parbarrier{p} > 0); 
      if ~isempty(idx); exist_analytic_term = 1; end
   end
   if (par.printlevel>=2)
      fprintf('\n num. of constraints = %2.0d',length(b));      
      if blkdim(1); 
         fprintf('\n dim. of sdp    var  = %2.0d,',blkdim(1)); 
         fprintf('   num. of sdp  blk  = %2.0d',numblk(1)); 
      end
      if blkdim(2); 
         fprintf('\n dim. of socp   var  = %2.0d,',blkdim(2)); 
         fprintf('   num. of socp blk  = %2.0d',numblk(2)); 
      end
      if blkdim(3); fprintf('\n dim. of linear var  = %2.0d',blkdim(3)); end
      if blkdim(4); fprintf('\n dim. of free   var  = %2.0d',blkdim(4)); end
   end
%%
%%-----------------------------------------
%% initial iterate
%%-----------------------------------------
%%
   if (nargin <= 5) | (isempty(X0) | isempty(y0) | isempty(Z0)); 
      par.startpoint = 1; 
      [X0,y0,Z0] = infeaspt(blk,At,C,b); 
   else
      par.startpoint = 2; 
      if ~iscell(X0);  X0 = {X0}; end;
      if ~iscell(Z0);  Z0 = {Z0}; end;
      y0 = real(y0); 
      if (length(y0) ~= length(b));
         error('sdpt3: length of b and y0 not compatible'); 
      end   
      [X0,Z0] = validate_startpoint(blk,X0,Z0,par.spdensity,iscmp); 
   end
%%
%%-----------------------------------------
%% detect unrestricted blocks in linear blocks
%%-----------------------------------------
%%   
   user_supplied_schurfun = 0; 
   for p = 1:size(blk,1)
      if ~isempty(par.schurfun{p}); user_supplied_schurfun = 1; end
   end
   if (user_supplied_schurfun == 0)
      [blk2,At2,C2,ublkinfo,parbarrier2,X02,Z02] = ...
      detect_ublk(blk,At,C,parbarrier,X0,Z0,par.printlevel);
   else
      blk2 = blk; At2 = At; C2 = C; 
      parbarrier2 = parbarrier; X02 = X0; Z02 = Z0; 
      ublkinfo = cell(size(blk2,1),1); 
   end
   ublksize = blkdim(4); 
   for p = 1:size(ublkinfo,1)
      ublksize = ublksize + length(ublkinfo{p}); 
   end
%%
%%-----------------------------------------
%% detect diagonal blocks in semidefinite blocks
%%-----------------------------------------
%%
   if (user_supplied_schurfun==0)
      [blk3,At3,C3,diagblkinfo,diagblkchange,parbarrier3,X03,Z03] = ...
      detect_lblk(blk2,At2,C2,b,parbarrier2,X02,Z02,par.printlevel); 
   else
      blk3 = blk2; At3 = At2; C3 = C2; 
      parbarrier3 = parbarrier2; X03 = X02; Z03 = Z02; 
      diagblkchange = 0;  
      diagblkinfo = cell(size(blk3,1),1); 
   end
%%
%%-----------------------------------------
%% main solver
%% model = 1: use sqlp
%%       = 2: use HSDsqlp
%%-----------------------------------------
%%
   if (par.vers == 0); 
      if blkdim(1); par.vers = 1; else; par.vers = 2; end
   end
   par.blkdim   = blkdim;
   par.ublksize = ublksize;
   if (exist_analytic_term); 
      model = 1; 
   else
      if (model == 0)
         if (ublksize > 0); model = 2; else; model = 1; end
      end
   end
   if (model==1)
      [obj,X3,y,Z3,info,runhist] = ...
      sqlpmain(blk3,At3,C3,b,par,parbarrier3,X03,y0,Z03);
   elseif (model==2)
      if (nargin <= 5) | (isempty(X0) | isempty(y0) | isempty(Z0)); 
         if (max([ops(At3,'norm'),ops(C3,'norm'),norm(b)]) > 1e2)
            [X03,y03,Z03] = infeaspt(blk3,At3,C3,b,1);
         else
            [X03,y03,Z03] = infeaspt(blk3,At3,C3,b,2,1);
         end
      end
      [obj,X3,y,Z3,info,runhist] = ...
      HSDsqlpmain(blk3,At3,C3,b,par,X03,y0,Z03);  
   end
   info.model = model; 
%%
%%-----------------------------------------
%% recover semidefinite blocks from linear blocks
%%-----------------------------------------
%%
   if any(diagblkchange)
      X2 = cell(size(blk2,1),1); Z2 = cell(size(blk2,1),1); 
      count = 0; 
      for p = 1:size(blk2,1)
         pblk = blk2(p,:); 
         n = sum(pblk{2}); 
         blkno      = diagblkinfo{p,1};
         idxdiag    = diagblkinfo{p,2}; 
         idxnondiag = diagblkinfo{p,3}; 
         if ~isempty(idxdiag)
            len = length(idxdiag); 
            Xtmp = [idxdiag,idxdiag,X3{end}(count+[1:len]); n, n, 0];
            Ztmp = [idxdiag,idxdiag,Z3{end}(count+[1:len]); n, n, 0];
            if ~isempty(idxnondiag)
               [ii,jj,vv] = find(X3{blkno});
               Xtmp = [Xtmp; idxnondiag(ii),idxnondiag(jj),vv];
               [ii,jj,vv] = find(Z3{blkno}); 
               Ztmp = [Ztmp; idxnondiag(ii),idxnondiag(jj),vv];
            end
            X2{p} = spconvert(Xtmp); 
            Z2{p} = spconvert(Ztmp); 
	    count = count + len; 
         else
	    X2(p) = X3(blkno); Z2(p) = Z3(blkno);    
         end
      end
   else
      X2 = X3; Z2 = Z3;
   end
%%
%%-----------------------------------------
%% recover linear block from unrestricted block
%%-----------------------------------------
%%
   numblk = size(blk,1); 
   numblknew = numblk; 
   X = cell(numblk,1); Z = cell(numblk,1); 
   for p = 1:numblk
      n = blk{p,2}; 
      if isempty(ublkinfo{p,1})
         X{p} = X2{p}; Z{p} = Z2{p}; 
      else
	 Xtmp = zeros(n,1); Ztmp = zeros(n,1); 
         Xtmp(ublkinfo{p,1}) = max(0,X2{p}); 
         Xtmp(ublkinfo{p,2}) = max(0,-X2{p});
         Ztmp(ublkinfo{p,1}) = max(0,Z2{p}); 
         Ztmp(ublkinfo{p,2}) = max(0,-Z2{p});
         if ~isempty(ublkinfo{p,3})
            numblknew = numblknew + 1; 
            Xtmp(ublkinfo{p,3}) = X2{numblknew}; 
            Ztmp(ublkinfo{p,3}) = Z2{numblknew}; 
         end
	 X{p} = Xtmp; Z{p} = Ztmp;
      end     
   end
%%
%%-----------------------------------------
%% recover complex solution
%%-----------------------------------------
%%
   if (iscmp)
      for p = 1:numblk
         pblk = blk(p,:); 
         n = sum(pblk{2})/2;
         if strcmp(pblk{1},'s'); 
            X{p} = X{p}(1:n,1:n) + sqrt(-1)*X{p}(n+[1:n],1:n); 
            Z{p} = Z{p}(1:n,1:n) + sqrt(-1)*Z{p}(n+[1:n],1:n); 
            X{p} = 0.5*(X{p}+X{p}');  
            Z{p} = 0.5*(Z{p}+Z{p}');  
         end
      end
   end
   if (isemptyAtb)
      X = X(1:end-1); Z = Z(1:end-1); 
   end
%%*****************************************************************************
