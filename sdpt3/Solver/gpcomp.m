%%*********************************************************************
%% gpcomp: Compute tp=1/gp in Proposition 2 of the paper: 
%% 
%% R.M. Freund, F. Ordonez, and K.C. Toh,    
%% Behavioral measures and their correlation with IPM iteration counts 
%% on semi-definite programming problems,  
%% Mathematical Programming, 109 (2007), pp. 445--475.
%%
%% [gp,info,Xfeas,blk2,At2,C2,b2] = gpcomp(blk,At,C,b,OPTIONS,solveyes);
%%
%% Xfeas = a feasible X for the primal problem if gp is finite. 
%%         That is, 
%%         norm(b-AXfun(blk,At,[],Xfeas)) 
%%         should be small
%%*********************************************************************

  function [gp,info,Xfeas,blk2,At2,C2,b2] = gpcomp(blk,At,C,b,OPTIONS,solveyes);

  if (nargin < 6); solveyes = 1; end
  if (nargin < 5)
     OPTIONS = sqlparameters; 
     OPTIONS.vers   = 1; 
     OPTIONS.gaptol = 1e-10;
     OPTIONS.printlevel = 3; 
  end
  if isempty(OPTIONS); OPTIONS = sqlparameters; end
  if ~isfield(OPTIONS,'solver'); OPTIONS.solver = 'sqlp'; end
  if ~isfield(OPTIONS,'printlevel'); OPTIONS.printlevel = 3; end
  if ~iscell(C); tmp = C; clear C; C{1} = tmp; end
%%
  m = length(b); 
  blk2 = blk;
  At2 = At; 
  C2 = cell(size(blk,1),1); 
  b2 = zeros(m,1);
%%
%% 
%%  
  dd = ones(1,m); 
  ee = zeros(1,m); EE = cell(size(blk,1),1);
  exist_ublk = 0;
  nn = zeros(size(blk,1),1); 
  for p = 1:size(blk,1)
     pblk = blk(p,:); 
     n = sum(pblk{2}); 
     if strcmp(pblk{1},'s')
        ee = ee + svec(pblk,speye(n),1)'*At{p}; 
        C2{p,1} = sparse(n,n); 
        EE{p} = speye(n); 
        nn(p) = n; 
     elseif strcmp(pblk{1},'q')
        eq = zeros(n,1); 
        idx1 = 1+[0,cumsum(pblk{2})]; 
        idx1 = idx1(1:length(idx1)-1);          
        eq(idx1) = ones(length(idx1),1);
        ee = ee + 2*eq'*At{p}; 
        C2{p,1} = zeros(n,1);
        EE{p} = eq;
        nn(p) = length(pblk{2});  
     elseif strcmp(pblk{1},'l')
        ee = ee + ones(1,n)*At{p}; 
        C2{p,1} = zeros(n,1); 
        EE{p} = ones(n,1); 
        nn(p) = n; 
     elseif strcmp(pblk{1},'u')
        C2{p,1} = zeros(n,1);         
        exist_ublk = 1; 
        EE{p} = sparse(n,1); 
        nn(p) = n; 
     end
     dd = dd + sqrt(sum(At{p}.*At{p}));
  end
  dd = 1./min(1e4,max(1,dd));
  ee = ee.*dd; 
  b  = b.*dd'; 
%%
%% scale data
%%
  D = spdiags(dd',0,m,m); 
  for p = 1:size(blk,1)   
     pblk = blk(p,:); 
     At2{p} = At2{p}*D;
  end
%%
%% New variables in primal problem: 
%% [x; tt; theta].
%%
   numblk = size(blk,1); 
   blk2{numblk+1,1} = 'l'; blk2{numblk+1,2} = 2; 
   At2{numblk+1,1} = [ee; -b']; 
   C2{numblk+1,1} = [-1; 0]; 
%%
%% 3 additional inequality constraints in primal problem. 
%%
   ss = 0; 
   for p = 1:size(blk,1)
      pblk = blk(p,:); 
      n = sum(pblk{2}); 
      if strcmp(pblk{1},'s')
         n2 = sum(pblk{2}.*(pblk{2}+1))/2; 
         At2{p} = [At2{p}, svec(pblk,speye(n,n),1), sparse(n2,2)]; 
         ss = ss + n;
      elseif strcmp(pblk{1},'q')
         eq = zeros(n,1); 
         idx1 = 1+[0,cumsum(pblk{2})]; 
         idx1 = idx1(1:length(idx1)-1);          
         eq(idx1) = ones(length(idx1),1);
         At2{p} = [At2{p}, sparse(eq), sparse(n,2)]; 
         ss = ss + 2*length(pblk{2}); 
      elseif strcmp(pblk{1},'l')
         At2{p} = [At2{p}, sparse(ones(n,1)), sparse(n,2)]; 
         ss = ss + n;
      elseif strcmp(pblk{1},'u')
         At2{p} = [At2{p}, sparse(n,3)]; 
      end
   end
   At2{numblk+1} = sparse([At2{numblk+1}, [ss;0], [0;1], [1;-1]]);
   b2 = [b2; 1; 1; 0];  
%%
%% Add in the linear block corresponding to the 3 slack variables.
%%
   blk2{numblk+2,1} = 'l'; blk2{numblk+2,2} = 3;   
   At2{numblk+2,1} = [sparse(3,m), speye(3,3)]; 
   C2{numblk+2,1}  = zeros(3,1); 
%%
%% Solve SDP
%%
   gp = []; info = []; Xfeas = []; 
   if (solveyes)
      if strcmp(OPTIONS.solver,'sqlp')
         [X0,y0,Z0] = infeaspt(blk2,At2,C2,b2,2,100);    
         [obj,X,y,Z,info] = sqlp(blk2,At2,C2,b2,OPTIONS,X0,y0,Z0); 
      else
         [obj,X,y,Z,info] = HSDsqlp(blk2,At2,C2,b2,OPTIONS); 
      end
      obj = -obj;
      tt = X{numblk+1}(1); theta = X{numblk+1}(2); 
      Xfeas = ops(ops(X(1:numblk),'+',EE(1:numblk),tt),'/',theta); 
%%
      if (obj(1) > 0) | (abs(obj(1)) < 1e-8)
         gp = 1/abs(obj(1));
      elseif (obj(2) > 0)
         gp = 1/obj(2);
      else
         gp = 1/exp(mean(log(abs(obj))));
      end
      err = max(info.dimacs([1,3,6])); 
      if (OPTIONS.printlevel)
         fprintf('\n ******** gp = %3.2e, err = %3.1e\n',gp,err); 
         if (err > 1e-6);
            fprintf('\n----------------------------------------------------')
            fprintf('\n gp problem is not solved to sufficient accuracy');
            fprintf('\n----------------------------------------------------\n')
         end
      end
   end
%%*********************************************************************
