%%*********************************************************************
%% gdcomp: Compute gd = 1/td in Equation (15) of the paper:
%%
%% R.M. Freund, F. Ordonez, and K.C. Toh,    
%% Behavioral measures and their correlation with IPM iteration counts 
%% on semi-definite programming problems,  
%% Mathematical Programming, 109 (2007), pp. 445--475.
%%
%% [gd,info,yfeas,Zfeas,blk2,At2,C2,b2] = gdcomp(blk,At,C,b,OPTIONS,solveyes);
%%
%% yfeas,Zfeas: a dual feasible pair when gd is finite.
%%              That is, if
%%              Aty = Atyfun(blk,At,[],[],yfeas); 
%%              Rd = ops(C,'-',ops(Zfeas,'+',Aty)); 
%%              then
%%              ops(Rd,'norm') should be small. 
%%
%%*********************************************************************

  function [gd,info,yfeas,Zfeas,blk2,At2,C2,b2] = gdcomp(blk,At,C,b,OPTIONS,solveyes);

  if (nargin < 6); solveyes = 1; end
  if (nargin < 5)
     OPTIONS = sqlparameters; 
     OPTIONS.vers   = 1; 
     OPTIONS.gaptol = 1e-10;
     OPTIONS.printlevel = 3; 
  end
  if isempty(OPTIONS); OPTIONS = sqlparameters; end
  if ~isfield(OPTIONS,'solver'); OPTIONS.solver = 'HSDsqlp'; end
  if ~isfield(OPTIONS,'printlevel'); OPTIONS.printlevel = 3; end
  if ~iscell(C); tmp = C; clear C; C{1} = tmp; end
%%
  m = length(b); 
  blk2 = blk;
  At2 = cell(size(blk,1),1); 
  C2 = cell(size(blk,1),1); 
  EE = cell(size(blk,1),1); 
%%
%% 
%%
  dd = zeros(1,m); 
  alp = 0; 
  beta = 0; 
  for p = 1:size(blk,1)
     pblk = blk(p,:); 
     n = sum(pblk{2}); 
     if strcmp(pblk{1},'s')
        C2{p,1} = sparse(n,n); 
     else
        C2{p,1} = zeros(n,1);
     end
     dd = dd + sqrt(sum(At{p}.*At{p}));
     beta = beta + norm(C{p},'fro'); 
     alp = alp + sqrt(n); 
  end
  alp  = 1./max(1,alp);   
  beta = 1./max(1,beta); 
  dd   = 1./max(1,dd);  
%%
%% New multipliers in dual problem: 
%% [v; tt; theta].
%%
   D = spdiags(dd',0,m,m); 
   
   ss = 0; cc = 0; aa = zeros(1,m); 
   exist_ublk = 0; 
   for p = 1:size(blk,1)
      pblk = blk(p,:); 
      n = sum(pblk{2}); 
      if strcmp(pblk{1},'s')
         At2{p} = [At{p}*D, svec(pblk,alp*speye(n,n),1), -svec(pblk,beta*C{p},1)];
         ss = ss + n; 
         cc = cc + trace(C{p}); 
         aa = aa + svec(pblk,speye(n),1)'*At{p}; 
         EE{p} = speye(n,n); 
      elseif strcmp(pblk{1},'q')
         eq = zeros(n,1); 
         idx1 = 1+[0,cumsum(pblk{2})]; 
         idx1 = idx1(1:length(idx1)-1);          
         eq(idx1) = ones(length(idx1),1);
         At2{p} = [At{p}*D, 2*sparse(alp*eq), -sparse(beta*C{p})];          
         ss = ss + 2*length(pblk{2}); 
         cc = cc + sum(C{p}(idx1)); 
         aa = aa + eq'*At{p}; 
         EE{p} = eq; 
      elseif strcmp(pblk{1},'l')
         el = ones(n,1); 
         At2{p} = [At{p}*D, sparse(alp*el), -sparse(beta*C{p})]; 
         ss = ss + n;
         cc = cc + el'*C{p}; 
         aa = aa + el'*At{p}; 
         EE{p} = el; 
      elseif strcmp(pblk{1},'u')
         At2{p} = [At{p}*D, sparse(n,1), -sparse(beta*C{p})]; 
         exist_ublk = 1; 
         EE{p} = sparse(n,1); 
      end
   end
   aa = aa.*dd;
   cc = cc*beta; 
%%
%% 4 additional inequality constraints in dual problem.
%%
   numblk = size(blk,1); 
   blk2{numblk+1,1} = 'l'; blk2{numblk+1,2} = 4; 
   C2{numblk+1,1}  = [1; 1; 0; 0]; 
   At2{numblk+1,1} = [-aa,         0,   cc; 
		     zeros(1,m),   0,   beta;
		     zeros(1,m),  alp, -beta
                     zeros(1,m), -alp,  0];
   At2{numblk+1} = sparse(At2{numblk+1}); 
   b2 = [zeros(m,1); alp; 0];
%%
%% Solve SDP
%%
   gd = []; info = []; yfeas = []; Zfeas = [];
   if (solveyes)
      if strcmp(OPTIONS.solver,'sqlp')
         [X0,y0,Z0] = infeaspt(blk2,At2,C2,b2,2,100);    
         [obj,X,y,Z,info] = sqlp(blk2,At2,C2,b2,OPTIONS,X0,y0,Z0); 
      else
         [obj,X,y,Z,info] = HSDsqlp(blk2,At2,C2,b2,OPTIONS); 
      end
      tt = alp*abs(y(m+1)); theta = beta*abs(y(m+2)); 
      yfeas = D*y(1:m)/theta; 
      Zfeas = ops(ops(Z(1:numblk),'+',EE,tt),'/',theta); 
      %%
      if (obj(2) > 0) | (abs(obj(2)) < 1e-8)
         gd = 1/abs(obj(2));
      elseif (obj(1) > 0)
         gd = 1/obj(1);
      else
         gd = 1/exp(mean(log(abs(obj))));
      end
      err = max(info.dimacs([1,3,6])); 
      if (OPTIONS.printlevel)
         fprintf('\n ******** gd = %3.2e, err = %3.1e\n',gd,err); 
         if (err > 1e-6);
            fprintf('\n----------------------------------------------------')
            fprintf('\n gd problem is not solved to sufficient accuracy');
            fprintf('\n----------------------------------------------------\n')
         end
      end
   end
%%*********************************************************************


