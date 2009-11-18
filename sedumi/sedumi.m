function [x,y,info] = sedumi(A,b,c,K,pars)
%      [x,y,info] = sedumi(A,b,c,K,pars)
%
% SEDUMI  Self-Dual-Minimization/ Optimization over self-dual homogeneous
%         cones.
%
% >  X = SEDUMI(A,b,c) yields an optimal solution to the linear program
%      MINIMIZE c'*x SUCH THAT A*x = b, x >= 0
%      x is a vector of decision variables.
%      If size(A,2)==length(b), then it solves the linear program
%      MINIMIZE c'*x SUCH THAT A'*x = b, x >= 0
%
% >  [X,Y,INFO] = SEDUMI(A,b,c) also yields a vector of dual multipliers Y,
%      and a structure INFO, with the fields INFO.pinf, INFO.dinf and
%      INFO.numerr.
%
%    (1) INFO.pinf=INFO.dinf=0: x is an optimal solution (as above)
%      and y certifies optimality, viz.\ b'*y = c'*x and c - A'*y >= 0.
%      Stated otherwise, y is an optimal solution to
%      MAXIMIZE b'*y SUCH THAT c-A'*y >= 0.
%      If size(A,2)==length(b), then y solves the linear program
%      MAXIMIZE b'*y SUCH THAT c-A*y >= 0.
%
%    (2) INFO.pinf=1: there cannot be x>=0 with A*x=b, and this is certified
%      by y, viz. b'*y > 0 and A'*y <= 0. Thus y is a Farkas solution.
%
%    (3) INFO.dinf=1: there cannot be y such that c-A'*y >= 0, and this is
%      certified by x, viz. c'*x <0, A*x = 0, x >= 0. Thus x is a Farkas
%      solution.
%
%    (I)   INFO.numerr = 0: desired accuracy achieved (see PARS.eps).
%    (II)  INFO.numerr = 1: numerical problems warning. Results are accurate
%          merely to the level of PARS.bigeps.
%    (III) INFO.numerr = 2: complete failure due to numerical problems.
%
%    INFO.feasratio is the final value of the feasibility indicator. This
%    indicator converges to 1 for problems with a complementary solution, and
%    to -1 for strongly infeasible problems. If feasratio in somewhere in
%    between, the problem may be nasty (e.g. the optimum is not attained),
%    if the problem is NOT purely linear (see below). Otherwise, the reason
%    must lie in numerical problems: try to rescale the problem.
%
% >  [X,Y,INFO] = SEDUMI(A,b,0) or SEDUMI(A,b) solves the feasibility problem
%    FIND x>=0 such that A*x = b
%
% >  [X,Y,INFO] = SEDUMI(A,0,c) or SEDUMI(A,c) solves the feasibility problem
%    FIND y such that A'*y <= c
%
% >  [X,Y,INFO] = SEDUMI(A,b,c,K) instead of the constraint "x>=0", this
%      restricts x to a self-dual homogeneous cone that you describe in the
%      structure K. Up to 5 fields can be used, called K.f, K.l, K.q, K.r and
%      K.s, for Free, Linear, Quadratic, Rotated quadratic and Semi-definite.
%      In addition, there are fields K.xcomplex, K.scomplex and K.ycomplex
%      for complex-variables.
%
%    (1) K.f is the number of FREE, i.e. UNRESTRICTED primal components.
%      The dual components are restricted to be zero. E.g. if
%      K.f = 2 then x(1:2) is unrestricted, and z(1:2)=0.
%      These are ALWAYS the first components in x.
%
%    (2) K.l is the number of NONNEGATIVE components. E.g. if K.f=2, K.l=8
%      then x(3:10) >=0.
%
%    (3) K.q lists the dimensions of LORENTZ (quadratic, second-order cone)
%      constraints. E.g. if K.l=10 and K.q = [3 7] then
%          x(11) >= norm(x(12:13)),
%          x(14) >= norm(x(15:20)).
%      These components ALWAYS immediately follow the K.l nonnegative ones.
%      If the entries in A and/or c are COMPLEX, then the x-components in
%      "norm(x(#,#))" take complex-values, whenever that is beneficial.
%       Use K.ycomplex to impose constraints on the imaginary part of A*x.
%
%    (4) K.r lists the dimensions of Rotated LORENTZ
%      constraints. E.g. if K.l=10, K.q = [3 7] and K.r = [4 6], then
%          2*x(21)x(22) >= norm(x(23:24))^2,
%          2*x(25)x(26) >= norm(x(27:30))^2.
%      These components ALWAYS immediately follow the K.q ones.
%      Just as for the K.q-variables, the variables in "norm(x(#,#))" are
%      allowed to be complex, if you provide complex data. Use K.ycomplex
%      to impose constraints on the imaginary part of A*x.
%
%    (5) K.s lists the dimensions of POSITIVE SEMI-DEFINITE (PSD) constraints
%      E.g. if K.l=10, K.q = [3 7] and K.s = [4 3], then
%          mat( x(21:36),4 ) is PSD,
%          mat( x(37:45),3 ) is PSD.
%      These components are ALWAYS the last entries in x.
%
%    (a) K.xcomplex lists the components in f,l,q,r blocks that are allowed
%     to have nonzero imaginary part in the primal. For f,l blocks, these
%    (b) K.scomplex lists the PSD blocks that are Hermitian rather than
%      real symmetric.
%    (c) Use K.ycomplex to impose constraints on the imaginary part of A*x.
%
%    The dual multipliers y have analogous meaning as in the "x>=0" case,
%    except that instead of "c-A'*y>=0" resp. "-A'*y>=0", one should read that
%    c-A'*y resp. -A'*y are in the cone that is described by K.l, K.q and K.s.
%    In the above example, if z = c-A'*y and mat(z(21:36),4) is not symmetric/
%    Hermitian, then positive semi-definiteness reflects the symmetric/
%    Hermitian parts, i.e. Z + Z' is PSD.
%
%    If the model contains COMPLEX data, then you may provide a list
%    K.ycomplex, with the following meaning:
%      y(i) is complex if ismember(i,K.ycomplex)
%      y(i) is real otherwise
%    The equality constraints in the primal are then as follows:
%          A(i,:)*x = b(i)      if imag(b(i)) ~= 0 or ismember(i,K.ycomplex)
%          real(A(i,:)*x) = b(i)  otherwise.
%    Thus, equality constraints on both the real and imaginary part
%    of A(i,:)*x should be listed in the field K.ycomplex.
%
%    You may use EIGK(x,K) and EIGK(c-A'*y,K) to check that x and c-A'*y
%    are in the cone K.
%
% >  [X,Y,INFO] = SEDUMI(A,b,c,K,pars) allows you to override the default
%      parameter settings, using fields in the structure `pars'.
%
%    (1) pars.fid     By default, fid=1. If fid=0, then SeDuMi runs quietly,
%      i.e. no screen output. In general, output is written to a device or
%      file whose handle is fid. Use fopen to assign a fid to a file.
%
%    (2) pars.alg     By default, alg=2. If alg=0, then a first-order wide
%      region algorithm is used, not recommended. If alg=1, then SeDuMi uses
%      the centering-predictor-corrector algorithm with v-linearization.
%      If alg=2, then xz-linearization is used in the corrector, similar
%      to Mehrotra's algorithm. The wide-region centering-predictor-corrector
%      algorithm was proposed in Chapter 7 of
%        J.F. Sturm, Primal-Dual Interior Point Approach to Semidefinite Pro-
%        gramming, TIR 156, Thesis Publishers Amsterdam, 1997.
%
%    (3) pars.theta, pars.beta   By default, theta=0.25 and beta=0.5. These
%      are the wide region and neighborhood parameters. Valid choices are
%      0 < theta <= 1 and 0 < beta < 1. Setting theta=1 restricts the iterates
%      to follow the central path in an N_2(beta)-neighborhood.
%
%    (4) pars.stepdif, pars.w. By default, stepdif = 2 and w=[1 1]. 
%       This implements an adaptive heuristic to control ste-differentiation.
%       You can enable primal/dual step length differentiation by setting stepdif=1 or 0.
%      If so, it weights the rel. primal, dual and gap residuals as
%      w(1):w(2):1 in order to find the optimal step differentiation.
%
%    (5) pars.eps     The desired accuracy. Setting pars.eps=0 lets SeDuMi run
%      as long as it can make progress. By default eps=1e-8.
%
%    (6) pars.bigeps  In case the desired accuracy pars.eps cannot be achieved,
%     the solution is tagged as info.numerr=1 if it is accurate to pars.bigeps,
%     otherwise it yields info.numerr=2.
%
%    (7) pars.maxiter Maximum number of iterations, before termination.
%
%    (8) pars.stopat  SeDuMi enters debugging mode at the iterations specified in this vector.
%
%    (9) pars.cg      Various parameters for controling the Preconditioned conjugate
%     gradient method (CG), which is only used if results from Cholesky are inaccurate.
%    (a) cg.maxiter   Maximum number of CG-iterates (per solve). Theoretically needed
%          is |add|+2*|skip|, the number of added and skipped pivots in Cholesky.
%          (Default 49.)
%    (b) cg.restol    Terminates if residual is a "cg.restol" fraction of duality gap.
%          Should be smaller than 1 in order to make progress (default 5E-3).
%    (c) cg.refine    Number of refinement loops that are allowed. The maximum number
%          of actual CG-steps will thus be 1+(1+cg.refine)*cg.maxiter. (default 1)
%    (d) cg.stagtol  Terminates if relative function progress less than stagtol (5E-14).
%    (e) cg.qprec    Stores cg-iterates in quadruple precision if qprec=1 (default 0).
%
%    (10) pars.chol   Various parameters for controling the Cholesky solve.
%     Subfields of the structure pars.chol are:
%    (a) chol.canceltol: Rel. tolerance for detecting cancelation during Cholesky (1E-12)
%    (b) chol.maxu:   Adds to diagonal if max(abs(L(:,j))) > chol.maxu otherwise (5E5).
%    (c) chol.abstol: Skips pivots falling below abstol (1e-20).
%    (d) chol.maxuden: pivots in dense-column factorization so that these factors
%      satisfy max(abs(Lk)) <= maxuden (default 5E2).
%
%    (11) pars.vplot  If this field is 1, then SeDuMi produces a fancy
%      v-plot, for research purposes. Default: vplot = 0.
% 
%    (12) pars.errors  If this field is 1 then SeDuMi outputs some error
%    measures as defined in the Seventh DIMACS Challenge. For more details
%    see the User Guide.
%
% Bug reports can be submitted at http://sedumi.mcmaster.ca.
%
% See also mat, vec, cellK, eyeK, eigK

% This file is part of SeDuMi 1.21 by Imre Polik and Oleksandr Romanko
% Copyright (C) 2009 Imre Polik, Lehigh University, Bethlehem, PA, USA (since 1.21)
%
% Copyright (C) 2006 McMaster University, Hamilton, CANADA  (since 1.1)
%
% Copyright (C) 2001 Jos F. Sturm (up to 1.05R5)
%   Dept. Econometrics & O.R., Tilburg University, the Netherlands.
%   Supported by the Netherlands Organization for Scientific Research (NWO).
%
% Affiliation SeDuMi 1.03 and 1.04Beta (2000):
%   Dept. Quantitative Economics, Maastricht University, the Netherlands.
%
% Affiliations up to SeDuMi 1.02 (AUG1998):
%   CRL, McMaster University, Canada.
%   Supported by the Netherlands Organization for Scientific Research (NWO).
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA

% J.F. Sturm, "Using SeDuMi 1.02, a MATLAB toolbox for optimization over
% symmetric cones," Optimization Methods and Software 11-12 (1999) 625-653.
% http://sedumi.mcmaster.ca

cputime0=cputime;
% ************************************************************
% INITIALIZATION
% ************************************************************
% ----------------------------------------
% Check input
% ----------------------------------------
if (nargin < 5)
    pars.fid = 1;
    if nargin < 3
        if nargin < 2
            error('Should have at least (A,b) or (A,c) arguments')
        end
        if length(b) == max(size(A))
            c = b; b = 0;           % given (A,c): LP feasibility problem
        else
            c = 0;                  % given (A,b): LP feasibility problem
        end
    end
    if isstruct(c)
        if nargin == 4
            pars = K;      % given (A,b,K,pars) or (A,c,K,pars)
        end
        K = c;          % given (A,b,K) or (A,c,K): cone feasibility problem
        if length(b) == max(size(A))
            c = b; b = 0;
        else
            c = 0;
        end
    elseif nargin < 4
        K.l = max(size(A));    % given (A,b,c): default to LP problem.
    end
end
% ----------------------------------------
% Bring data in (strict) internal SeDuMi form, check parameters,
% and print welcome.
% ----------------------------------------
[A,b,c,K,prep,origcoeff] = pretransfo(A,b,c,K,pars);
[N,m]=size(A);
if issparse(A)
    if sprank([A;b'])>sprank(A)
        info.pinf=1;
        x=[];
        %A dual improving direction could be computed easily
        y=[A;b']\[zeros(N,1);1];

        warning('SeDuMi:PrimalInfeasible','The problem is primal infeasible. there is no x such that Ax=b.')
        return
    elseif sprank(A)<m
        error('The coefficient matrix is not full row rank.')
    end
else
    if rank([A;b'])>rank(A)
        info.pinf=1;
        x=[];
        y=[A;b']\[zeros(N,1);1];
        warning('SeDuMi:PrimalInfeasible', 'The problem is primal infeasible, there is no x such that Ax=b.')
        return
    elseif rank(A)<m
        error('The coefficient matrix is not full row rank.')
    end
end
if prep.cpx.dim>0
    origcoeff=[];      % No error measures for complex problems.
end
lponly = (K.l == length(c));
pars = checkpars(pars,lponly);
% ----------------------------------------
% Print welcome and statistics of cone-problem
% ----------------------------------------
my_fprintf(pars.fid,'SeDuMi 1.21 by AdvOL, 2005-2008 and Jos F. Sturm, 1998-2003.\n');
switch pars.alg
    case 0
        my_fprintf(pars.fid,'Alg = 0: No corrector, ');
    case  1
        my_fprintf(pars.fid,'Alg = 1: v-corrector, ');
    case 2
        my_fprintf(pars.fid,'Alg = 2: xz-corrector, ');
end
switch pars.stepdif
    case 0
    case 1
        my_fprintf(pars.fid,'Step-Differentiation, ');
    case 2
        my_fprintf(pars.fid,'Adaptive Step-Differentiation, ');
end
my_fprintf(pars.fid,'theta = %5.3f, beta = %5.3f\n',pars.theta,pars.beta);
% --------------------------------------------------
% Print preprocessing information
% --------------------------------------------------
if pars.prep==1
    if isfield(prep,'sdp')
        blockcount=0;
        varcount=0;
        for sdpind=1:length(prep.sdp)
            if prep.sdp{sdpind}(1)==1
                blockcount=blockcount+1;
                varcount=varcount+prep.sdp{sdpind}(2);
            end
        end
        if blockcount>0
            my_fprintf(pars.fid,'Detected %i diagonal SDP block(s) with %i linear variables\n',blockcount,varcount);
        end
    end
    if isfield(prep,'freeblock1') && ~isempty(prep.freeblock1)
        my_fprintf(pars.fid,'Detected %i free variables in the linear part\n',length(prep.freeblock1));
    end
    if isfield(prep,'Kf') && prep.Kf>0
        switch pars.free
            case 0
                my_fprintf(pars.fid,'Split %i free variables\n',prep.Kf);
            case 1
                my_fprintf(pars.fid,'Put %i free variables in a quadratic cone\n',prep.Kf);
        end
    end
end
% --------------------------------------------------
% Remove dense columns (if any)
% --------------------------------------------------
Ablkjc = partitA(A,K.mainblks);
[dense,DAt.denq] = getdense(A,Ablkjc,K,pars);
if ~isempty(dense.cols)
    dense.A = A(dense.cols,:)';
    A(dense.cols,:) = 0.0;
    Ablkjc = partitA(A,K.mainblks);
else
    dense.A = sparse(length(b),0);
end
% ----------------------------------------
% Order constraints from sparse to dense, and find corresponding
% incremental nonzero pattern "Aord.dz" of At*dy in PSD part.
% ----------------------------------------
Aord.lqperm = sortnnz(A,[],Ablkjc(:,3));        % Sparse LP+Lorentz
DAt.q = findblks(A,Ablkjc,2,3,K.qblkstart);     % Lorentz ddotA-part
if ~isempty(DAt.q)
    DAt.q(dense.q,:) = 0.0;
    DAt.q = DAt.q + spones(extractA(A,Ablkjc,1,2,K.mainblks(1),K.mainblks(2)));
    Aord.qperm = sortnnz(DAt.q,[],[]);
else
    Aord.qperm = (1:length(b))';
end
[Aord.sperm, Aord.dz] = incorder(A,Ablkjc(:,3),K.mainblks(3));  % PSD
% ----------------------------------------
% Get nz-pattern of ADA.
% ----------------------------------------
ADA = getsymbada(A,Ablkjc,DAt,K.sblkstart);
% ----------------------------------------
% Ordering and symbolic factorization of ADA.
% ----------------------------------------
L = symbchol(ADA);
% --------------------------------------------------
% Symbolic fwsolve dense cols: L\[dense.A, dense.blkq],
% sparse ordering for dense column factorization
% --------------------------------------------------
symLden = symbcholden(L,dense,DAt);
% ----------------------------------------
% Initial solution
% ----------------------------------------
[d, v,vfrm, y,y0, R] = sdinit(A,b,c,dense,K,pars);
n = length(vfrm.lab);                         % order of K
merit = (sum(R.w) + max(R.sd,0))^2 * y0 / R.b0;  % Merit function
my_fprintf(pars.fid,'eqs m = %g, order n = %g, dim = %g, blocks = %g\n',...
    length(b),n,length(c),1 + length(K.q) + length(K.s));
my_fprintf(pars.fid,'nnz(A) = %d + %d, nnz(ADA) = %d, nnz(L) = %d\n',nnz(A),nnz(dense.A), nnz(ADA), nnz(L.L));
if ~isempty(dense.cols)
    my_fprintf(pars.fid,'Handling %d + %d dense columns.\n',...
        length(dense.cols),length(dense.q));
end
my_fprintf(pars.fid,' it :     b*y       gap    delta  rate   t/tP*  t/tD*   feas cg cg  prec\n');
my_fprintf(pars.fid,'  0 :            %8.2E %5.3f\n',merit,0);
% ----------------------------------------
% Initialize iterative statistics
% ----------------------------------------
STOP = 0;
err.maxb = 0.0;                 % Estimates the need to recompute residuals
iter = 0;
if pars.vplot == 1
    vlist = [vfrm.lab];
    ratelist = [];
end
wr.delta = 0.0;
wr.desc = 1;
%Seems unnecessary
%rate = 1.0;
feasratio = 0.0;
cputime1 = cputime;
% ************************************************************
% MAIN PREDICTOR-CORRECTOR LOOP
% ************************************************************
while STOP == 0
    iter = iter+1;
    if any(iter == pars.stopat)
        keyboard
    end
    
    if pars.stepdif==2 && ...
            (iter>20 || (iter>1 && (err.kcg + Lsd.kcg>3)) || ...
            (iter>5 && abs(1-feasratio)<0.05) )
        pars.stepdif=1;
    end
    % --------------------------------------------------
    % Compute ADA
    % --------------------------------------------------
    DAt = getDAtm(A, Ablkjc, dense, DAt.denq, d, K);
    ADA = getada1(ADA, A, Ablkjc(:,3), Aord.lqperm, d, K.qblkstart);
    ADA = getada2(ADA, DAt, Aord, K);
    [ADA,absd] = getada3(ADA, A, Ablkjc(:,3), Aord, invcholfac(d.u, K, d.perm), K);
    % ------------------------------------------------------------
    % Block Sparse Cholesky: ADA(L.perm,L.perm) = L.L*diag(L.d)*L.L'
    % ------------------------------------------------------------
    [L.L,L.d,L.skip,L.add] = blkchol(L,ADA,pars.chol,absd);
    % ------------------------------------------------------------
    % Factor dense columns
    % ------------------------------------------------------------
    [Lden, L.d] = deninfac(symLden, L,dense,DAt,d,absd, K.qblkstart,pars.chol);
    % ----------------------------------------
    % FACTORIZATION of self-dual embedding
    % ----------------------------------------
    Lsd = sdfactor(L,Lden, dense,DAt, d,v,y, A,c,K,R, y0,pars);
    % ------------------------------------------------------------
    % Compute and take IPM-step
    % from (v,y,v, y0) --> (xscl,y,zscl,y0)
    % ------------------------------------------------------------
    y0Old = y0;
    [xscl,yNxt,zscl,y0Nxt, w,relt, dxmdz,err, wr] = ...
        wregion(L,Lden,Lsd,...
        d,v,vfrm,A,DAt,dense, R,K,y,y0,b, pars, wr);
    % ------------------------------------------------------------
    % Evaluate the computed step.
    % ------------------------------------------------------------
    if y0Nxt > 0
        R.b = R.b + err.b / y0Nxt;
        R.sd = R.sd + err.g / y0Nxt;
        R.b0 = R.b0 + err.db0 / y0Nxt;
        y0 = y0Nxt;
    else
        R.b = (y0Nxt * R.b + err.b) / y0Old;      % In fact, we should have y0=0.
        R.sd = (y0Nxt * R.sd + err.g) / y0Old;
        R.b0 = (y0Nxt * R.b0 + err.db0) / y0Old;
        R.w(2) = abs(y0Nxt/y0Old) * R.w(2);        %=0: dual feasible
        R.c = (y0Nxt/y0Old) * R.c;
        R.maxRc = norm(R.c,inf);
        y0 = y0Old;
    end
    R.maxRb = norm(R.b,inf);                % Primal residual
    R.w(1) = 2 * pars.w(1) * R.maxRb / (1+R.maxb);
    meritOld = merit;
    merit = (sum(R.w) + max(R.sd,0))^2 * y0 / R.b0;
    rate = merit / meritOld;
    if (rate >= 0.9999) && (wr.desc == 1)
        % ------------------------------------------------------------
        % STOP = -1  --> Stop due to numerical problems
        % ------------------------------------------------------------
        STOP = -1;                  % insuf. progress in descent direction.
        iter = iter - 1;
        y0 = y0Old; %#ok
        my_fprintf(pars.fid,'Run into numerical problems.\n');
        break
    end
    feasratio = dxmdz(1) / v(1);            % (deltax0/x0) - (deltaz0/z0)
    % --------------------------------------------------
    % Primal-Dual transformation
    % --------------------------------------------------
    y = yNxt;
    by = full(sum(b.*y));
    [d,vfrm] = updtransfo(xscl,zscl,w,d,K);
    v = frameit(vfrm.lab,vfrm.q,vfrm.s,K);
    x0 = sqrt(d.l(1)) * v(1);
    % ----------------------------------------
    % SHOW ITERATION STATISTICS
    % ----------------------------------------
    my_fprintf(pars.fid,' %2.0f : %10.2E %8.2E %5.3f %6.4f %6.4f %6.4f %6.2f %2d %2d  ',...
        iter,by/x0,merit,wr.delta,rate,relt.p,relt.d,feasratio,err.kcg, Lsd.kcg);
    if pars.vplot == 1
        vlist = [vlist vfrm.lab/sqrt((R.b0*y0)/n)];
        ratelist = [ratelist rate];
    end
    % ----------------------------------------
    % If we get in superlinear region of LP,
    % try to guess optimal solution:
    % ----------------------------------------
    if lponly && (rate < 0.05)
        [xsol,ysol] = optstep(A,b,c, y0,y,d,v,dxmdz, ...
            K,L,symLden,dense, Ablkjc,Aord,ADA,DAt, feasratio, R,pars);
        if ~isempty(xsol)
            STOP = 2;                   % Means that we guessed right !!
            feasratio = 1 - 2*(xsol(1)==0);
            break
        end
    elseif (by > 0) && (abs(1+feasratio) < 0.05) && (R.b0*y0 < 0.5)
        if max(eigK(full(qreshape(Amul(A,dense,y,1),1,K)),K)) <= pars.eps * by
            STOP = 3;                   % Means Farkas solution found !
            break
        end
    end
    % --------------------------------------------------
    % OPTIMALITY CHECK: stop if y0*resid < eps * (x0+z0).
    % For feas. probs, we should divide the residual by x0, otherwise by z0.
    % Before stopping, recompute R.norm, since it may have changed due to
    % residual updates (the change should be small though).
    % --------------------------------------------------
    r0 = sum(R.w);
    cx = by + y0*R.sd - x0 / d.l(1);
    rgap = max(cx-by,0) / max([abs(cx),abs(by),1e-3 * x0]);
    precision1=y0*r0/(1+x0);
    precision2=(y0 * r0 + rgap)/x0;
    my_fprintf(pars.fid,'%1.1E\n',max(precision1,precision2));
    if precision1 < pars.eps       % P/D residuals small
        if precision2 < pars.eps    %Approx feasible and optimal
            STOP = 1;
            break
        elseif y0 * R.maxRb + x0 * R.maxb < -pars.eps * cx   % Approx Farkas
            STOP = 1;
            break
        elseif y0 * R.maxRc + x0 * R.maxc < pars.eps * by    % Approx Farkas
            STOP = 1;
            break;
        end
    end
    if iter >= pars.maxiter
        my_fprintf(pars.fid,'Maximum number of iterations reached.\n');
        STOP = -1;
    end
end % while STOP == 0.
my_fprintf(pars.fid,'\n');
clear ADA 
nnzLadd=nnz(L.add);
nnzLskip=nnz(L.skip);
normLL=full(max(max(abs(L.L))));
clear L
% ************************************************************
% FINAL TASKS:
% ************************************************************
cputime2=cputime;
info.iter = iter;
info.feasratio = feasratio;
info.pinf = 0; info.dinf = 0;
info.numerr = 0;
% ------------------------------------------------------------
% Create x = D*v.
% ------------------------------------------------------------
if STOP == 2                % Exact optimal solution found (LP)
    x = xsol;
    y = ysol;
elseif STOP == 3            % Farkas solution y found (in early stage)
    x = zeros(length(c),1);
else
    x = [sqrt(d.l).*v(1:K.l); asmDxq(d,v,K); psdscale(d,v,K,1)];
end
% --------------------------------------------------
% Compute cx, Ax, etc.
% --------------------------------------------------
x0 = x(1);
cx = full(sum(c.*x)); 
abscx = sum(abs(c).*abs(x));
by = full(sum(b.*y));
Ax = Amul(A,dense,x,0);
Ay = full(Amul(A,dense,y,1));      % "full" since y may be scalar.
normy = norm(y);
normx = norm(x(2:end));
clear A

% ------------------------------------------------------------
% Determine infeasibility
% ------------------------------------------------------------

pinf = norm(x0*b-Ax);
z = qreshape(Ay-x0*c,1,K);
dinf = max(eigK(z,K));
if x0 > 0
    relinf = max(pinf / (1+R.maxb), dinf / (1+R.maxc)) / x0;
    % ------------------------------------------------------------
    % If infeasibility larger than epsilon, evaluate Farkas-infeasibility
    % ------------------------------------------------------------
    if relinf > pars.eps
        pdirinf = norm(Ax);
        ddirinf = max(eigK(qreshape(Ay,1,K),K));
        if cx < 0.0
            reldirinf = pdirinf / (-cx);
        else
            reldirinf = inf;
        end
        if by > 0.0
            reldirinf = min(reldirinf, ddirinf / by);
        end
        % ------------------------------------------------------------
        % If the quality of the Farkas solution is good and better than
        % the approx. feasible soln, set x0=0: Farkas solution found.
        % ------------------------------------------------------------
        if (reldirinf < pars.eps) || (relinf > max(pars.bigeps, reldirinf))
            x0 = 0.0;
            pinf = pdirinf;
            dinf = ddirinf;
        end
    end % relinf too large
end % x0 > 0
% ------------------------------------------------------------
% Interpret the solution as feasible:
% ------------------------------------------------------------
info.r0 = Inf;
if x0 > 0
    x = x / x0;
    y = y / x0;
    pinf = pinf /x0;
    dinf = dinf / x0;
    cx = cx/ x0;
    by = by / x0;
    normx = normx / x0;
    normy = normy / x0;
    if cx <= by                % zero or negative duality gap
        r0 = 0;
    elseif cx == 0.0           % Dual feasibility problem
        r0 = -by/(R.maxb*normy +1E-10 * x0);
    elseif by == 0.0           % Primal feasibility problem
        r0 = cx / (R.maxc*normx +1E-10 * x0);
    else                       % Optimization problem
        r0 = (cx-by)/(abs(by) + 1E-5 * (x0+abscx));
    end
    if r0 == 0,
        sigdig = Inf;
    else
        sigdig = -log10(r0);
    end
    my_fprintf(pars.fid,...
        'iter seconds digits       c*x               b*y\n');
    my_fprintf(pars.fid,'%3d %8.1f %5.1f %- 17.10e %- 17.10e\n',...
        iter,cputime2-cputime1,sigdig,cx,by);
    my_fprintf(pars.fid,'|Ax-b| = %9.1e, [Ay-c]_+ = %9.1E, |x|=%9.1e, |y|=%9.1e\n',...
        pinf,dinf,normx,normy);
    % ------------------------------------------------------------
    % Determine level of numerical problems with x0>0 (feasible)
    % ------------------------------------------------------------
    info.r0 = max([r0 ; pinf;dinf] ./ [1; 1+R.maxb+(1E-3)*R.maxRb;...
            1+R.maxc+(1E-3)*R.maxRc]);
    if STOP == -1
        if info.r0 > pars.bigeps
            my_fprintf(pars.fid, 'No sensible solution found.\n');
            info.numerr = 2;                          % serious numerical error
        elseif info.r0 > pars.eps
            info.numerr = 1;                          % moderate numerical error
        else
            info.numerr = 0;                          % achieved desired accuracy
        end
    else
        info.r0 = min( info.r0, pars.eps );
    end
else  % (if x0>0)
    % --------------------------------------------------
    % Infeasible problems: pinf==norm(Ax), dinf==max(eigK(At*y,K)).
    % --------------------------------------------------
    if pinf < -pars.bigeps * cx
        info.r0 = abs(pinf/cx);
        info.dinf = 1;
        abscx = -cx;
        pinf = pinf / abscx;
        normx = normx / abscx;
        x = x / abscx;
        my_fprintf(pars.fid, 'Dual infeasible, primal improving direction found.\n');
    end
    if dinf < pars.bigeps * by
        info.r0 = abs(dinf/by);
        info.pinf = 1;
        dinf = dinf / by;
        normy = normy / by;
        y = y / by;
        my_fprintf(pars.fid, 'Primal infeasible, dual improving direction found.\n');
    end
    my_fprintf(pars.fid,'iter seconds  |Ax|    [Ay]_+     |x|       |y|\n');
    my_fprintf(pars.fid,'%3d %8.1f %9.1e %9.1e %9.1e %9.1e\n',...
        iter,cputime2-cputime1,pinf,dinf,normx,normy);
    % --------------------------------------------------
    % Guess infeasible, but stopped due to numerical problems
    % --------------------------------------------------
    if info.pinf + info.dinf == 0
        my_fprintf(pars.fid, 'Failed: no sensible solution/direction found.\n');
        info.numerr = 2;
    elseif STOP == -1
        if (pinf > -pars.eps * cx) && (dinf > pars.eps * by)
            info.numerr = 1;
        else
            info.numerr = 0;
        end
    end
end
% ----------------------------------------
% - Bring xsol into the complex format of original (At,c),
% - transform q-variables into r-variables (Lorentz),
% - bring ysol into complex format, indicated by K.ycomplex.
% - at 0's in ysol where rows where removed.
% ----------------------------------------'
[x,y,K] = posttransfo(x,y,prep,K,pars); %#ok
% Detailed timing
%Preprocessing+IPM+Postprocessing
info.timing=[cputime1-cputime0 cputime2-cputime1 cputime-cputime2];
% Total time (for backward compatibility)
info.cpusec=sum(info.timing);
my_fprintf(pars.fid,'\nDetailed timing (sec)\n')
my_fprintf(pars.fid,'   Pre          IPM          Post\n')
my_fprintf(pars.fid,'%1.3E    ',info.timing)
my_fprintf(pars.fid,'\n')


% ----------------------------------------
% Make a fancy v-plot if desired
% ----------------------------------------
if pars.vplot == 1
    subplot(2,1,1)
    plot(0:iter,vlist,'o',[0 iter],[1 1],'b',...
        [0 iter],[pars.theta pars.theta],'g')
    title('Wide region v-plot')
    xlabel('iterations')
    ylabel('normalized v-values')
    subplot(2,1,2)
    plot(1:iter,ratelist)
    axis([0 iter 0 1])
    title('Reduction rates')
    xlabel('iterations')
    ylabel('reduction rate')
end
my_fprintf(pars.fid,'Max-norms: ||b||=%d, ||c|| = %d,\n',R.maxb,R.maxc);
my_fprintf(pars.fid,'Cholesky |add|=%d, |skip| = %d, ||L.L|| = %g.\n',...
    nnzLadd, nnzLskip, normLL);

% ----------------------------------------
% Compute error measures if needed
% ----------------------------------------
if ~isempty(origcoeff)
    %Reload the original coefficients
    s=(origcoeff.c)-(origcoeff.At)*sparse(y);       %To make s sparse
    cx=sum((origcoeff.c).*x);        %faster than c'*x
    by=sum((origcoeff.b).*y);
    xs=sum(x.*s);
    normb=norm(origcoeff.b,1);
    normc=norm(origcoeff.c,1);
    info.err=zeros(1,6);
    %     Error measures.
    %     Primal infeasibility
    info.err(1)=norm(x'*(origcoeff.At)-(origcoeff.b)',2)/(1+normb);
    %Let us get rid of the K.f part, since the free variables don't make
    %any difference in the cone infeasibility.
    %origcoeff.K.f=0;
    %     Primal cone infeasibility
    info.err(2)=max(0,-min(eigK(full(x(origcoeff.K.f+1:end)),origcoeff.K))/(1+normb));
    %     Dual infeasibility
    %info.err(3)=0.0; %s is not maintained explicitely
    %     Dual cone infeasibility
    info.err(4)=max(0,-min(eigK(full(s(origcoeff.K.f+1:end)),origcoeff.K))/(1+normc));
    %     Relative duality gap
    info.err(5)=(cx-by)/(1+abs(cx)+abs(by));
    %     Relative complementarity
    info.err(6)=xs/(1+abs(cx)+abs(by));

    my_fprintf(pars.fid,'\nDIMACS error measures\n')
    my_fprintf(pars.fid,'  PInf     PConInf     DInf    DConInf    RelGap    RelComp\n')
    my_fprintf(pars.fid,'%2.2E  ',info.err)
    my_fprintf(pars.fid,'\n')
end
