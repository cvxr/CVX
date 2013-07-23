function [y,k, DAy] = loopPcg(L,Lden,At,dense,d, DAt,K, b,p,ssqrNew,cgpars, restol)
% [y,k, DAy] = loopPcg(L,Lden,At,dense,d, DAt,K, b,p,ssqrNew,cgpars, restol)
%
% LOOPPCG Solve y from AP(d)A' * y = b
% using PCG-method and Cholesky L as conditioner.
% If L is sufficiently accurate, then only 1 CG-step is needed.
% It assumes that the previous step was p, with
% ssqrNew = bOld'*inv(L*THETA*L')*bOld, and bOld the residual before
% the step p was taken. If p = [], then the PCG is started from scratch.
%
%  k = #(CG-iterations).
%  DAy = D*A'*y with D'*D = P(d).
%
% Warning: if the scaling operation P(d) gets ill-conditioned, the
% precision in y may be insufficient to compute DAy satisfactory.
% In this case, one should allow refinement, see wrapPcg.
%
% **********  INTERNAL FUNCTION OF SEDUMI **********
%
% See also sedumi, wrapPcg

% This file is part of SeDuMi 1.1 by Imre Polik and Oleksandr Romanko
% Copyright (C) 2005 McMaster University, Hamilton, CANADA  (since 1.1)
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
%

% --------------------------------------------------
% INITIALIZE:
%  y=0 with r = b.
% --------------------------------------------------
k = 0;
STOP = 0;
r = b;
finew = 0;
y = [];        % means all-0
normrmin = norm(r,inf);
ymin = [];
% --------------------------------------------------
% Conjugate Gradient until convergence
% --------------------------------------------------
while STOP == 0
    % --------------------------------------------------
    % P r e - C o n d i t i o n i n g:
    % Solve L*Lr = r, L*THETA*tmp = r.
    % p=[] ==> initialize p to solve L*THETA*L'*p = r.
    % --------------------------------------------------
    Lr = fwdpr1(Lden,sparfwslv(L, r));
    tmp = Lr ./ L.d;
    if isempty(p)
        ssqrNew = Lr'*tmp;
        p = sparbwslv(L, bwdpr1(Lden,tmp));
    else
        % --------------------------------------------------
        % General iterate: make p conjugate to previous iterate(s):
        % --------------------------------------------------
        ssqrOld = ssqrNew;
        ssqrNew = Lr'*tmp;
        p = (ssqrNew/ssqrOld) * p;
        p = p + sparbwslv(L, bwdpr1(Lden,tmp));
    end
    % --------------------------------------------------
    % SCALING OPERATION AND MATRIX*VECTOR.
    % Let DDAp = P(d)*A'*p and ssqrDAp = ||P(d)^{1/2}*A'*p||^2.
    % Set alpha = ssqrNew / ssqrDAp
    % --------------------------------------------------
    Ap = vecsym(Amul(At,dense,p,1), K);
    [DDAp, DApq, DAps, ssqrDAp] = PopK(d,Ap,K);
    if ssqrDAp > 0.0
        k = k + 1;
        % --------------------------------------------------
        % Take step:  y := y + alpha*p
        %--------------------------------------------------
        alpha = ssqrNew / ssqrDAp;
        if ~isempty(y)
            if isstruct(y)
                [y.hi,y.lo] = quadadd(y.hi,y.lo,alpha*p);
            else
                y = y + alpha * p;
            end
        elseif cgpars.qprec > 0
            y.hi = alpha * p; y.lo = zeros(length(p),1);
        else
            y = alpha * p;
        end
        % --------------------------------------------------
        % Update residual r := r - alpha * A*[P(d)*Ap]. MATRIX*VECTOR.
        % --------------------------------------------------
        tmp = Amul(At,dense,DDAp)+ DAt.q'*DApq(:,1) +...
            DAt.denq*DApq(dense.q,1);
        r = r - alpha * tmp;
        % --------------------------------------------------
        % Convergence check (HEURISTIC)
        % --------------------------------------------------
        fiprev = finew;
        if isstruct(y)
            finew = (b+r)'*y.hi + (b+r)'*y.lo;
        else
            finew = (b+r)'*y;
        end
        normr = norm(r,inf);
        if normr < normrmin
            ymin = y;
            normrmin = normr;
        end
        if normr < restol
            STOP = 1;
        elseif (finew-fiprev < cgpars.stagtol * fiprev)
            STOP = 2;
        elseif k >= cgpars.maxiter
            STOP = 2;
        end
    else
        %     my_fprintf('Warning: DAp = 0 in PCG\n');
        STOP = 1;         % If DAp == 0 then can't go on.
    end
end
% --------------------------------------------------
% OUTPUT: return projection D * At*y on request.
% --------------------------------------------------
if STOP == 2
    y = ymin;      % Take best so far
end
if isempty(y)
    DAy = [];
    return
end
if nargout >= 3
    if k == 1
        DAy = alpha*[sqrt(d.l).*Ap(1:K.l); asmDxq(d,Ap,K,DApq); DAps];
    else
        if isstruct(y)
            Ap = vecsym(Amul(At,dense,y.hi,1), K);
        else
            Ap = vecsym(Amul(At,dense,y,1), K);
        end
        DAy = [sqrt(d.l).*Ap(1:K.l); asmDxq(d,Ap,K); psdscale(d,Ap,K)];
        if isstruct(y)
            Ap = vecsym(Amul(At,dense,y.lo,1), K);
            DAy = DAy + [sqrt(d.l).*Ap(1:K.l); asmDxq(d,Ap,K); psdscale(d,Ap,K)];
        end
    end
end
if isstruct(y)
    y = y.hi;
end