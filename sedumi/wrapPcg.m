function [y,dx,k, r] = wrapPcg(L,Lden,At,dense,d, DAt,K, rb,rv,cgpars, y0)
% [y,r,k, DAy] = normeqPcg(L,Lden,At,dense,d, DAt,K, b, cgpars, y0,rhs)
%
% WRAPPCG  Solve y from AP(d)A' * y = b
% using PCG-method and Cholesky L as conditioner.
% If L is sufficiently accurate, then only 1 CG-step is needed.
% In general, proceeds until ||DAy - DA'(AD^2A')^{-1}b||
% has converged (to zero). k = #(CG-iterations).
%
% **********  INTERNAL FUNCTION OF SEDUMI **********
%
% See also sedumi, loopPcg

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

% --------------------------------------------------
% INITIALIZE:
% --------------------------------------------------
restol = y0 * cgpars.restol;
dx = [sqrt(d.l).*rv(1:K.l); asmDxq(d,rv,K); psdscale(d,rv,K,1)];
r = Amul(At,dense,dx);
if ~isempty(rb)
    r = r+rb;
end
% --------------------------------------------------
% Pre - Conditioning/ Direct Solving:
% Solve L*THETA*L'*p = r, set ssqrNew = r'*inv(L*THETA*L')*r.
% --------------------------------------------------
p = fwdpr1(Lden,sparfwslv(L, r));
y = p ./ L.d;
ssqrNew = p'*y;
p = sparbwslv(L, bwdpr1(Lden,y));
% --------------------------------------------------
% SCALING OPERATION AND MATRIX*VECTOR.
% Let dx = D*A'*p
% Set alpha = ssqrNew / ||dx||^2
% --------------------------------------------------
x = vecsym(Amul(At,dense,p,1), K);
dx = [sqrt(d.l).*x(1:K.l); asmDxq(d,x,K); psdscale(d,x,K)];
ssqrdx = norm(dx)^2;
if ssqrdx <= 0.0
    y = zeros(length(r),1);
    k = 0;
    dx = rv;
    return
end
% --------------------------------------------------
% Take 1st step:  y *= alpha, alpha:=ssqrNew / ssqrdx,
% dx = rv - alpha * dx.
% --------------------------------------------------
k = 1;
alpha = ssqrNew / ssqrdx;
y = alpha * p;
dx = rv - alpha * dx;
% --------------------------------------------------
% Compute 1st residual r := b + A*D*dx.   (MATRIX*VECTOR).
% --------------------------------------------------
x = [sqrt(d.l).*dx(1:K.l); asmDxq(d,dx,K); psdscale(d,dx,K,1)];
r = Amul(At,dense,x);
if ~isempty(rb)
    r = r + rb;
end
normr = norm(r,inf);
if normr < restol
    return                % No CG nor refinement needed.
end
STOP = 0;
trial = 0;
k = 1;
% --------------------------------------------------
% Conjugate Gradient until convergence
% --------------------------------------------------
while ~STOP
    [dy,dk, x] = loopPcg(L,Lden,At,dense,d, DAt,K, r,p,ssqrNew,...
        cgpars, restol);
    if isempty(dy)
        return
    else
        % ------------------------------------------------------------
        % Add step from PCG-method. If P(d) is highly ill-conditioned, then
        % "x=DA*dy" as returned from PCG is inaccurate, and the residual
        % based on x can differ from the internal PCG residual. Therefore,
        % we recompute r based on x, and refine if needed & allowed.
        % ------------------------------------------------------------
        k = k + dk;
        y = y + dy;
        dx = dx - x;
        x = [sqrt(d.l).*dx(1:K.l); asmDxq(d,dx,K); psdscale(d,dx,K,1)];
        r = Amul(At,dense,x);
        if ~isempty(rb)
            r = r+rb;
        end
        normr = norm(r,inf);
        if normr < restol
            return
        elseif trial >= cgpars.refine
            return
        else
            p = [];
            trial = trial + 1;              % Refine
        end
    end
end