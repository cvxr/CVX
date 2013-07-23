function [x,y] = optstep(A,b,c, y0,y,d,v,dxmdz, K,L,symLden,...
    dense,Ablkjc,Aord,ADA,DAt, feasratio, R,pars)
% [x,y] = optstep(A,b,c, y0,y,d,v,dxmdz, K,L,symLden,dense,Ablkjc,Aord,...
%                 ADA,DAt, feasratio, R,pars)
%
% OPTSTEP Implements Mehrotra-Ye type optimality projection for
%  IPM-LP solver.
%
% **********  INTERNAL FUNCTION OF SEDUMI **********
%
% See also sedumi

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

% ------------------------------------------------------------
% For LP-problems, we PROJECT onto OPTIMAL FACE (x0 > 0)  *OR*
%   onto a direction (x0 = 0) with the same sign for c'*x and b'*y.
% ------------------------------------------------------------
if abs(abs(feasratio)-1) < 0.1
    x0 = sqrt(d.l(1)) * v(1);
    z0 = x0 / d.l(1);
    %This value is never used.
    %deptol = 1E-10 * max(x0,z0);
    if (feasratio < -0.5) && (x0 < z0*z0)
        x0 = 0;                          % Try project onto direction.
    end
    % ------------------------------------------------------------
    % Set d(Non-basic-LP) = 0
    % Nonbasic is where dx < dz i.e. dxmdz < 0, ideally: dx_N = -v_N, dz_N = 0.
    % ------------------------------------------------------------
    lpNB = find(dxmdz < 0);
    d.l(lpNB) = 0;
    % --------------------------------------------------
    % Compute ADA, with d[N] = 0. Hence A[B]*D[B]^2*A[B]'.
    % --------------------------------------------------
    DAt = getDAtm(A,Ablkjc,dense,DAt.denq,d,K);
    if sum(K.s)==0
        %ADA is global already
        absd=getada(A,K,d,DAt);
    else
        ADA = getada1(ADA, A, Ablkjc(:,3), Aord.lqperm, d, K.qblkstart);
        ADA = getada2(ADA, DAt, Aord, K);
        [ADA,absd] = getada3(ADA, A, Ablkjc(:,3), Aord, invcholfac(d.u, K, d.perm), K);
    end
    % ------------------------------------------------------------
    % Block Sparse Cholesky: ADA(L.perm,L.perm) = L.L*diag(L.d)*L.L'
    % ------------------------------------------------------------
    [L.L,L.d,L.skip,L.add] = blkchol(L,ADA,pars.chol,absd);
    clear ADA;
    % ------------------------------------------------------------
    % Factor dense columns
    % ------------------------------------------------------------
    [Lden,L.d] = deninfac(symLden, L,dense,DAt,d,absd,K.qblkstart,pars.chol);
    % ------------------------------------------------------------
    % Solve ADAt*psi = -x0*b+A*D*v, dx = v-D*At*psi.  LEAST SQUARES.
    % ------------------------------------------------------------
    [psi,dx,err.kcg,err.b] = wrapPcg(L,Lden,A,dense,d, DAt,K,...
        (-x0) * b,v, pars.cg,pars.eps / pars.cg.restol); %#ok
    x = sqrt(d.l) .* dx;
    % ----------------------------------------
    % CHECK WHETHER x[B] >= 0 AND WHETHER RESIDUAL DID NOT DETERIORATE.
    % ----------------------------------------
    if (min(x) < 0.0) || (norm(err.b,inf) > 2 * max(max(y0,1e-10 * x0) * R.maxb, y0 * R.maxRb))
        x = [];  % Incorrect guess of LP-basis
        return
    else
        % ==================== DUAL PART (LP ONLY) ====================
        % ------------------------------------------------------------
        % Solve A[B]'*dy = x0*c[B]-A[B]'*y, so that zB=0. We solve
        % the equivalent system obtained after pre-multiplying with A[B]*P(d[b]).
        % ------------------------------------------------------------
        rhs = sqrt(d.l) .* (x0*c - Amul(A,dense,y,1));
        % ------------------------------------------------------------
        % Solve ADA*dy = A*D*rhs.   THIS IS DEBATABLE !
        % ------------------------------------------------------------
        dy = wrapPcg(L,Lden,A,dense,d, DAt,K,...
            zeros(length(b),1),rhs, pars.cg,pars.eps / pars.cg.restol);
        y = y+dy;
        % ------------------------------------------------------------
        % CHECK WHETHER Z[N] >= 0 AND RESID ON Z[B]=0 DID NOT DETERIORATE.
        % ALSO, we need z0 >= 0 and x0+z0 > 0.
        % ------------------------------------------------------------
        z = x0*c - Amul(A,dense,y,1);
        z(1) = 0;
        zB = z;
        zB(lpNB) = 0;
        normzB = norm(zB,inf);
        cx = c'*x;
        by = b'*y;
        z0 = by - cx; %[JFS 9/2003: changed condition below
        if (~isempty(lpNB) && (min(z(lpNB)) < 0.0)) || ...
                normzB > 5 * max(1E-10 * (x0+(x0==0)) * norm(c), min(y0,1e-8) * norm(R.c))
            x = [];  % Incorrect guess of LP-basis
            return
        end
        if x0 == 0
            if z0 <= 0
                x = [];  % z0 <= 0 --> P/D Direction not improving anymore
                return
            end
        elseif z0 < - (5E-8) * (1+abs(by) +max(abs(b)))
            x = [];        % x0>0, z0 << 0 then not optimal.
            return
        end
    end
else
    % ------------------------------------------------------------
    % If optimality step impossible , then return emptyset.
    % ------------------------------------------------------------
    x = []; y = [];
end