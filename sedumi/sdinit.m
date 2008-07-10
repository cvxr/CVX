%                               [d, v,vfrm,y,y0, R] = sdinit(At,b,c,dense,K,pars)
% SDINIT  Initialize with identity solution, for self-dual model.
%
% **********  INTERNAL FUNCTION OF SEDUMI **********
%
% See also sedumi

function [d, v,vfrm,y,y0, R] = sdinit(At,b,c,dense,K,pars)
%
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
% Initialize m = #eqns, n = order(K).
% --------------------------------------------------
m = length(b);
n = K.l + 2 * length(K.q) + K.rLen + K.hLen;
% ------------------------------------------------------------
% Include some inf-norms to determine scaling/magnitude of data.
% ------------------------------------------------------------
R.maxb = norm(b,inf);
R.maxc = norm(c,inf);
% ----------------------------------------
% Choose trivial soln on central path:
% y = 0;  v = mu * identity with mu := pars.mu * sqrt(norm(b)*norm(c)).
% ----------------------------------------
y = zeros(m,1);
mu = pars.mu * sqrt((1+R.maxb)*(1+R.maxc));
id = qreshape(eyeK(K),0,K);
v = mu * id;
y0 = n * mu;   % b0 * y0 = norm(v)^2.
R.b0 = mu;
% ------------------------------------------------------------
% Determine primal/dual scaling "d0", so that x = d0*v, z=v/d0.
% sd-var different: x0 = 1, z0 = mu^2.
% ------------------------------------------------------------
d0 = sqrt((1+R.maxb) / (1+R.maxc));
x0 = pars.mu; z0 = mu^2 / x0;
cx = d0 * (c'*v);
R.sd = (z0 + cx)/ y0;
% ----------------------------------------
% Since x=d0^2*z, we have d = d0 * identity
% ----------------------------------------
d.l = d0^2 * ones(K.l,1);                                  % LP
d.l(1) = x0/z0;    % =1/mu^2.
d.det = d0^2 * ones(length(K.q),1);                        % LORENTZ
d.q1 = (sqrt(2) * d0) * ones(length(K.q),1);   % trace part of d0 * identity
d.q2 = zeros(K.mainblks(3)-K.mainblks(2),1);
d.auxdet = sqrt(2*d.det);
d.auxtr = sqrt(2)*(d.q1 + d.auxdet);
d.u = sqrt(d0) * id(K.lq+1:end);        % sqrt(d0) * identity   (PSD)
d.perm = [];
% ----------------------------------------
% Create Jordan frame of v = mu * identity
% vfrm.q denotes Lorentz frame (norm 1/sqrt(2) blocks).
% vfrm.s = I in Householder product form,
% ----------------------------------------
vfrm.lab = mu*ones(n,1);         % identity
vfrm.q = d.q2;
vfrm.s = qrK(d.u,K);
% ----------------------------------------
% Residuals R.b, R.c, R.sd, and size-parameter R.b0.
%         A x - x0 b - y0 R.b      = 0
%  -A'y       + x0 c + y0 R.c - z  = 0
%   b'y - c'x        + y0 R.sd - z0 = 0
%  Rb'y -Rc'x - x0 Rsd             = -b0
% ----------------------------------------
R.b = d0 * Amul(At,dense,v,0);          % x = d0 *v
R.b = (R.b - x0*b) / y0;
R.c = vecsym(v/d0-x0*c,K) / y0; % z = v/d0
R.c(1) = 0.0;               % for artificial (x0,z0)
% ------------------------------------------------------------
% Some inf-norms to determine scaling/magnitude of Residuals.
% Note that if R.sd < 0 then this residual may be absorbed into z0.
% ------------------------------------------------------------
R.maxRb = max(1e-6, norm(R.b,inf));
R.maxRc = max(1e-6, norm(R.c,inf));
R.norm = max([R.maxRb, R.maxRc, R.sd]);
R.w = 2 * pars.w .* [R.maxRb;R.maxRc] ./ [1+R.maxb;1+R.maxc];