function y = psdscale(ud,x,K,transp)
% y = psdscale(ud,x,K [,transp])
%
% PSDSCALE  Computes length lenud (=sum(K.s.^2)) vector y.
%   !transp (default) then y[k] = vec(Ldk' * Xk * Ldk)
%   transp == 1 then y[k] = vec(Udk' * Xk * Udk)
%   Uses pivot ordering ud.perm if available and nonempty.
%
% **********  INTERNAL FUNCTION OF SEDUMI **********
%
% See also scaleK, factorK.

% This file is part of SeDuMi 1.3 by Imre Polik and Oleksandr Romanko
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

%The function is quite dirty as the type and dimension of the inputs may
%change, that is why we have so many subcases.

Ks = K.s;
if isempty(Ks),
    y = [];
    return
end
if nargin < 4,
    transp = false;
end
Kq = Ks .* Ks;
nr = K.rsdpN;
nc = length(Ks);
N  = sum(Kq) + sum(Kq(nr+1:end));
y  = zeros(N,1);
xi = length(x) - N;
yi = 0;
ui = 0;
if isstruct(ud)
    perm = ud.perm;
    if isempty(perm),
        prep = false;
        postp = false;
    else
        prep = ~transp;
        postp = transp;
        pi = 0;
    end
    ud = ud.u;
else
    prep  = false;
    postp = false;
end
for i = 1 : nc,
    ki = Ks(i);
    qi = Kq(i);
    TT = ud(ui+1:ui+qi); ui=ui+qi;
    if i > nr,
        TT = TT + 1i*ud(ui+1:ui+qi); ui=ui+qi;
    end
    TT = reshape(TT,ki,ki);
    if transp,
        TT = triu(TT);
    else
        TT = tril(TT);
    end
    XX = x(xi+1:xi+qi); xi=xi+qi;
    if i > nr,
        XX = XX + 1i*x(xi+1:xi+qi); xi=xi+qi;
    end
    XX = reshape(XX,ki,ki);
    if prep,
        PP = perm(pi+1:pi+ki); pi=pi+ki;
        if any(diff(PP)~=1),
            XX = XX(PP,PP);
        end
    end
    if nnz(XX) < 0.1 * qi,
        XX = sparse(XX);
    end
    XX = TT' * XX * TT;
    if postp,
        PP = perm(pi+1:pi+ki); pi=pi+ki;
        if any(diff(PP)~=1),
            XX(PP,PP) = XX;
        end
    end
    y(yi+1:yi+qi) = real(XX); 
    yi = yi+qi;
    if i > nr,
        XX = imag(XX);
        % Needed, otherwise psdfactor() will sometimes fail.
        XX(1:ki+1:end) = 0;
        y(yi+1:yi+qi) = XX; 
        yi = yi+qi;
    end
end
