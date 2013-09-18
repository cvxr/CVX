function [lab,q] = psdeig(x,K)
% [lab,q] = psdeig(x,K)
%
% PSDEIG  Computes spectral coefficients of x w.r.t. K
%   Arguments "q" is optional - without it's considerably faster.
%   FLOPS indication: 1.3 nk^3 versus 9.0 nk^3 for nk=500,
%                     1.5 nk^3        9.8 nk^3 for nk=50.
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

Ks = K.s;
if isempty(Ks),
    lab = [];
    return
end
Kq  = Ks .* Ks;
nr  = K.rsdpN;
nc  = length(Ks);
N   = sum(Kq) + sum(Kq(nr+1:end));
xi  = length(x) - N;
ei  = 0;
lab = zeros(sum(Ks),1);
needv = nargout > 1;
if needv,
    q = zeros(N,1);
    vi = 0;
end
for i = 1 : nc,
    ki = Ks(i);
    qi = Kq(i);
    XX = x(xi+1:xi+qi); 
    xi = xi+qi;
    if i > nr,
        XX = XX + 1i*x(xi+1:xi+qi); 
        xi = xi+qi;
    end
    XX = reshape(XX,ki,ki);
    XX = XX + XX';
    try
        if needv,
            [QQ,DD] = eig(XX);
            DD = diag(DD);
        else
            DD = eig(XX);
        end
    catch
        % If eig() fails to converge, fall back onto svd(). This costs
        % more, so we don't want to use it every time.
        [QQ,DD,VV] = svd(XX);
        DD = diag(DD).*sign(real(sum(conj(QQ).*VV)'));
    end
    lab(ei+1:ei+ki) = 0.5*DD;
    ei = ei + ki;
    if needv,
        q(vi+1:vi+qi) = real(QQ);
        vi = vi + qi;
        if i > nr,
            q(vi+1:vi+qi) = imag(QQ);
            vi = vi + qi;
        end
    end
end
