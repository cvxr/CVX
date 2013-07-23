function z = qjmul(x,y,K)
% z = qjmul(x,y,K)
%
% QJMUL  Implements Jordan product for Lorentz cones
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

if isempty(K.q)
    z = zeros(0,1);
else
    % ------------------------------------------------------------
    % Let i1, i2 such that x(i1:i2-1) = "x1", i.e. Lorentz trace part.
    % ------------------------------------------------------------
    if length(x) ~= length(y)
        error('x,y size mismatch');
    end
    ix = K.mainblks;
    if length(x) < K.lq        % Only q-part given?
        ix = (1-ix(1)) + ix;
        if length(x) ~= length(y)
            error('x, y size mismatch')
        end
    end
    % ------------------------------------------------------------
    % Let z1 = x'*y/sqrt2, z2 = (x1*y2+y1*x2)/sqrt2
    % ------------------------------------------------------------
    z1 = x(ix(1):ix(2)-1).*y(ix(1):ix(2)-1)...
        + ddot(x(ix(2):ix(3)-1),y,K.qblkstart);
    z = [z1; qblkmul(x(ix(1):ix(2)-1),y,K.qblkstart)...
        + qblkmul(y(ix(1):ix(2)-1),x,K.qblkstart)] / sqrt(2);
end