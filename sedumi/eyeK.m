function x = eyeK(K)

% eyeK    Identity w.r.t. symmetric cone.
%
%    x = eyeK(K) produces the identity solution w.r.t. the symmetric cone,
%    that is described by the structure K. This is the vector for which
%    eigK(x) is the all-1 vector.
%
% See also eigK.

% Complete rewrite for SeDuMi 1.3 by Michael C. Grant
% Copyright (C) 2013 Michael C. Grant.
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

% The existence of rsdpN is code for 'is this the internal format?'
is_int = isfield(K,'rsdpN');
if is_int,
    N = K.N;
else
    N = 0;
    if isfield(K,'f'), N = N + K.f; end
    if isfield(K,'l'), N = N + K.l; end
    if isfield(K,'q'), N = N + sum(K.q); end
    if isfield(K,'r'), N = N + sum(K.r);  end
    if isfield(K,'s'), N = N + sum(K.s.^2); end
    if isfield(K,'z'), N = N + sum(K.z.^2); end
end
x = zeros(N,1);
xi = 0;
if ~is_int && isfield(K,'f'),
    xi = xi + K.f;
end
if isfield(K,'l'),
    x(xi+1:xi+K.l) = 1;
    xi = xi + K.l;
end
if isfield(K,'q') && ~isempty(K.q),
    if is_int,
        % Internal version: all of the x0 values are at the front, and the
        % vectors are stacked in order after that.
        x(xi+1:xi+length(K.q)) = sqrt(2.0);
    else
        tmp = K.q(1:end-1);
        x(K.f+K.k+cumsum([1;tmp(:)])) = sqrt(2.0);
    end
    xi = xi + sum(K.q);
end
if ~is_int && isfield(K,'r') && ~isempty(K.r),
    tmp = K.r(1:end-1);
    tmp = cumsum([1;tmp(:)]);
    x([tmp;tmp+1]) = 1;
    xi = xi + sum(K.r);
end
if isfield(K,'s') && ~isempty(K.s),
    nc = length(K.s);
    if is_int,
        nr = K.rsdpN;
    else
        nr = nc;
    end
    for i = 1 : nc,
        ki = K.s(i);
        qi = ki * ki;
        x(xi+1:ki+1:xi+qi) = 1.0;
        xi = xi + ((1+(i>nr))*qi);
    end
end
if ~is_int && isfield(K,'z') && ~isempty(K.z),
    for i = 1 : length(K.z),
        ki = K.z(i);
        qi = ki * ki;
        x(xi+1:ki+1:xi+qi) = 1.0;
        xi = xi + qi;
    end
end    

