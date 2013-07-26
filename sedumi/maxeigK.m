function lab = maxeigK(x,K)

% lab = maxeigK(x,K)
%
% MAXEIGK  Computes the maximum eigenvalue of a vector x with respect to a 
%          self-dual homogenous cone K.
%
% See also sedumi, mat, vec, eyeK.

% New function by Michael C. Grant
% Copyright (C) 2013 Michael C. Grant.
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

% The existence of rsdpN is code for 'is this the internal format?'
if isfield(K,'rsdpN'),
    % The internal SeDuMi cone format.
    is_int = true;
    nf = 0;
    nl = K.l;
    nq = length(K.q);
    nr = 0;
    ns = length(K.s);
    nrsdp = K.rsdpN;
else
    % The external SeDuMi cone format.
    is_int = false;
    if isfield(K,'f'), nf = K.f; else nf = 0; end
    if isfield(K,'l'), nl = K.l; else nl = 0; end
    if isfield(K,'q'), nq = length(K.q); else nq = 0; end
    if isfield(K,'r'), nr = length(K.r); else nr = 0; end
    if isfield(K,'s'), ns = length(K.s); else ns = 0; end
end
xi = nf;
lab = max([-Inf;x(xi+1:xi+nl)]);
xi = xi + nl;
if nq,
    tmp = sqrt(0.5);
    if is_int,
        % Internal version: all of the x0 values are at the front, and the
        % vectors are stacked in order after that.
        zi = xi;
        xi = xi + nq;
        for i = 1:nq,
            kk = K.q(i) - 1;
            x0 = x(zi+i);
            lab = max(lab,tmp*(x0+norm(x(xi+1:xi+kk))));
            xi = xi + kk;
        end
    else
        for i = 1:length(K.q)
            kk = K.q(i);
            x0 = x(xi+1);
            lab = max(lab,tmp*(x0+norm(x(xi+2:xi+kk))));
            xi = xi + kk;
        end
    end
end
for i = 1:nr,
    % Only the external format need be implemented here
    ki = K.r(i);
    x1 = xx(xi+1);
    x2 = xx(xi+2);
    lab = max(lab,0.5*(x1+x2+norm([x1-x2;2*x(xi+3:xi+ki)])));
    xi = xi + ki;
end
if ns,
    for i = 1 : ns,
        ki = K.s(i);
        qi = ki * ki;
        XX = x(xi+1:xi+qi); 
        xi = xi + qi;
        if i > nrsdp,
            XX = XX + 1i*x(xi+1:xi+qi); 
            xi = xi + qi;
        end
        XX = reshape(XX,ki,ki);
        XX = XX + XX';
        if ki > 500,
            if nnz(XX) < 0.1 * numel(XX), XX = sparse(XX); end
            [v,val,flag] = eigs(XX,1,'LA',struct('issym',true)); %#ok
            if flag, val = max(eig(XX)); end
        else
            val = max(eig(XX));
        end
        lab = max(lab,0.5*val);
    end
end


