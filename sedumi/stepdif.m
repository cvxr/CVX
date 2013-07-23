function [t,rcdx] = stepdif(d,R,y0,x,y,z,dy0,dx,dy,dz,b,mint,tpmtd)
% [t,rcdx] = stepdif(d,R,y0,x,y,z,dy0,dx,dy,dz,b,mint,tpmtd)
%
% STEPDIF  Implements Primal-Dual Step-Differentiation for self-dual model
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
% Self-Dual Primal/Dual Step length differentiation:
% ------------------------------------------------------------
d0 = sqrt(d.l(1));
% ------------------------------------------------------------
% Compute (rdy0, rcdx) := ((dy0,cdx) - rdx0 * (y0,cx))/y0 using that
% cx - by + z0 = y0*R.sd.
% ------------------------------------------------------------
rdx0 = dx(1) / x(1);
rdy0 = dy0 / y0 - rdx0;
rcdx = full(b'*dy - rdx0*(b'*y)) - (dz(1)-rdx0*z(1))/d0;
rcdx = rdy0 * R.sd + rcdx / y0;
gap = R.b0 * y0;
if tpmtd > 0
    del1 = (z'*dx) / gap;
    dRg = rdx0 * R.sd + rcdx;
else
    del1 = (x'*dz) / gap;
    dRg = (dy0/y0) * R.sd - rcdx;
end
% ------------------------------------------------------------
% CASE Rsd > 0: gap constraint locally part of merit.
% CASE Rsd < 0: gap constraint locally not part of merit.
% ------------------------------------------------------------
usegap = (R.sd > 0) | (R.sd == 0 & dRg > 0);
if usegap
    r0 = R.w(1)+R.w(2)+R.sd;
    beta = (rdy0 * R.w(1) + rcdx)/ r0;
else
    r0 = R.w(1)+R.w(2);
    beta = rdy0 * R.w(1)/ r0;
end
% ------------------------------------------------------------
% CASE tp > td: merit is multiplied by 1 + t * (rdx0 + beta)
% CASE tp < td: merit is multiplied by 1 + t * (dy0/y0 - beta)
% ------------------------------------------------------------
if tpmtd > 0
    beta = rdx0 + beta;
else
    beta = (dy0/y0) - beta;
end
% ------------------------------------------------------------
% Compute the minimizer t of the local merit function
% ------------------------------------------------------------
c = 2*[beta;rdx0*del1] - (rdx0 + del1)*[1;beta];
if c(1) <= 0
    % ------------------------------------------------------------
    % If merit decreasing at t=0, i.e. c(1)<0, then consider t > 0
    % ------------------------------------------------------------
    if c(2) >= 0
        t = abs(tpmtd);        % merit decreasing for all feasible t>=0
    else
        t = min(abs(tpmtd),c(1) / c(2));  % Set derivative c(1)-t*c(2) to zero
    end
else
    % ------------------------------------------------------------
    % If merit increasing at t=0, i.e. c(1)>0, then consider t < 0
    % ------------------------------------------------------------
    if c(2) >= 0
        t = mint;        % c(1)-t*c(2) > 0 for all feasible t<=0
    else
        t = max(mint,c(1) / c(2));  % Set derivative c(1)-t*c(2) to zero
    end
end
% ------------------------------------------------------------
% Compute the break point, where R.sd + t*dRg = 0.
% ------------------------------------------------------------
if dRg ~= 0
    tg = -R.sd / dRg;
else
    tg = t;    % no break point, keep optimal t
end
if tg <= 0
    if t > 0
        tg = t;    % no positive break point, keep optimal t > 0.
    end
else
    if t < 0
        tg = t;    % break point positive, keep optimal t < 0.
    end
end
% ------------------------------------------------------------
% If t beyond break point, then minimize merit for t >= tg.
% ------------------------------------------------------------
if abs(t) > abs(tg)
    if usegap             % Now, EXCLUDE gap for |t|>|tg|:
        beta = rdy0 * R.w(1)/ r0;
        alpha = 1 - R.sd / r0;   % Remove contribution from Rsd
    else                  % Now, INCLUDE gap for |t|>|tg|:
        beta = (rdy0 * R.w(1) + rcdx)/ r0;
        alpha = 1 + R.sd / r0;   % Include contribution from Rsd
    end
    if tpmtd > 0
        beta = rdx0*alpha + beta;
    else
        beta = (dy0/y0)*alpha - beta;
    end
    % ------------------------------------------------------------
    % Compute the minimizer |t| >= |tg| of the merit function beyond tg
    % ------------------------------------------------------------
    c = 2*[beta;rdx0*del1] - (rdx0 + del1)*[1;beta];
    % ------------------------------------------------------------
    % If t >= tg >= 0:
    % ------------------------------------------------------------
    if t >= 0
        if c'*[1;-tg] <= 0
            if c(2) >= 0
                t = abs(tpmtd);        % c(1)-t*c(2)<=0 for all feasible t>=tg
            else
                t = min(abs(tpmtd),c(1) / c(2));  % Set c(1)-t*c(2) to zero
            end
        else
            t = tg;     % merit is increasing at t=tg >= 0.
        end
        % ------------------------------------------------------------
        % If t <= tg <= 0:
        % ------------------------------------------------------------
    else
        if c'*[1;-tg] >= 0
            if c(2) >= 0
                t = mint;        % c(1)-t*c(2)>=0 for all feasible t<=tg
            else
                t = max(mint,c(1) / c(2));  % Set c(1)-t*c(2) to zero
            end
        else
            t = tg;     % merit is decreasing at t=tg <= 0.
        end
    end
end % |t|>|tg|
% ------------------------------------------------------------
% Observe the break-point of Rp (tpmtd>0) or Rd (tpmtd<0), where
% y0+t*dy0 = 0. Hitting this makes (P) resp. (D) feasible.
% ------------------------------------------------------------
if y0+t*dy0 <= 0
    t = -y0/dy0;        % For simplicity, we don't search beyond break point.
end
rcdx = y0*rcdx;
