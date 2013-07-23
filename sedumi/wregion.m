function [xscl,y,zscl,y0, w,relt, dxmdz,err, wr] = ...
    wregion(L,Lden,Lsd,d,v,vfrm,A,DAt,dense, R,K,y,y0,b, pars, wr)
%  [xscl,y,zscl,y0, w,relt, dxmdz,err, wr] = wregion(L,Lden,Lsd,...
%                             d,v,vfrm,A,DAt,dense, R,K,y,y0,b, pars, wr)
%
% WREGION  Implements Sturm-Zhang Wide-region Interior Point Method.
%
% ******************** INTERNAL FUNCTION OF SEDUMI ********************
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
%

STOP = 0;
n = length(vfrm.lab);
% ----------------------------------------
% INITIAL CENTERING
% ----------------------------------------
if wr.delta > 0.0           % initial centering needed ?
    vTAR = (1-wr.alpha) * max(wr.h,vfrm.lab);
    pv = 2*(vTAR-vfrm.lab);
    pMode = 1;
    [dx, dy, dz, dy0, errc] = sddir(L,Lden,Lsd,pv,...
        d,v,vfrm,A,DAt,dense, R,K,y,y0,b, pars,pMode);
    % ------------------------------------------------------------
    % Take initial centering step
    % ------------------------------------------------------------
    xc = v + dx;              %scaled
    zc = v + dz;
    yc = y + dy;
    y0c = y0+dy0;
    uxc.tdet = tdet(xc,K);
    uzc.tdet = tdet(zc,K);
    [uxc.u,xispos] = psdfactor(xc,K);
    [uzc.u,zispos] = psdfactor(zc,K);
    critval = max(y0, sqrt(min(d.l(1),1/d.l(1)))*v(1));    % >= y0
    critval = max(1E-3, pars.cg.restol) * critval * R.maxRb;
    if (~xispos) || (~zispos) || (errc.maxb > critval) ...
            || ( (~isempty(uxc.tdet)) && (min(uxc.tdet) <= 0.0)) ...
            || ( (~isempty(uzc.tdet)) && (min(uzc.tdet) <= 0.0))
        STOP = -1;   % Reject and terminate
        dxmdz = [];
        err = errc;
    end
    pv = -vTAR;
else
    vTAR = vfrm.lab;
    xc = v;                    % No initial centering needed.
    ix = K.mainblks;
    uxc.tdet = 2*vfrm.lab(ix(1):ix(2)-1) .* vfrm.lab(ix(2):2*ix(2)-ix(1)-1);
    [uxc.u] = psdfactor(xc,K);
    zc = v; uzc = uxc;          %zscl=xscl=v
    yc = y;
    y0c = y0;
    errc.b = sparse(length(y),1);
    errc.maxb = 0; errc.db0 = 0;
    pv = [];              % Means pv = -v, the predictor direction.
    pMode = 2;
end
% ----------------------------------------
% PREDICTOR
% ----------------------------------------
if STOP ~= -1
    [dx,dy,dz,dy0, err] = sddir(L,Lden,Lsd,pv,...
        d,v,vfrm,A,DAt,dense, R,K,y,y0,b, pars,pMode);
    dxmdz = dx-dz;
    if (pars.alg ~= 0)
        pMode = 3;
        gd1 = [dxmdz(1:K.l)./vTAR(1:K.l);qinvjmul(vTAR,vfrm.q,dxmdz, K);...
            psdinvjmul(vTAR,vfrm.s,dxmdz, K)];
        maxt1 = min(maxstep(dx,xc,uxc,K), maxstep(dz,zc,uzc,K));
        % ----------------------------------------
        % 2ND ORDER CORRECTOR
        % alg == 1 :  2nd order expansion of v = sqrt(eig(D(x)z)), like Sturm-Zhang.
        % alg == 2 :  2nd order expansion of vsqr = eig(D(x)z), like Mehrotra.
        % ----------------------------------------
        switch pars.alg
            case 1 % alg 1: v-expansion
                tTAR = 1-(1-maxt1);
                pv = tTAR^2*[gd1(1:K.l).*dxmdz(1:K.l); qjmul(gd1,dxmdz,K);...
                    psdjmul(gd1,dxmdz,K)];
                pv2 = 2*tTAR*(1-tTAR)*((sum(vTAR)/n)*ones(n,1) - vTAR) -(2*tTAR)*vTAR;
            case 2 % alg 2: v^2, like Mehrotra.
                tTAR = 1-(1-maxt1)^3;
                pv = (tTAR / 4) * [gd1(1:K.l).*dxmdz(1:K.l); qjmul(gd1,dxmdz,K);...
                    psdjmul(gd1,dxmdz,K)];
                pv2 = ((1-tTAR)*tTAR*R.b0*y0/n)./vTAR - (1+tTAR/4)*vTAR;
        end
        pv = pv + frameit(pv2, vfrm.q, vfrm.s, K);
        [dx,dy,dz,dy0, err] = sddir(L,Lden,Lsd,pv,...
            d,v,vfrm,A,DAt,dense, R,K,y,y0,b, pars,pMode);
    end
    % ----------------------------------------
    % The steplength should be such that
    % y0(t) * maxRb + t*errb <= (1+phi*t*|dy0/y0|) * y0(t)*maxRb
    % i.e. errb <= phi*(-dy0) * (1+t*dy0/y0) * maxRb
    % ----------------------------------------
    PHI = 0.5;
    if dy0 < 0 && (PHI*dy0^2*R.maxRb)~=0
        critval =  - (PHI * dy0*R.maxRb + err.maxb)*y0c / (PHI*dy0^2*R.maxRb);
    else
        critval = 1;
    end
    if critval <= 0
        STOP = -1;
    else
        % ----------------------------------------
        % Compute boundary step length "maxt".
        % ----------------------------------------
        tp = maxstep(dx,xc,uxc,K);
        td = maxstep(dz,zc,uzc,K);
        if dy0 < 0                     % keep (t/y0Next)*errmaxb <= R.maxRb.
            tp = min(tp,critval);
        end
        if xc(1)+ td * dx(1) < 0      % Cannot step beyond x0=0.
            td = xc(1) / (-dx(1));
        end
        maxt = min(tp,td);
        % ----------------------------------------
        % CREGION NBRHD STEP OF PRED+COR
        % ----------------------------------------
        [t,wr,w] = widelen(xc,zc,y0c, dx,dz,dy0,0, maxt,pars,K);
        % ----------------------------------------
        % Take step (in scaled space)
        % ----------------------------------------
        xscl = xc + t*dx;
        y    = yc + t*dy;
        zscl = zc + t*dz;
        y0   = y0c + t*dy0;
        if pars.stepdif == 1
            [tdif,rcdx]=stepdif(d,R,y0,xscl,y,zscl,dy0,dx,dy,dz,b,-t,tp-td);
            if tdif ~= 0
                rdx0 = dx(1) / xscl(1);
                mu = 1+tdif*rdx0;
                if tp > td
                    newx = xscl + tdif*dx;
                    newz = mu*zscl;
                else
                    newx = mu*xscl;
                    newz = zscl+tdif*dz;
                end
                [tdif,wr,w] = trydif(tdif,wr,w, newx,newz,pars,K);
            end
        else
            tdif = 0;
        end
        if tdif ~= 0
            rdy0 = dy0 - rdx0 * y0;
            zscl = newz;
            xscl = newx;
            if tp > td
                y = mu*y;
                y0 = mu*y0;
                err.b = (tdif*rdy0) * R.b + errc.b + (t+tdif) * err.b;
                err.g = tdif * rcdx;
                relt.p = (t+tdif)/tp;
                relt.d = t/td;
            else   % td > tp
                y = y + tdif * dy;
                y0 = y0 + tdif * dy0;
                err.b = -(tdif*rdy0)*R.b + mu*(errc.b + t * err.b);
                err.g = -tdif * rcdx;
                relt.p = t/tp;
                relt.d = (t+tdif)/td;
            end
        else
            err.b = errc.b + t * err.b;
            err.g = 0;
            relt.p = t / maxt;
            relt.d = relt.p;
        end
        wr.tpmtd = tp-td;
        err.maxb = (errc.maxb + t*err.maxb);   % to be divided by y0
        err.db0 = xscl'*zscl - y0*R.b0;        % to be divided by y0
    end % critval > 0
end % STOP ~= -1
if STOP == -1
    relt.p = 0;
    relt.d = 0;
    w = [];
    xscl = [];
    zscl = [];
    err.b=zeros(size(b));
    err.db0 = 0;
    err.g = 0;
end
