function [delta,h,alpha] = iswnbr(w,thetaSQR)
% [delta,h,alpha] = iswnbr(vSQR,thetaSQR)
%
% ISWNBR  Checks feasibility w.r.t. wide region/neighborhood of Sturm-Zhang.
%   vTAR:= (1-alpha)*max(h,v) projection v onto theta-central region
%   delta = (sqrt(n)/theta) * norm(vTAR - v) / norm(v)
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
%

%Indicate to the user Matlab cannot find the SeDuMi binaries
sedumi_binary_error();

% ----------------------------------------
% r = n/thetaSQR
% hSQR = sumwNT/(r-|T|),  hubSQR = sumwNT/(r-|T| - |Q|)
% sumdifv = h*|T| - sumvT   (sumvT = sum(v_T), growing )
% sumdifw = hSQR*|T| - sumwT
% alpha = sumdifv / (r*h)
% deltaSQR = r * ( 2*alpha-alpha^2 - (1-alpha)^2 * sumdifw/gap )
% WE UPDATE sumdifv AND sumdifw IN A STABLE WAY
% ----------------------------------------
n = length(w); gap = sum(w);
sumwNT = gap;
r = n / thetaSQR;
cardT = 0; wQ = []; sumdifv = 0; sumdifw = 0;
cardQ = n;
hSQR = sumwNT / (r - cardT);  hubSQR = sumwNT / (r-(n-1));
for j = 1:n
    wj = w(j);
    if wj >= hubSQR              % wj >= hubSQR ==> not in T
        cardQ = cardQ - 1;
        hubSQR = sumwNT / (r-cardT-cardQ);
    elseif wj < hSQR             % wj < hSQR ==> in T
        cardT = cardT + 1;
        cardQ = cardQ - 1;
        hubSQR = (1-wj/sumwNT) * hubSQR;
        sumwNT = sumwNT - wj;
        oldhSQR = hSQR;
        hSQR = sumwNT / (r - cardT);
        sumdifw = sumdifw + (oldhSQR-wj) + cardT * (hSQR-oldhSQR);
        sumdifv = sumdifv + (sqrt(oldhSQR)-sqrt(wj)) + ...
            cardT * (sqrt(hSQR)-sqrt(oldhSQR));
    else                    % Inconclusive: j in Q
        wQ = [wQ;wj]; %#ok
    end % if
end % for
% ----------------------------------------
% The same treatment for the Q set, but we
% sort the (presumably short) wQ first.
% ----------------------------------------
if ~isempty(wQ)
    sort(wQ);
    STOP = 0; j = 1;
    while ~STOP
        wj = wQ(j);
        if wj >= hSQR
            STOP = 1;
        else
            cardT = cardT + 1;
            sumwNT = sumwNT - wj;
            oldhSQR = hSQR;
            hSQR = sumwNT / (r - cardT);
            sumdifw = sumdifw + (oldhSQR-wj) + cardT * (hSQR-oldhSQR);
            sumdifv = sumdifv + (sqrt(oldhSQR)-sqrt(wj)) + ...
                cardT * (sqrt(hSQR)-sqrt(oldhSQR));
            j = j+1;
            if j > length(wQ)
                STOP = 1;
            end
        end
    end
end % treatment Q
% ----------------------------------------
% alpha = sumdifv/(r*h)
% deltaSQR = r * ( 2*alpha-alpha^2 - (1-alpha)^2 * sumdifw/gap )
%  (THE ABOVE DIFFERENCE SHOULD NOT BE NUMERICALLY DANGEROUS,
%    SINCE alpha IS *SIGNIF* BIGGER THAN sumdifw/gap )
% ----------------------------------------
h = sqrt(hSQR);
alpha = sumdifv/ (r*h);
deltaSQR = alpha*(2-alpha) - (1-alpha)^2 * sumdifw/gap;
delta = sqrt(r*deltaSQR);