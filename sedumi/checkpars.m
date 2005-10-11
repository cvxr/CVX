%                                          pars = checkpars(pars,lponly)
% CHECKPARS  Fills in defaults for missing fields in "pars" structure.
%
% **********  INTERNAL FUNCTION OF SEDUMI **********
%
% See also sedumi


function pars = checkpars(pars,lponly)
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


% --------------------------------------------------
% Algorithm selection parameters
% --------------------------------------------------
if ~isfield(pars,'alg') | sum([0 1 2] == pars.alg) == 0
    pars.alg = 2;
end
if ~isfield(pars,'beta')       % 0.1 <= beta <= 0.9 (theoretically in (0,1))
    pars.beta = 0.5;
elseif pars.beta > 0.9
    pars.beta = 0.9;
elseif pars.beta < 0.1
    pars.beta = 0.1;
end
if ~isfield(pars,'theta')      % 0.01 <= theta <= 1.0 (theoretically in (0,1])
    if pars.alg == 0
        pars.theta = 1;
    else
        pars.theta = 0.25;
    end
elseif pars.theta > 1.0
    pars.theta = 1.0;
elseif pars.theta < 0.01
    pars.theta = 0.01;
end
if ~isfield(pars,'stepdif')
    pars.stepdif = 2;  %adaptive step-differentiation
end
if ~isfield(pars,'w')
    pars.w = [1;1];
elseif length(pars.w)~=2
    warning('pars.w should be vector of order 2')
    pars.w = [1;1];
else
    pars.w = max(pars.w,1e-8);    % positive weights
end
pars.w = reshape(pars.w,2,1);
% -------------------------------------------------
% Preprocessing
% -------------------------------------------------
% How to treat free variables: 0: split them, 1: put them in Q-cone
% (default)
if ~isfield(pars,'free')
    pars.free=1;
end
if ~isfield(pars,'sdp')
    pars.sdp=1; %Enable/disable SDP preprocessing
end
% --------------------------------------------------
% Initialization
% --------------------------------------------------
if ~isfield(pars,'mu') | pars.mu <= 0
    pars.mu = 1;
end
% --------------------------------------------------
% Stopping and reporting criteria
% --------------------------------------------------
if ~isfield(pars,'fid')                   % iter output file handle
    pars.fid = 1;
end
if ~isfield(pars,'eps')                   % stopping tolerance
    pars.eps = 1E-8;
end
if ~isfield(pars,'bigeps')                   % threshold for numerr=1 vs 2.
    pars.bigeps = 1E-3;
end
if ~isfield(pars,'maxiter')
    pars.maxiter = 150;
end
% --------------------------------------------------
% Debugging and algorithmic/diagnostic analysis
% --------------------------------------------------
if ~isfield(pars,'vplot')
    pars.vplot = 0;
end
if ~isfield(pars,'stopat')
    pars.stopat = -1;
end
if ~isfield(pars,'errors') %Print and return the DIMACS errors
    pars.errors = 1;
end
if ~isfield(pars,'prep') %Print preprocessing information
    pars.prep=1;
end
% --------------------------------------------------
% Dense column handling
% --------------------------------------------------
if ~isfield(pars,'denq')
    pars.denq = 0.75;
end
if ~isfield(pars,'denf')
    pars.denf=10;
end
% --------------------------------------------------
% Numerical control
% --------------------------------------------------
if ~isfield(pars,'numtol')           % Criterion for refinement
    pars.numtol = 5E-7;
end
if ~isfield(pars,'bignumtol')        % Criterion for failure
    pars.bignumtol = 0.9;             % will be multiplied by y0.
end
if ~isfield(pars,'numlvl')
    pars.numlvl = 0;
end
if isfield(pars,'chol')
    cholpars = pars.chol;
    if ~isfield(cholpars,'skip')
        cholpars.skip = 1;
    end
    if ~isfield(cholpars,'abstol')
        cholpars.abstol = 1E-20;
    end
    if ~isfield(cholpars,'canceltol')         % Cancelation in LDL computation
        cholpars.canceltol = 1E-12;
    end
    if ~isfield(cholpars,'maxu')
        cholpars.maxu = 5E5;
    end
    if ~isfield(cholpars,'maxuden')
        cholpars.maxuden = 5E2;
    end
else
    cholpars.skip = 1;
    cholpars.abstol = 1e-20;
    cholpars.canceltol = 1E-12;
    cholpars.maxu = 5E5;
    cholpars.maxuden = 5E2;
end
pars.chol = cholpars;
if isfield(pars,'cg')           % Criterion Conjugate Gradient
    cgpars = pars.cg;
    if ~isfield(cgpars,'qprec')    % Use quadruple precision on main iters
        cgpars.qprec = 1;
    end
    if ~isfield(cgpars,'restol')    % Required relative residual accuracy
        cgpars.restol = 5E-3;
    end
    if ~isfield(cgpars,'stagtol')   % Stagnation in function progress
        cgpars.stagtol = 5E-14;
    end
    if ~isfield(cgpars,'maxiter')   % Max number of cg steps
        cgpars.maxiter = 49;
    end
    if ~isfield(cgpars,'refine')
        cgpars.refine = 1;
    end
else
    cgpars.qprec = 1;
    cgpars.restol = 5E-3;
    cgpars.stagtol = 5E-14;
    cgpars.maxiter = 49;
    cgpars.refine = 1;
end
pars.cg = cgpars;