function [At,b,c,K,prep,origcoeff] = pretransfo(At,b,c,K,pars)

% [At,b,c,K,prep] = pretransfo(At,b,c,K)
%
% PRETRANSFO  Checks data and then transforms into internal SeDuMi format.
%
% **********  INTERNAL FUNCTION OF SEDUMI **********
%
% See also sedumi

% Nearly complete rewrite for SeDumI 1.4
% Copyright (c) 2013 Michael C. Grant
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

% -------------------------------------------------------------------------
% Make sure that all fields exist in K, and verify that they are valid
% -------------------------------------------------------------------------

if ~isfield(K,'f') || isempty(K.f)
    K.f = 0;
elseif numel(K.f) ~= 1 || K.f ~= floor(K.f) || K.f < 0 || ~isreal(K.f)
    error('K.f should be nonnegative integer')
end
if ~isfield(K,'l') || isempty(K.l)
    K.l = 0;
elseif K.l ~= floor(K.l) || K.l < 0 || ~isreal(K.l)
    error('K.l should be nonnegative integer')
end
if ~isfield(K,'q') || ~nnz(K.q)
    K.q = zeros(1,0);
else
    K.q = K.q(:)';
    if any(K.q ~= floor(K.q)) || any(K.q<2) || ~isreal(K.q)
        error('K.q should contain only integers bigger than 1')
    end
end
if ~isfield(K,'r') || ~nnz(K.r)
    K.r = zeros(1,0);
else
    K.r = K.r(:)';
    if any(K.r ~= floor(K.r)) || any(K.r<3) || ~isreal(K.r)
        error('K.r should contain only integers bigger than 2')
    end
end
if ~isfield(K,'s') || ~nnz(K.s)
    K.s = zeros(1,0);
else
    K.s = K.s(:)';
    if any(K.s ~= floor(K.s)) || any(K.s<1) || ~isreal(K.s)
        error('K.s should contain only positive integers')
    end
end
% As an alternative to the 'scomplex' flag, I've added a 'z' parameter
% containing a list of Hermitian semidefinite cone sizes. For now, this
% is just translated to 'scomplex' for you. In the future, we may use
% this internally *instead* of scomplex or rsdpN.
if ~isfield(K,'z') || ~nnz(K.z)
    K.z = zeros(1,0);
else
    K.z = K.z(:)';
    if any(K.z ~= floor(K.z)) || any(K.z<1) || ~isreal(K.z)
        error('K.z should contain only positive integers')
    end
end

N_f    = K.f;
N_l    = K.l;
N_fl   = N_f + N_l;
L_q    = length(K.q);
N_q    = sum(K.q);
L_r    = length(K.r);
N_r    = sum(K.r);
N_qr   = N_q + N_r;
L_qr   = L_q + L_r;
L_s    = length(K.s);
L_z    = length(K.z);
L_sz   = L_s + L_z;
N_s    = sum((K.s).^2);
N_z    = sum((K.z).^2);
N_sz   = N_s + N_z;
L_qrsz = L_qr + L_sz;
N_flqr = N_fl + N_qr;
N      = N_flqr + N_sz;

if ~isfield(K,'ycomplex') || isempty(K.ycomplex)
    K.ycomplex = zeros(1,0);
else
    K.ycomplex = sort(K.ycomplex(:))';
    K.ycomplex(find(~diff(K.ycomplex))+1) = [];
    if any(K.ycomplex ~= floor(K.ycomplex)) || any(K.ycomplex<1)
        error('K.ycomplex should contain only positive integers');
    elseif any(K.ycomplex > numel(b))
        error('Elements of K.ycomplex are out of range');
    end
end
if ~isfield(K,'xcomplex') || isempty(K.xcomplex)
    K.xcomplex = zeros(1,0);
else
    K.xcomplex = sort(K.xcomplex(:))';
    K.xcomplex(find(~diff(K.xcomplex))+1) = [];
    if any(K.xcomplex ~= floor(K.xcomplex)) || any(K.xcomplex<1)
        error('K.xcomplex should contain only positive integers');
    elseif any(K.xcomplex>N_flqr)
        error('Elements of K.xcomplex are out of range');
    end
end
if ~isfield(K,'scomplex') || isempty(K.scomplex)
    K.scomplex = zeros(1,0);
else
    K.scomplex = sort(K.scomplex(:))';
    K.scomplex(find(~diff(K.scomplex))+1) = [];
    if any(K.scomplex~=floor(K.scomplex)) || any(K.scomplex<1)
        error('K.scomplex should contain only positive integers');
    elseif any(K.scomplex>L_s)
        error('Elements of K.xcomplex are out of range');
    end
end
if L_z,
    K.s = [ K.s, K.z ];
    K.scomplex = [ K.scomplex, L_s + 1 : L_sz ];
    K.z = zeros(1,0);
    L_s = L_s + L_z;
    N_s = N_s + N_z;
    L_z = 0; %#ok
    N_z = 0; %#ok
end

% -------------------------------------------------------------------------
% Verify the size and validity of At, b, and c
% N = # variables
% -------------------------------------------------------------------------

% SeDuMi assumes that if At is not consistent with K but At' is, that the
% user supplied the transpose of the coefficient matrix. This introduces a
% rare ambiguity in the case where At happens to be square. Past versions 
% of SeDuMi would not have allowed this, instead rejecting the case where
% m >= N; however, with complex variables, the situation is not quite as
% clear, and there may technically be cases where m >= N is acceptable.

if ndims(At) > 2 %#ok
    error('A must be a matrix');
elseif nnz(isnan(At)) || nnz(isinf(At)),
    error('A contains NaN or Inf');
elseif size(At,1) == N,
    % nothing
elseif size(At,2) == N,
    At = At';
else
    error('(At,K) size mismatch');
end
if all(size(b)>1)
    error('Parameter b must be a vector');
elseif any(isnan(b)) || any(isinf(b)),
    error('b contains NaN or Inf');
elseif length(b) ~= size(At,2),
    error('(At,b) size mismatch');
else
    b = b(:);
end
if all(size(c)>1)
    error('Parameter c must be a vector');
elseif any(isnan(c)) || any(isinf(c))
    error('c contains NaN or Inf');
elseif length(c) ~= N
    error('(c,K) size mismatch');
else
    c = c(:);
end

% -------------------------------------------------------------------------
% Save the standardized data for further use if needed
% -------------------------------------------------------------------------

if isfield(pars,'errors') && pars.errors==1
    origcoeff.At=At;
    origcoeff.c=c;
    origcoeff.b=b;
    origcoeff.K=K;
else
    origcoeff=[];
end

% -------------------------------------------------------------------------
% Flag diagonal SDP blocks for removal
% -------------------------------------------------------------------------

% There is some serious MATLAB trickery here (if I do say so myself) that
% merits explanation. "spattern" contains the indices of symbolically 
% nonzero elements of the dual variable z = c - At * y. This is a single
% vector across all LMI constraints, so "sblk" tells us which indices 
% belong to which constraints. The "rem" statement is what determines if
% a particular element is off-diagonal. If there is even one off-diagonal
% element in an SDP, then its corresponding element of "sdiag" is false.
%
% This new version of the analysis replaces the old preprocessSDP(), and
% seems inexpensive enough to apply to all LMIs regardless of size/count.
% The previous version was applied more sparingly.
%
% In theory one could do a more complex analysis of the block structure of
% an LMI, potentially breaking a larger into smaller ones. But this would
% likely be significantly more expensive.

if L_s && ( ~isfield(pars,'sdp') || pars.sdp ),
    ssiz        = (K.s).^2;
    strt        = cumsum([1,ssiz(1:end-1)]);
    sblk        = zeros(1,N_s);
    sblk(strt)  = 1;
    sblk        = cumsum(sblk);
    spattern    = find(c(N_flqr+1:N)~=0|any(At(N_flqr+1:N,:),2))';
    sblk        = sblk(spattern);
    sblk        = sblk(rem(spattern-strt(sblk),K.s(sblk)+1)~=0);
    sdiag       = true(1,L_s);
    sdiag(sblk) = false;
else
    % Even if we disable SDP processing, we're going to move 1x1 SDPs into
    % the nonnegative variable block. It just doesn't make sense to deploy
    % all of that SDP machinery for nonnegative variables.
    sdiag = K.s == 1;
end

% -------------------------------------------------------------------------
% Handle K.ycomplex by splitting apart the complex constraints into pairs
% of real constraints
% -------------------------------------------------------------------------

if ~isempty(K.ycomplex)
    b  = [ real(b) ; imag(b(K.ycomplex)) ];
    At = [ At, 1j * At(:,K.ycomplex) ];
else
    b  = real(b);
end

% -------------------------------------------------------------------------
% Find the locations of the the complex data, so we can convert into
% SeDuMi's internal format, which uses only MATLAB's real representation.
% -------------------------------------------------------------------------
% Strictly speaking, nonnegative variables, the first variable in a Lorentz 
% cone, and the first two variables in a rotated Lorentz cone are real. But
% But SeDuMi has allowed them all to be specified as complex; the imaginary
% portions are interpreted as free variables. We have kept that behavior.

% This code replaces whichcpx.c in its entirety. It actually fixes a bug in
% rotated Lorentz cone handling that was probably never exercised.
if isempty(K.xcomplex) && isempty(K.scomplex),
    K.fcplx = zeros(1,0);
    K.qcplx = zeros(1,0);
    K.rcplx = zeros(1,0);
    scplx   = false(1,L_s);
    sreal   = ~sdiag;
    K.rsdpN = L_s;
    N_fc = 0;
    K.cdim = 0;
else
    xc = K.xcomplex;
    tt = xc <= N_fl;
    K.fcplx = xc(tt);
    xc = xc(~tt) -  N_fl;
    tt = xc <= N_q;
    K.qcplx = xc(tt);
    xc = xc(~tt) - N_q;
    tt = xc <= N_r;
    K.rcplx = xc(tt);
    if ~isempty(K.qcplx),
        ndxs = cumsum([1,K.q(1:end-1)]);
        t2 = any(bsxfun(@eq,K.qcplx,ndxs'),1);
        K.fcplx = [ K.fcplx, K.qcplx(t2) + N_fl ];
        K.qcplx(t2) = [];
        t2 = sum(bsxfun(@gt,K.qcplx,ndxs'),1);
        K.q = K.q + full(sparse(1,t2,1,1,L_q));
        K.qcplx = K.qcplx + (1:length(K.qcplx));
        N_q = N_q + length(K.qcplx);
    end
    if ~isempty(K.rcplx),
        ndxs = cumsum([1,K.r(1:end-1)]);
        t2 = any(bsxfun(@eq,K.rcplx,[ndxs,ndxs+1]'),1);
        K.fcplx = [ K.fcplx, K.rcplx(t2) + N_fl + N_q ];
        K.rcplx(t2) = [];
        t2 = sum(bsxfun(@gt,K.rcplx,ndxs'),1);
        K.r = K.r + full(sparse(1,t2,1,1,L_r));
        % This 2*t2 offset is required because QR makes two accesses
        % each of the first two elements of a rotated Lorentz cone.
        K.rcplx = K.rcplx + (1:length(K.rcplx)) + 2 * t2;
        N_r = N_r + length(K.rcplx);
    end
    N_fc = length(K.fcplx);
    N_f  = N_f + N_fc;
    N_fl = N_f + N_l; %#ok
    N_qr = N_q + N_r;
    scplx = false(1,L_s);
    scplx(K.scomplex&~sdiag) = true;
    sreal = ~scplx & ~sdiag;
    K.rsdpN = nnz(sreal);
    K.cdim = length(K.xcomplex) + sum(K.s(scplx).^2);
end

% -------------------------------------------------------------------------
% We have significantly rewritten this section of the code. This section
% constructs a sparse matrix that represents the following transformations:
%   --- Free variables split into differences of nonnegative variables; OR
%   --- Free variables placed in a Lorentz-cone
%   --- Rotated lorentz cones translated to standard Lorentz cones
%   --- Lorentz cones rearranged to trace block + norm-bound blocks
%   --- Conversion of diagonal SDPs to nonnegative variables
%   --- SDP coefficients moved to the lower triangle for increased sparsity
% This code replaces the rotlorenz, qreshape, and vectril MEX files.
% -------------------------------------------------------------------------

newL = 0;
newQ = zeros(1,0);
ii = {}; jj = {}; vv = {};

% Split free variables into the difference of nonnegatives
if ~isfield( pars, 'free' ) || pars.free == 2 && L_qrsz,
    pars.free = 1;
end
if N_f && ~pars.free,
    jt = [ 1 : K.f, K.fcplx ; 1 : K.f, K.fcplx ];
    vt = [ ones(1,K.f), -1j*ones(1,N_fc) ; -ones(1,K.f), 1j*ones(1,N_fc) ];
    ii{end+1} = 1 : 2 * N_f;
    jj{end+1} = jt(:)';
    vv{end+1} = vt(:)';
    newL = 2 * N_f;
    prep.freeL = N_f;
end

% Copy nonnegative variables without change
if K.l,
    ii{end+1} = newL + 1 : newL + K.l;
    jj{end+1} = K.f + 1 : K.f + K.l;
    vv{end+1} = ones(1,K.l);
    newL = newL + K.l;
end

% Convert diagonal SDPs to nonnegative variables
if any(sdiag),
    dsize = K.s(sdiag);
    sdpL = sum(dsize);
    prep.sdiag = dsize;
    jstrt = cumsum([N_flqr+1,K.s(1:end-1).^2]);
    jstrt = jstrt(sdiag);
    istrt = cumsum([1,dsize(1:end-1)]);
    dsize = dsize + 1;
    dblks = cumsum(full(sparse(1,istrt,1,1,sdpL)));
    ii{end+1} = newL + 1 : newL + sdpL;
    jj{end+1} = jstrt(dblks) + dsize(dblks) .* ( ( 1 : sdpL ) - istrt(dblks) );
    vv{end+1} = ones(1,sdpL);
    newL = newL + sdpL;
end

% Stuff free variables into a Lorentz cone
tr_off = newL;
nb_off = newL + L_qr;
if N_f && pars.free,
    tr_off = tr_off + 1;
    nb_off = nb_off + 1;
    ii{end+1} = nb_off + 1 : nb_off + N_f;
    jj{end+1} = [ 1 : K.f, K.fcplx ];
    vv{end+1} = [ ones(1,K.f), -1j*ones(1,N_fc) ];
    nb_off = nb_off + K.f;
    newQ = N_f + 1;
end

% Rearrange Lorentz cones to trace block + norm-bound blocks
if N_q,
    ndxs      = cumsum([1,K.q(1:end-1)]);
    it        = zeros(1,N_q);
    it(ndxs)  = tr_off + 1 : tr_off + L_q;
    it(it==0) = nb_off + 1 : nb_off + ( N_q - L_q );
    jt = K.f + K.l + 1 : K.f + K.l + N_q;
    vt = ones(1,N_q);
    if ~isempty(K.qcplx),
        jt = jt - cumsum(full(sparse(1,K.qcplx,1,1,N_q)));
        vt(K.qcplx) = -1j;
    end
    ii{end+1} = it(:)';
    jj{end+1} = jt;
    vv{end+1} = vt;
    tr_off    = tr_off + L_q;
    nb_off    = nb_off + N_q - L_q;
end

% Transform rotated Lorentz cones to standard Lorentz cones, and rearrange
% to trace block + norm-bound blocks.
if N_r,
    ndxr       = cumsum([1,K.r(1:end-1)]);
    ndxp       = ndxr + 2*(0:L_r-1); 
    it         = zeros(1,N_r+2*L_r);
    it(ndxp)   = tr_off + 1 : tr_off + L_r;
    it(ndxp+1) = -1;
    it(ndxp+2) = it(ndxp);
    it(it==0)  = nb_off + 1 : nb_off + ( N_r - L_r );
    it(ndxp+1) = it(ndxp+3);
    jt = [ K.f + K.l + N_q + 1, ones(1,N_r+2*L_r-1) ];
    jt([ndxp+1,ndxp+3]) = 0;
    vt = ones(1,N_r+2*L_r);
    vt([ndxp,ndxp+1,ndxp+2]) = sqrt(0.5);
    vt(ndxp+3) = -sqrt(0.5);
    if ~isempty(K.rcplx),
        jt(K.rcplx) = 0;
        vt(K.rcplx) = -1j;
    end
    ii{end+1} = it;
    jj{end+1} = cumsum(jt);
    vv{end+1} = vt;
    nb_off = nb_off + N_r - L_r;
end

% Replace non-diagonal real SDP coefficients with tril(X) + tril(X',-1).
% This cuts the number of nonzeros approximately in half.
if K.rsdpN,
    dsize = K.s(sreal);
    sdpL  = sum(dsize.^2);
    jstrt = cumsum([N_flqr+1,K.s(1:end-1).^2]);
    jstrt = jstrt(sreal);
    istrt = cumsum([1,dsize(1:end-1).^2]);
    dblks = cumsum(full(sparse(1,istrt,1,1,sdpL)));
    istrt = istrt + nb_off;
    dsize = dsize(dblks);
    istrt = istrt(dblks);
    jndxs = ( nb_off + 1 : nb_off + sdpL ) - istrt;
    cols  = floor(jndxs ./ dsize);
    rows  = jndxs - dsize .* cols;
    ii{end+1} = max(rows,cols) + min(rows,cols) .* dsize + istrt;
    jj{end+1} = jndxs + jstrt(dblks);
    vv{end+1} = ones(1,sdpL);
    nb_off = nb_off + sdpL;
    clear dsize jstrt istrt dblks jndxs rows cols
end

% Replace Hermitian SDP coefficients with tril(X) + tril(X',-1). This one's
% a bit trickier because we have the real and complex values interleaved, 
% and the imaginary values along the diagonal are zero.
if K.rsdpN < length(K.s),
    dsize = K.s(scplx);
    jsize = dsize .^ 2;
    sdpL  = 2 * sum(jsize);
    jstrt = cumsum([N_flqr+1,K.s(1:end-1).^2]);
    jstrt = jstrt(scplx);
    bstrt = cumsum([1,2*jsize(1:end-1)]);
    dblks = cumsum(full(sparse(1,bstrt,1,1,sdpL)));
    istrt = bstrt + nb_off;
    dsize = dsize(dblks);
    istrt = istrt(dblks);
    bndxs = ( nb_off + 1 : nb_off + sdpL ) - istrt;
    cols  = floor( bndxs ./ dsize );
    rows  = bndxs - dsize .* cols;
    imgv  = cols >= dsize;
    cols  = cols - imgv .* dsize;
    indxs = max(rows,cols) + min(rows,cols) .* dsize + imgv .* jsize(dblks) + istrt;
    vals  = ( 1 - 2 * ( cols > rows ) ) .* ( 1 - ( 1 + 1j ) .* imgv );
    keep  = ~imgv | ( rows ~= cols );
    jndxs = rows + cols .* dsize + jstrt(dblks);
    ii{end+1} = indxs(keep);
    jj{end+1} = jndxs(keep);
    vv{end+1} = vals(keep);
    clear dsize jsize jstrt bstrt istrt bndxs rows cols vals imgv keep
end

% Update free, nonnegative, and Lorentz variable counts
K.f = 0;
K.l = newL;
K.q = [ newQ, K.q, K.r ];
K.r = zeros(0,1);
K.s = [ K.s(:,~scplx&~sdiag), K.s(scplx&~sdiag) ];
K.rsdpN = nnz(~scplx&~sdiag);
K.N = K.l + sum(K.q) + sum(K.s(1:K.rsdpN).^2)+2*sum(K.s(K.rsdpN+1:end).^2);

% Create the artificial (x0,z0) variable for the self-dual model by
% appending a zero row to At and c. This is accomplished in QR by adding
% 1 to all of the row indices created above.
K.N  = K.N + 1;
K.l  = K.l + 1;
K.m  = length(b);

% Transform At, c
% The transformation matrix QR does not satisfy QR'*QR=I. But it does, in 
% fact, serve as a reverse transformation:
%    --- For Lorentz cones, QR applies a permutation; QR' reverses it.
%    --- For rotated Lorentz cones, QR applies a unitary, self-adjoint
%        rotation on the first two variables, so QR' reverses it.
%    --- For split free variables, QR creates the positive and negative
%        parts; QR' combines them back together.
%    --- For free variables placed in a Lorentz cone, QR moves them to the
%        Lorentz cone block, adding an extra epigraph variable. QR' moves
%        them back and drops the extra variable.
%    --- For semidefinte cones, QR adds the strict upper triangle to the
%        strict lower triangle; QR' copies the strict lower triangle to the
%        strict upper triangle, ensuring symmetry.
%    --- QR adds a row for the self-dual variable; QR' removes it.
[dummy,ndxs] = sort(cellfun(@(x) x(1),jj)); %#ok
QR = sparse( horzcat(ii{ndxs})+1, horzcat(jj{ndxs}), horzcat(vv{ndxs}), K.N, length(c) );
At = real( sparse( QR * At ) );
c  = real( sparse( QR * c  ) );
b  = sparse( b );
prep.QR = QR;
clear ii jj vv

% -------------------------------------------------------------------------
% Now K has field K.{l,q,s}
% Generate a more detailed description of cone K:
% Let K.blkstart(i):K.blkstart(i+1)-1 give the index range of block i.
% Compute maxN=max(order), len=sum(order) for LORENTZ, real PSD, herm PSD
% yields: K.{l,q,s,rsdpN,blkstart,rLen,hLen,qMaxn,rMaxn,hMaxn}
% -------------------------------------------------------------------------
Ksr = K.s(1:K.rsdpN);
Ksc = K.s(K.rsdpN+1:end);
K.blkstart = cumsum([K.l+1,length(K.q)+length(K.r),K.q-1,Ksr.^2,2*Ksc.^2]);
K.rLen = sum(Ksr);
K.hLen = sum(Ksc);
K.qMaxn = max([0,K.q]);
K.rMaxn = max([0,Ksr]);
K.hMaxn = max([0,Ksc]);
K.mainblks = K.blkstart(cumsum([1 1 length(K.q)]));
K.qblkstart = K.blkstart(2:2+length(K.q));  % Also include blkend
K.sblkstart = K.blkstart(2+length(K.q):end);
K.lq = K.mainblks(end)-1;
