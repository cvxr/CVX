function zD = cvx_newcnstr( zb, iseq, need_dual )

global cvx___
try
    pstr = cvx___.problems(end);
catch
    error( 'CVX:NoModel', 'No CVX model is present.' );
end
if isa( zb, 'cvx' ),
    zb = cvx_basis( zb );
end
if iseq,
    cmode = 'full';
else
    cmode = 'magnitude';
end
[ zR, zb ] = cvx_bcompress( zb, cmode );
mN = size( zb, 2 );
if ~iseq,
    ndim = cvx_pushnneg( mN );
    zb(ndim(end),1) = 0;
    zb = zb + sparse( ndim, 1:mN, -1 );
end
cvx___.equalities{end+1} = zb;
cvx___.inequality(end+1) = ~iseq;
mO = cvx___.n_equality;
cvx___.n_equality = mO + mN;
if pstr.n_variable > 1 && pstr.complete,
    nv = min(pstr.n_variable,size(zb,1));
    if nnz(zb(2:nv,:)), 
        pstr.complete = false; 
        cvx___.problems(end) = pstr;
    end
end
if nargout
    if nargin == 3 && need_dual
        zD = cvx_invert_structure( zR )';
        zD = sparse( mO + 1 : mO + mN, 1 : mN, 1 ) * zD;
    else
        zD = [];
    end
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
