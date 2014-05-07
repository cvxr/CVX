function zD = cvx_pushcnstr( zb, iseq, need_dual )

global cvx___
if isa( zb, 'cvx' ),
    zb = cvx_basis( zb );
end
mO = cvx___.n_equality;
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
cvx___.n_equality = cvx___.n_equality + mN;
if nargout
    if nargin == 3 && need_dual
        zD = cvx_invert_structure( zR )';
        zD = sparse( mO + 1 : mO + mN, 1 : mN, 1 ) * zD;
    else
        zD = [];
    end
end
y = any( zb, 2 );
pstr = cvx___.problems( end );
v  = pstr.t_variable;
nv = numel( v );
if nv <= 1, return; end
ny = numel( y );
if ny < nv,
    v( find( y ) ) = 1; %#ok
else
    v( y( 1 : nv ) ) = 1;
end
cvx___.problems( end ).t_variable = v;

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
