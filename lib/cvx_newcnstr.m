function zD = cvx_newcnstr( zb, iseq, need_dual )

global cvx___
try
    pstr = cvx___.problems(end);
catch
    cvx_throw( 'No CVX model is present.' );
end
if isa( zb, 'cvx' ),
    zb = cvx_basis( zb );
end
[nN,mN] = size( zb );
isr = isreal( zb );
if ~isr,
    mN = 2 * mN;
    zb = reshape( [ real(zb) ; imag(zb) ], nN, mN );
end
if mN,
    if ~iseq,
        ndim = cvx_pushnneg( mN );
        zb(ndim(end),1) = 0;
        zb = zb + sparse( ndim, 1:mN, -1 );
    end
    cvx___.equalities{end+1} = zb;
    cvx___.inequality(end+1) = ~iseq;
    mO = cvx___.n_equality;
    cvx___.n_equality = mO + mN;
    cp = pstr.checkpoint(4);
    if cp > 1,
        x = find( any( zb, 2 ), 2, 'first' );
        if any( x > 1 & x <= cp ),
            cvx___.problems(end).checkpoint(4) = x(1+(x(1)==1)) - 1;
        end
    end
end
if nargout
    if nargin == 3 && need_dual
        zD = sparse( mO + 1 : mO + mN, 1 : mN, 1 );
        if ~isr, zD = zD(:,1:2:end) + 1j * zD(:,2:2:end); end
    else
        zD = [];
    end
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
