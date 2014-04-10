function x = subsref( x, S )

%   Disciplined convex/geometric programming information for SUBSREF:
%      The use of subscripts to extract elements or "slices" of any CVX
%      variable is identical to their use with numeric arrays. All 
%      conventions are preserved, including the colon ':' and 'end'.

sx = x.size_;
nx = prod( sx );
if ~isequal( S.subs, {':'} ),
    try
        ndxs = builtin( 'subsref', reshape( 1 : nx, sx ), S );
    catch errmsg
        throw( errmsg );
    end
    x = cvx( size( ndxs ), x.basis_( :, ndxs ) );
elseif sx(1) ~= nx,
    x = cvx( [ nx, 1 ], x.basis_ );
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
