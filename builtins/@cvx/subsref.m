function x = subsref( x, S )

%   Disciplined convex/geometric programming information for SUBSREF:
%      The use of subscripts to extract elements or "slices" of any CVX
%      variable is identical to their use with numeric arrays. All 
%      conventions are preserved, including the colon ':' and 'end'.

sx = x.size_;
if numel( S ) == 1 && isequal( S.type, '()' ) && all( strcmp( S.subs, ':' ) ),
    ns = numel(S.subs);
    sx( end + 1 : numel(ns) ) = 1;
    sx = [ sx(1:ns-1), prod(sx(ns:end)) ];
    x = cvx( sx, x.basis_ );
else
    try
        nx = prod( sx );
        ndxs = builtin( 'subsref', reshape( 1 : nx, sx ), S );
    catch errmsg
        throw( errmsg );
    end
    x = cvx( size( ndxs ), x.basis_( :, ndxs ) );
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
