function y = permute( x, order )

%   Disciplined convex/geometric programming information for PERMUTE:
%      PERMUTE(A,ORDER) imposes no convexity restrictions on A. ORDER
%      must be constant.

%
% Determine the permutation
%

s = x.size_;
ndxs = 1 : prod( s );
try
    ndx2 = permute( reshape( ndxs, s ), order );
catch
    error( lasterror );
end

%
% Permute the data
%

b = x.basis_;
try
    b = x.basis_( :, ndx2 );
catch
    ndxs( ndxs2( : ).' ) = ndxs;
    [ r, c, v ] = find( b );
    b = sparse( r, ndxs( c ), v, size( b, 1 ), size( b, 2 ) );
    clear r c v
end
y = cvx( size( ndx2 ), b );

% Copyright 2008 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
