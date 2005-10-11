function ans = cvx_basis( x, ndx )
error( cvx_verify( x ) );
ans = x.basis_;
if nargin == 2,
    if ndx < size( ans, 2 ),
        ans = reshape( ans( :, ndx + 1 ), x.size_ );
    else,
        ans = reshape( sparse( [], [], [], prod( x.size_ ), 1 ), x.size_ );
    end
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
