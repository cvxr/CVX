function ans = cvx_basis( x, ndx )
error( cvx_verify( x ) );
ans = x.basis_;
if nargin == 2,
    s = x.size_;
    if ndx < size( ans, 2 ),
        ans = ans( :, ndx + 1 );
    else,
        ans = sparse( prod(s), 1 );
    end
    if length( s ) > 2,
        ans = full( ans );
    elseif issparse( ans ),
        c = 2*~isreal(ans);
        if (3+c)*nnz(ans) >= ((2+c)*s(1)-1)*s(2)-1,
            ans = full( ans );
        end
    end
    ans = reshape( ans, s );
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
