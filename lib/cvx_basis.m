function y = cvx_basis( x )

if isempty( x )
    y = sparse( 1, 0 );
else
    y = sparse( reshape( x, 1, numel(  x  ) ) );
end

% Copyright 2010 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
