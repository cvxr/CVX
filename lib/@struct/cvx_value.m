function x = cvx_value( x )
x = cell2struct( cvx_value( struct2cell( x ) ), fieldnames( x ), 1 );

% Copyright 2010 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
