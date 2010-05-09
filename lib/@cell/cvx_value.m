function y = cvx_value( x )
global cvx___
if cvx___.hcellfun,
    y = cellfun( @cvx_value, x, 'UniformOutput', false );
else
    sx = size( x );
    y = cell( sx );
    for k = 1 : prod( sx ),
        y{k} = cvx_value( x{k} );
    end
end

% Copyright 2010 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
