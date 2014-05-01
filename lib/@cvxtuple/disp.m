function disp( x, prefix )
if nargin < 2,
    prefix = '';
end
disp( [ prefix, 'CVX tuple object: ' ] );
prefix = [ prefix, '   ' ];
strs = cvx_dispvar( x.value_, '', false );
for k = 1 : length( strs ),
    disp( [ prefix, strs{k} ] );
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
