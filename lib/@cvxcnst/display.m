function display( x )
long = ~isequal(get(0,'FormatSpacing'),'compact');
if long, disp( ' ' ); end
disp(x);
if long, disp( ' ' ); end

% Copyright 2012 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
