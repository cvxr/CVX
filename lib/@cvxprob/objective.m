function y = objective( prob )
global cvx___
p = cvx___.problems( index( prob ) ).objective;

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
