function disp( x, prefix, iname )
if nargin < 2, prefix = ''; end
if nargin < 3, iname = ''; end
nm = cvx_subs2str( x.name_ );
nm = nm(2:end);
str = [ prefix, 'CVX dual variable' ];
if ~isequal( nm, iname ),
    str = [ str, ' ', nm ];
end
tp = type( x );
disp( [ str, ': ', tp ] );
if isequal( tp, 'composite' ),
    global cvx___ %#ok
    tv = subsref( cvx___.problems( x.problem_ ).duals, x.name_ );
    strs = cvx_dispvar( tv, '', true );
    for k = 1 : numel( strs ),
        disp( [ prefix, '   ', strs{k} ] );
    end
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
