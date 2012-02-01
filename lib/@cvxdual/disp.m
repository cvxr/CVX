function disp( x, prefix, iname )
if nargin < 2, prefix = ''; end
if nargin < 3, iname = ''; end
nm = cvx_subs2str( x.name_ );
nm = nm(2:end);
if ~isequal( nm, iname ),
    disp( [ prefix, 'cvx dual variable ', nm, ' (', type( x ), ')' ] );
else
    disp( [ prefix, 'cvx dual variable (', type( x ), ')' ] );
end

% Copyright 2012 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
