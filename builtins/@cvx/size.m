function varargout = size( x, dim )

%   Disciplined convex/geometric programming information for SIZE:
%       SIZE imposes no convexity restrictions on its first argument.
%       The second argument, if supplied, must be a positive integer.

if nargin == 1,
    dim = [];
end
s = x.size_;
if ~isempty( dim ),
    if nargout > 1,
        error( 'Too many output arguments.' );
    elseif ~isnumeric( dim ) | length( dim ) ~= 1 | dim <= 0 | dim ~= floor( dim ),
        error( 'Dimension argument must be a positive integer scalar.' );
    elseif dim > length( s ),
        s = 1;
    else
        s = s( dim );
    end
end
if nargout <= 1,
    varargout{1} = s;
else
    if nargout > length( s ),
        s( nargout ) = 1;
    elseif nargout < length( s ),
        s( nargout ) = prod( s( nargout + 1 : end ) );
    end
    for k = 1 : nargout,
        varargout{k} = s(k);
    end
end

% Copyright 2007 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
