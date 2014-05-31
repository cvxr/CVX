function [ s, varargout ] = size( x, dim )

%   Disciplined convex/geometric programming information for SIZE:
%       SIZE imposes no convexity restrictions on its first argument.
%       The second argument, if supplied, must be a positive integer.

s = x.size_;
if nargin > 1,
    if nargout > 1,
        cvx_throw( 'Too many output arguments.' );
    elseif ~isnumeric( dim ) || length( dim ) ~= 1 || dim <= 0 || dim ~= floor( dim ),
        cvx_throw( 'Dimension argument must be a positive integer scalar.' );
    elseif dim > length( s ),
        s = 1;
    else
        s = s(dim);
    end
elseif nargout > 1,
    no = nargout;
    s( end + 1 : no ) = 1;
    s( no ) = prod( s( no : end ) );
    for k = 2 : no,
        varargout{k-1} = s(k); %#ok
    end
    s = s(1);
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
