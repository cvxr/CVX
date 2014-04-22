function x = reshape( x, varargin )

%   Disciplined convex/geometric programming information for RESHAPE:
%       RESHAPE imposes no convexity restrictions on its arguments.

switch nargin,
case {0,1},
    error( 'CVX:ArgError', 'Not enough input arguments.' );
case 2,
    fnan = 0;
    sz = varargin{1};
    if ~( isnumeric(sz) && numel(sz) == length(sz) && all(sz>=0) && all(sz==floor(sz)) )
        error( 'CVX:ArgError', 'Size argument must be a vector of nonnegative integers.' );
    end
otherwise,
    nel = nargin - 1;
    sz = zeros( 1, nel );
    fnan = 0;
    for k = 1 : nel,
        s = varargin{k};
        if isempty( s ),
            if ~fnan,
                error( 'CVX:ArgError', 'Size can only have one unknown dimension.' );
            end
            sz(k) = 1;
            fnan = k;
        elseif isnumeric( s ) && length( s ) == 1,
            sz(k) = s;
        else
            error( 'CVX:ArgError', 'Size arguments must be nonnegative integers.' );
        end
    end
end

sx = x.size_;
px = prod(sx);
if fnan,
    s = px / prod(sz);
    if isnan(s) || s ~= floor(s),
        error( 'Product of known dimensions, %d, not divisible into total number of elements, %d.', prod(sz), prod(sx) );
    end
    sz(fnan) = s;
end
if isequal( sx, sz ),
    return;
elseif px ~= prod( sz ),
    error( 'To RESHAPE the number of elements must not change.' );
end

x = cvx( sz, x.basis_ );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
