function x = reshape( x, varargin )

%   Disciplined convex/geometric programming information for RESHAPE:
%       RESHAPE imposes no convexity restrictions on its arguments.

sx = x.size_;
switch nargin,
case 2,
    sz = varargin{1};
    if isequal( sx, sz ), return; end
    fnan = 0;
    nel = numel(sz);
    if ~( isnumeric(sz) && size(sz,2) == nel && all(sz>=0) && all(sz==floor(sz)) )
        cvx_throw( 'Size argument must be a row vector with integer elements.' );
    elseif nel < 2,
        cvx_throw( 'Size vector must have at least two elements.' );
    end
case {0,1},
    cvx_throw( 'Not enough input arguments.' );
otherwise,
    nel = nargin - 1;
    sz = zeros( 1, nel );
    fnan = 0;
    for k = 1 : nel,
        s = varargin{k};
        if isempty( s ),
            if ~fnan,
                cvx_throw( 'Size can only have one unknown dimension.' );
            end
            sz(k) = 1;
            fnan = k;
        elseif isnumeric( s ) && length( s ) == 1,
            sz(k) = s;
        else
            cvx_throw( 'Size arguments must be nonnegative integers.' );
        end
    end
    if isequal( sx, sz ), return; end
end
px = prod(sx);
pz = prod(sz);
if fnan
    if rem( px, pz ) ~= 0,
        cvx_throw( 'Product of known dimensions, %d, not divisible into total number of elements, %d.', prod(sz), prod(sx) );
    end
    sz(fnan) = px / pz;
elseif px ~= pz
    cvx_throw( 'To RESHAPE the number of elements must not change.' );
elseif nel > 2 && sz(nel) == 1
    nel = max([2,find(sz~=1,1,'last')]);
    sz = sz(1:nel);
end
x = cvx( sz, x.basis_ );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
