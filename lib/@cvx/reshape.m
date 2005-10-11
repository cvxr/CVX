function x = reshape( x, varargin )
error( cvx_verify( x ) );

%
% Check size arguments
%

switch nargin,
    case {0,1},
        error( 'Not enough input arguments.' );
    case 2,
        [ temp, sz ] = cvx_check_dimlist( varargin{1}, true );
        if ~temp,
            error( 'Second argument must be a valid dimension list.' );
        end
    otherwise,
        [ temp, sz ] = cvx_check_dimlist( varargin, true );
        if ~temp,
            error( 'Second and subsequent arguments must be nonnegative integers.' );
        end
end

%
% Quick exit if the size remains the same
%

sx = size( x );
if isequal( sx, sz ),
    return;
end

%
% Confirm compatible reshape
%

if prod( sx ) ~= prod( sz ),
    error( 'To RESHAPE the number of elements must not change.' );
end

%
% Perform the resize
%

x = cvx( problem( x ), sz, cvx_basis( x ) );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
