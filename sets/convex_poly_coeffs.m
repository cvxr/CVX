function cvx_optpnt = convex_poly_coeffs( deg, mm )
error( nargchk( 1, 2, nargin ) );

%
% Check argument
%

if ~cvx_check_dimension( deg, true ),
    error( 'Argument must be a nonnegative integer.' );
elseif rem( deg, 2 ) ~= 0 & deg ~= 1,
    error( 'Degree must be 0, 1, or even.' );
end

% Check range argument
%

if nargin < 2 | isempty( mm ),
    mm = [ -Inf, +Inf ];
else
    if ~isa( mm, 'double' ) | ~isreal( mm ) | ndims( mm ) > 2 | numel( mm ) ~= 2 & size( mm, 2 ) ~= 2,
        error( 'Second argument must be a range [ xmin xmax ] or a matrix of them.' );
    end
    mm = reshape( mm, 0.5 * numel( mm ), 2 );
    m1 = mm(:,1);
    m2 = mm(:,2);
    if any( ( m1 == m2 ) & isinf( m1 ) ),
        error( 'Intervals [-Inf,-Inf] and [+Inf,+Inf] are not accepted.' );
    end
end

%
% Construct set
%

cvx_begin_set
    variable coeffs(deg+1);
    if deg >= 2,
        ((deg:-1:2).*(deg-1:-1:1))'.*coeffs(1:end-2,:) == nonneg_poly_coeffs(deg-2,mm);
    end
cvx_end_set

% Copyright 2007 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
