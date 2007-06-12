function coeffs = nonneg_poly_coeffs( deg, trig, mm )
error( nargchk( 1, 3, nargin ) );

 
%
% Check degree argument
%

if ~cvx_check_dimension( deg, true ),
    error( 'Argument must be a nonnegative integer.' );
end

%
% Check trig argument
%

switch nargin,
    case 1,
        mm = [];
        trig = false;
    case 2,
        if isequal( trig, 'trig' ),
            trig = true;
            mm = [];
        else
            mm = trig;
            trig = false;
        end
    case 3,
        if ~isequal( trig, 'trig' ),
            error( 'Second argument must be a range, or ''trig''.' );
        end
        trig = true;
end

%
% Check range argument
%
    
if isempty( mm ),
    mm = [ -Inf, +Inf ];
elseif ~isa( mm, 'double' ) | ~isreal( mm ) | numel( mm ) ~= 2,
    error( 'Second argument, if supplied, must be a range [ xmin xmax ].' );
elseif any( mm(1) == mm(2) & isinf( mm(1) ) ),
    error( 'Intervals [-Inf,-Inf] and [Inf,Inf] are not accepted.' );
end

%
% Construct set
%

if trig,
    
    %
    % Trigonometric
    %
    
    cvx_begin set
        variable coeffs(deg+1) complex;
        if mm(1) == mm(2),
            % Positive at a single point
            polyval_trig( coeffs, mm(1) ) >= 0;
        elseif mm(2) > mm(1) + 2 * pi,
            % Positive over the entire unit circle
            [ii,jj,vv] = find(hermitian_semidefinite(deg+1));
            coeffs == sparse( deg+1-abs(ii-jj), 1, vv );
        elseif mm(2) > mm(1),
            % Positive over a subset of the unit circle
            a = exp( j * ( 0.5 * ( mm(2) + mm(1) ) ) );
            b = cos( 0.5 * ( mm(2) - mm(1) ) );
            coeffs1 = nonneg_poly_coeffs(deg,'trig');
            coeffs2 = nonneg_poly_coeffs(deg-1,'trig');
            coeffs == coeffs1 ...
                + [ (0.5*a)*coeffs2 ; 0 ] ...
                + [ 0 ; 0 ; (0.5*conj(a))*coeffs2(1:end-1) ] ...
                + [ zeros(deg-1,1) ; (0.5*a)*conj(coeffs2(end)) ; 0 ] ...
                - [ 0 ; b*coeffs2 ];
        end
    cvx_end
    
else        
    
    %
    % Non-trigonometric
    %

    cvx_begin set
        variable coeffs(deg+1);
        if mm(1) == mm(2),
            % Positive at a single point
            polyval( coeffs, mm(1) ) >= 0;
        elseif mm(1) == -Inf,
            isodd = rem(deg,2);
            deg2 = floor(0.5*deg) + 1;
            [ii,jj,vv] = find(semidefinite(deg2));
            coeffs1 = sparse(ii+jj-1+isodd,1,vv);
            if mm(2) == +Inf,
                % [ -Inf, +Inf ]
                coeffs == coeffs1;
            else
                % [ -Inf, mm(2) ]
                coeffs2 = nonneg_poly_coeffs(deg-1);
                coeffs == coeffs1 - [ coeffs2 ; 0 ] + [ 0 ; mm(2) * coeffs2 ];
            end
        elseif mm(2) == +Inf,
            % [ mm(1), +Inf ]
            coefff1 = nonneg_poly_coeffs(deg);
            coeffs2 = nonneg_poly_coeffs(deg-1);
            coeffs == coeffs1 + [ coeffs2 ; 0 ] - [ 0 ; mm(1) * coeffs2 ];
        elseif mm(1) < mm(2),
            % [ mm(1), mm(2) ]
            coeffs1 = nonneg_poly_coeffs(deg);
            coeffs2 = nonneg_poly_coeffs(deg);
            [ 0 ; coeffs ] == ...
                + [ coeffs1 ; 0 ] - [ 0 ; mm(1) * coeffs1 ] ...
                - [ coeffs2 ; 0 ] + [ 0 ; mm(2) * coeffs2 ];
        end
    cvx_end

end

% Copyright 2007 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.