function y = pow_cvx( x, p, mode )

%POW_CVX   Internal cvx version.

%
% Expression type matrices
%

persistent remap_pos remap_p remap_abs remap_pwr remap_uniq
if isempty( remap_uniq ),
    % vertical dim: {<0,=0,0<p<1,=1,>1,>1 int,>1 even int,2}
    remap_x1  = cvx_remap( 'real' );
    remap_x1c = cvx_remap( 'complex' );
    remap_x1a = remap_x1 | remap_x1c;
    remap_x1x = ~ remap_x1a;
    
    remap_x2  = cvx_remap( 'l-valid'  ) & remap_x1x;
    remap_x2x = ~remap_x2 | remap_x1x;
    
    remap_x3n = cvx_remap( 'n-affine' ) & remap_x2x;
    remap_x3  = cvx_remap( 'r-affine' ) & remap_x2x;
    remap_x3p = cvx_remap( 'p-affine' ) & remap_x2x;
    remap_x3r = remap_x3n | remap_x3 | remap_x3p;
    remap_x3c = cvx_remap( 'c-affine' )  & ~remap_x1c;
    remap_x3a = remap_x3r | remap_x3c;
    remap_x3x = ~ remap_x3a | remap_x2x;
    
    remap_x4n = cvx_remap( 'n-concave' ) & remap_x3x;
    remap_x4  = cvx_remap( 'concave' )   & remap_x3x;
    remap_x4p = cvx_remap( 'p-concave' ) & remap_x3x;
    remap_x4a = remap_x4n | remap_x4 | remap_x4p;
    remap_x4x = ~ remap_x4a | remap_x3x;
    
    remap_x5n = cvx_remap( 'n-convex' ) & remap_x4x;
    remap_x5  = cvx_remap( 'convex'   ) & remap_x4x;
    remap_x5p = cvx_remap( 'p-convex' ) & remap_x4x;
    
    %
    % That's right, ladies and gentlemen, we have 25 different branches
    % to deal with here; 2 map to the same implementation.
    %
    
    % POW_P( X, P ) : X ^ P if X >  0, +infty otherwise  (P<0)   : cvx/decr
    %               : 1     if X >= 0, +-infty otherwise (P=0)   : cvx/decr, ccv/incr
    %               : X ^ P if X <  0, -infty otherwise  (0<P<1) : ccv/incr
    %               : X     if X >= 0, +-infty otherwise (P=1)   : cvx/sign, ccv/incr
    %               : X ^ P if X >  0, +infty otherwise  (P>1)   : cvx/sign
    
    remap_p   = [1 ;1 ; 1; 1; 1; 1; 1; 1] * remap_x1a + ...
                [2 ;2 ; 2; 2; 2; 2; 2; 2] * remap_x2  + ...
                [3 ;5 ; 8;10;10;10;10;10] * ( remap_x3n | remap_x4n ) + ...
                [4 ;6 ; 9;11;13;13;13;13] * remap_x3  + ...
                [4 ;6 ; 9;11; 0; 0; 0; 0] * remap_x4  + ...
                [4 ;7 ; 9;12;14;14;14;24] * remap_x3p + ...
                [4 ;7 ; 9;12; 0; 0; 0; 0] * remap_x4p + ...
                [0 ;7 ; 0;12;14;14;14;24] * remap_x5p;

    % POW_POS( X, P ) : X ^ P if X >= 0, 0 otherwise (P>=1) : cvx/incr
    % This is the only branch by definition.
            
    remap_pos = [0 ;0 ; 0; 1; 1; 1; 1; 1] * remap_x1a + ...
                [0 ;0 ; 0; 2; 2; 2; 2; 2] * remap_x2  + ...
                [0 ;0 ; 0;15;15;15;15;15] * ( remap_x3n | remap_x5n ) + ...
                [0 ;0 ; 0;16;17;17;17;25] * remap_x3  + ...
                [0 ;0 ; 0;12;14;14;14;24] * remap_x3p + ...
                [0 ;0 ; 0;16;14;14;14;25] * remap_x5  + ...
                [0 ;0 ; 0;12;14;14;14;24] * remap_x5p;
            
    % POW_ABS( X, P ) : |X| ^ P (P>=1) : cvx/sign
    % This is the only branch by definition.
            
    remap_abs = [0 ;0 ; 0; 1; 1; 1; 1; 1] * remap_x1a + ...
                [0 ;0 ; 0; 2; 2; 2; 2; 2] * remap_x2  + ...
                [0 ;0 ; 0;18;19;19;19;19] * remap_x3  + ...
                [0 ;0 ; 0;18;20;20;20;20] * remap_x3c + ...
                [0 ;0 ; 0;18;14;14;14;24] * ( remap_x5p | remap_x3p | remap_x3n | remap_x4n );

    % POWER( X, P ) : P < 0: cvx/decr if X > 0, ccv/decr if X < 0
    %               : P = 0: 1
    %               : POW_P( X, P ) if 0 < P < 1
    %               : P = 1: X
    %               : P > 1: cvx/incr if X > 0
    %               : P > 1, odd int: cvx/incr if X > 0, ccv/incr if X < 0
    %               : P > 1, even int: cvx/incr if X > 0, cvx/decr if X < 0
    
    remap_pwr = [0 ;21; 0;12; 0; 0; 0; 0] * cvx_remap( 'valid' ) + ...
                [1 ; 0; 1; 0; 1; 1; 1; 1] * remap_x1a + ...
                [2 ; 0; 2; 0; 2; 2; 2; 2] * remap_x2  + ...
                [22; 0; 8; 0; 0;23;14;24] * remap_x3n + ...
                [ 0; 0; 9; 0; 0; 0;19;24] * remap_x3  + ...
                [ 4; 0; 9; 0;14;14;14;24] * remap_x3p + ...
                [ 0; 0; 8; 0; 0;23;14;24] * remap_x4n + ...
                [ 0; 0; 9; 0; 0; 0; 0; 0] * remap_x4  + ...
                [ 4; 0; 9; 0; 0; 0; 0; 0] * remap_x4p + ...
                [22; 0; 0; 0; 0; 0; 0; 0] * remap_x5n + ...
                [ 0; 0; 0; 0;14;14;14;24] * remap_x5p;
            
    remap_uniq = logical( [0,0,0,1,0,1,0,0,0,0,1,1,1] );
end

%
% Argument check
%

error(nargchk(3,3,nargin)); %#ok
p = cvx_constant( p );
if nnz( isinf( p ) | isnan( p ) ),
    error( 'Second argument must be Inf or NaN.' );
end
if ~ischar( mode ) || size( mode, 1 ) ~= 1,
    error( 'Third argument must be a string.' );
end
switch mode,
    case 'power',   remap = remap_pwr;
    case 'pow_p',   remap = remap_p;
    case 'pow_pos', remap = remap_pos;
    case 'pow_abs', remap = remap_abs;
    otherwise, error( [ 'Invalid power mode: ', mode ] );
end

%
% Check sizes
%

sx = size( x ); xs = all( sx == 1 );
sp = size( p ); ps = all( sp == 1 );
if xs,
    sy = sp;
elseif ps || isequal( sx, sp ),
    sy = sx;
else
    error( 'Matrix dimensions must agree.' );
end

%
% Determine the expression types
%

v = 1 + ( p >= 0 ) + ( p > 0 ) + ( p >= 1 ) + ( p > 1 ) .* ( 1 + ( rem( p, 1 ) == 0 ) + ( rem( p, 2 ) == 0 ) + ( p == 2 ) );
v = remap( v(:)' + size(remap,1) * ( cvx_classify( x ) - 1 ) );
if ~ps,
    t = remap_uniq( v + 1 );
    if any( t ),
        [ pk, pi, pj ] = unique( p( t ) ); %#ok
        vt = v( t );
        v( t ) = vt + ( reshape( pj, size(vt) ) - 1 ) / length( pk );
    end
end
    
%
% Process each type of expression one piece at a times
%

xt = x;
pt = p;
vu = sort( v(:) );
vu = vu([true;diff(vu)~=0]);
nv = length( vu );
if nv ~= 1,
    if cvx_isconstant( x ),
        y = zeros( sy );
    else
        y = cvx( sy, [] );
    end
end
for k = 1 : nv,
    
    %
    % Select the category of expression to compute
    %
    
    vk = vu( k );
    if nv ~= 1,
        t = v == vk;
        if ~xs, xt = cvx_subsref( x, t ); end
        if ~ps, pt = cvx_subsref( p, t ); end
    end
    sz = size(xt);
    nd = length(sz)+1;
    sw = sz;
    sw(nd) = 2;
    
    %
    % Perform the computations
    %
    
    vk = floor(vk);
    switch vk,
        case 0,
            % Invalid
            pt = unique( pt );
            if length(pt) > 1,
                pt = sprintf( '%g,', unique(pt) );
                pt = [ '{', pt(1:end-1), '}' ];
            else
                pt = sprintf( '%g', pt );
            end
            if isequal( mode, 'power' ),
                error( 'Disciplined convex programming error:\n    Illegal operation: {%s} .^ %s\n    (Consider POW_P, POW_POS, or POW_ABS instead.)', cvx_class( xt, true, true ), pt ); 
            else
                error( 'Disciplined convex programming error:\n    Illegal operation: %s( {%s}, %s )', mode, cvx_class( xt, true, true ), pt ); 
            end
            
        case 1,
            % Constant result
            yt = cvx( feval( mode, cvx_constant( xt ), pt ) );
            
        case 2,
            % Log-convex/affine/concave
            yt = exp( log( xt ) .* pt );
            
        case 3,
            % pow_p( n-concave, p < 0 )
            error( 'Disciplined convex programming error:\n    Trivially infeasible: %s( {%s}, %g )', mode, cvx_class( xt, true, true ), pt ); 
            
        case 4,
            % 4: pow_p( concave, p < 0 )
            yt = [];
            if vk == 22, xt = -xt; end
            cvx_begin
                epigraph variable yt(sz)
                { cat( nd, cvx_accept_concave(xt), yt ), 1 } == geo_mean_cone( sw, nd, [-pt,1], 'func' ); %#ok
                cvx_setnneg(yt);
            cvx_end
            if vk == 22, yt = -yt; end
            
        case 5,
            % pow_p( n-concave, p = 0 )
            warning( 'Disciplined convex programming warning:\n    Almost certainly infeasible: %s( {%s}, %g )', mode, cvx_class( xt, true, true ), pt );
            yt = ones(sz);
            cvx_begin
                xt >= 0; %#ok
            cvx_end
            
        case 6,
            % pow_p( concave, p = 0 )
            yt = ones(sz);
            cvx_begin
                xt >= 0; %#ok
            cvx_end
            
        case {8,10},
            % pow_p( n-concave, p > 0 )
            warning( 'Disciplined convex programming warning:\n    Almost certainly infeasible: %s( {%s}, %g )', mode, cvx_class( xt, true, true ), pt );
            yt = zeros(sz);
            cvx_begin
                xt >= 0; %#ok
            cvx_end
            
        case 9,
            % pow_p( concave, 0 < p < 1 )
            yt = [];
            cvx_begin
                hypograph variable yt(sz)
                { cat( nd, cvx_accept_concave(xt), ones(sz) ), yt } == geo_mean_cone( sw, nd, [pt,1-pt], 'func' ); %#ok
                cvx_setnneg(yt);
            cvx_end
            
        case 11,
            % pow_p( concave, p = 1 )
            yt = xt;
            cvx_begin
                xt >= 0; %#ok
            cvx_end
            
        case 12,
            % pow_p( p-nonconst, p = 1 )
            % pow_pos( p-nonconst, p = 1 )
            % power( valid, p = 1 )
            yt = xt;
            
        case 13,
            % pow_p( affine, p > 1 )
            yt = [];
            cvx_begin
                epigraph variable yt(sz)
                { cat( nd, yt, ones(sz) ), xt } == geo_mean_cone( sw, nd, [1/pt,1-1/pt], 'pos' );  %#ok
                cvx_setnneg(yt);
            cvx_end
            
        case 14,
            % pow_p( p-convex, p > 1 )
            % pow_pos( p-convex, p > 1 )
            % pow_abs( n-concave, p > 1 )
            yt = [];
            cvx_begin
                epigraph variable yt(sz)
                { cat( nd, yt, ones(sz) ), cvx_accept_convex(abs(xt)) } == geo_mean_cone( sw, nd, [1/pt,1-1/pt], 'func' );  %#ok
                cvx_setnneg(yt);
            cvx_end
            
        case 15,
            % pow_pos( n-nonconst, p >= 1 )
            yt = zeros( sz );
            
        case 16,
            % pow_pos( convex, p = 1 )
            yt = max( xt, 0 );
            
        case 17,
            % pow_pos( convex, p > 1 )
            yt = [];
            cvx_begin
                epigraph variable yt(sz)
                variable zt(sz)
                { cat( nd, yt, ones(sz) ), zt } == geo_mean_cone( sw, nd, [1/pt,1-1/pt], 'pos' );  %#ok
                xt <= zt; %#ok
                cvx_setnneg(yt);
                cvx_setnneg(zt);
            cvx_end
            
        case 18,
            % pow_abs( affine/p-nonconst/n-nonconst, p = 1 )
            yt = abs( xt );
            
        case 19,
            % pow_abs( affine, p > 1 )
            % power( affine, p even )
            yt = [];
            cvx_begin
                epigraph variable yt(sz)
                { cat( nd, yt, ones(sz) ), xt } == geo_mean_cone( sw, nd, [1/pt,1-1/pt], 'abs' ); %#ok 
                cvx_setnneg(yt);
            cvx_end
            
        case 20,
            % pow_abs( c-affine, p > 1 )
            yt = [];
            cvx_begin
                epigraph variable yt(sz)
                { cat( nd, yt, ones(sz) ), xt } == geo_mean_cone( sw, nd, [1/pt,1-1/pt], 'cabs' ); %#ok 
                cvx_setnneg(yt);
            cvx_end
            
        case 21,
            % power( valid, p = 1 )
            yt = ones(sz);
            
        case 22,
            % pow( n-affine/n-convex, p < 0 )
            yt = [];
            cvx_begin
                hypograph variable yt(sz)
                { cat( nd, cvx_accept_concave(-xt), -yt ), 1 } == geo_mean_cone( sw, nd, [-pt,1], 'func' ); %#ok
                cvx_setnpos(yt);
            cvx_end
            
        case 23,
            % power( n-affine/n-concave, p odd integer )
            yt = [];
            cvx_begin
                hypograph variable yt(sz)
                { cat( nd, -yt, ones(sz) ), -cvx_accept_concave(xt) } == geo_mean_cone( sw, nd, [1/pt,1-1/pt], 'abs' ); %#ok 
                cvx_setnpos(yt);
            cvx_end

        case 24,
            % power( affine, 2 )
            % power/pow_abs/pow_p( p-convex, 2 )
            % power/pow_abs( n-concave, 2 )
            yt = square( xt );

        case 25,
            % pow_pos( affine, 2 )
            yt = square_pos( xt );
    end
    
    %
    % Store the results
    %
    
    if nv ~= 1,
        y = cvx_subsasgn( y, t, yt );
    elseif isequal(sz,sy),
        y = yt;
    else
        y = yt * ones(sy);
    end
    
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.

