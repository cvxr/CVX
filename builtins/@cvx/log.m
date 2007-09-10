function y = log( x )

% LOG    Natural logarithm.
%   LOG(X) is the natural logarithm of the elements of X.
%   Complex results are produced if X is not positive.
%
%   Disciplined geometric programming information:
%       LOG(X) currently supports only constants, monomials, and
%       posynomials. Support for affine and/or concave expressions
%       will come in a later revision.

%
% Determine the expression types
%

persistent remap
if isempty( remap ),
    remap_1 = cvx_remap( 'constant' );
    remap_2 = cvx_remap( 'real-affine' ) & ~remap_1;
    remap_3 = cvx_remap( 'monomial' );
    remap_4 = cvx_remap( 'posynomial' );
    remap   = remap_1 + 2 * remap_2 + 3 * remap_3 + 4 * remap_4;
end
v = remap( cvx_classify( x ) );

%
% Process each type of expression one piece at a time
%

vu = unique( v );
nv = length( vu );
if nv ~= 1,
    y = cvx( size( x ), [] );
end
for k = 1 : nv,

    %
    % Select the category of expression to compute
    %

    vk = vu( k );
    if nv == 1,
        xt = x;
    else
        t = v == vk;
        xt = cvx_subsref( x, t );
    end

    %
    % Perform the computations
    %

    switch vk,
        case 0,
            % Invalid
            error( sprintf( 'Disciplined convex programming error:\n    Illegal operation: log( {%s} ).', cvx_class( xt, true, true, true ) ) );
        case 1,
            % Constant
            yt = log( cvx_constant( xt ) );
        case 2,
            % Affine (invalid)
            error( sprintf( 'Disciplined convex programming error:\n    Logarithms of affine expressions are not yet supported.' ) );
        case 3,
            % Monomial
            nb = prod( xt.size_ );
            [ rx, cx, vx ] = find( xt.basis_ );
            nx = length( rx );
            global cvx___
            logs = cvx___.logarithm( rx, 1 );
            tt = vx ~= 1; nt = sum( tt );
            bx = sparse( [ ones( nt, 1 ) ; logs ], [ cx( tt ) ; cx ], [ log( vx( tt ) ) ; ones( nb, 1 ) ], max( logs ), size( xt.basis_, 2 ) );
            yt = cvx( xt.size_, bx );
        case 4,
            % Posynomial
            global cvx___
            sx = xt.size_;
            xt = xt.basis_;
            rc = sum( xt ~= 0, 1 );
            ru = unique( rc );
            nu = length( ru );
            if nu ~= 1,
                yt = cvx( sx, [] );
            end
            for k = 1 : nu,
                rk = ru( k );
                if nu == 1,
                    xtt = xt;
                else
                    tt  = rc == rk;
                    xtt = xt( :, tt );
                end
                [ rx, cx, vx ] = find( xtt );
                rx = rx( : ); cx = cx( : ); vx = vx( : );
                nq = length( vx );
                vx = log( vx );
                tz = rx ~= 1;
                rx = cvx___.logarithm( rx( tz ), 1 );
                vx = vx + cvx( nq, sparse( rx, find( tz ), 1, max( rx ), nq ) );
                vx = reshape( vx, rk, nq / rk );
                vx = logsumexp_sdp( vx, 1, cvx___.gptol );
                if nu == 1,
                    yt = reshape( vx, sx );
                else
                    yt = cvx_subsasgn( yt, tt, vx );
                end
            end
        otherwise,
            error( 'Shouldn''t be here.' );
    end

    %
    % Store the results
    %

    if nv == 1,
        y = yt;
    else
        y = cvx_subsasgn( y, t, yt );
    end

end

% Copyright 2007 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
