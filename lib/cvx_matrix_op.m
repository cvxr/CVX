function y = cvx_matrix_op( p, args )

if ischar( args{1} ) && isequal(args{1},'test'),
    y = run_test(p);
    return
end

narg = length(args);
if narg < p.nargs(end),
    if narg < p.nargs(1),
        cvx_throw( 'Not enough arguments.' );
    else
        args{end+1:p.nargs(end)} = [];
    end
elseif narg > p.nargs(end),
    cvx_throw( 'Too many arguments.' );
end
if ~isempty(p.args),
    [ args{:} ] = p.args( args{:} );
end

global cvx___
narg = length(args);
mult = false(1,narg);
cz = true;
oz = false;
sm = [];
nm = 1;
need_bc = false;
for k = 1 : narg,
    xt = args{k};
    if isa( xt, 'cvx' ),
        oz = true;
        if cvx_isconstant( xt ),
            xt = cvx_constant( xt );
            args{k} = xt;
        else
            cz = false;
        end
    end
    sx = size( xt );
    if ~cvx_isaffine( xt ),
        if narg == 1,
            cvx_dcp_error( p.name, 'unary', xt );
        else
            cvx_dcp_error( p.name, 'binary', args{1:2} );
        end
    end
    sa = sx(3:end);
    na = prod( sa );
    if na > 1,
        mult(k) = true;
        if nm == 1,
            nm = na;
            sm = sa;
        elseif ~isequal( sa, sm ),
            if cvx___.broadcast,
                need_bc = true;
                sa(end+1:length(sm)) = 1;
                sm(end+1:length(sa)) = 1;
                if any(sa(sa~=sm)~=1),
                    cvx_throw( 'N-D dimensions are not compatible.' );
                end
            else
                cvx_throw( 'N-D dimensions are not compatible.' );
            end
        end
    end
end

xt = args{1};
sx = size( xt );
sm(end+1:2) = 1;
if ~all( sx ),
    if ~all( sx(1:2) ),
        if isempty( p.empty ),
            cvx_throw( 'Matrix input must not be empty.' );
        else
            y = repmat( p.empty, sm );
        end
    else
        y = zeros( sa );
        if oz, y = cvx( y ); end
    end
    return
end
if ~isempty( p.structure ),
    switch p.structure,
        case { 'square', 'symm', 'eig', 'chol', 'psdeig' },
            if sx(1) ~= sx(2),
                cvx_throw( 'Matrix must be square.' );
            end
    end
    switch p.structure,
        case { 'symm', 'eig', 'psdchol', 'psdeig' },
            xt  = args{1};
            xtt = permute( xt, [2,1,3:length(sx)] );
            if ~isreal( xtt ), xtt = conj( xtt ); end
            err = xt - xtt;
            xt  = 0.5 * ( xt + xtt );
            err = max( sum(abs(cvx_basis(err)),1) ./ max(sum(abs(cvx_basis(xt)),1)) );
            if err > 2 * eps,
                cvx_throw( [ ...
                    'Matrix input %d expected to be symmetric (deviation: %g).\n', ...
                    '--- If this number is small (<1e-6), it may simply be due to roundoff error.\n', ...
                    '    This can be corrected by applying the SYM(X) function.\n', ...
                    '--- Otherwise, this is likely due to a modeling error. Did you declare the\n', ...
                    '    relevant matrix variables to be "symmetric" or "hermitian"?' ], k, full(err) );
            end
            args{1} = xt;
    end
    if cz,
        err = 0;
        success = true;
        sa = sx(3:end);
        na = prod(sa);
        switch p.structure,
            case { 'eig', 'psdeig' },
                zt = zeros([sx(1),1,sa]);
                for a = 1 : na,
                    zt(:,a) = eig(full(xt(:,:,a)));
                end
                if strncmp( p.structure, 'psd', 3 ),
                    err = max( -min(zt(:,:),1) ./ max(zt(:,:),1) );
                    success = err < 2 * eps;
                end
                args{1} = zt;
            case 'chol';
                for a = 1 : na,
                    xtt = xt(:,:,a);
                    [ R, q ] = chol(xtt);
                    if q ~= 0, 
                        err = max( err, norm( R' * R - xtt, 'fro' ) / norm( xtt, 'fro' ) );
                    end
                    xt(:,:,a) = R;
                end
                args{1} = xt;
                success = err < 2 * sx(1) * eps;
            case 'svd',
                zt = zeros([min(sx(1:2)),1,sa]);
                for a = 1 : na,
                    zt(:,a) = svd(full(xt(:,:,a)));
                end
                args{1} = zt;
        end
        if ~success,
            cvx_throw( [ ...
                'Matrix input %d expected to be positive definite (deviation: %g).\n', ...
                '--- If this number is small (<1e-8), it may simply be due to roundoff error,\n', ...
                '    which can be rectified by adding a small positive value to the diagonal.' ], k, err );
        end
    end
end

if need_bc,
    for a = 1 : narg,
        if mult(a),
            xt = args{xt};
            sx = size( xt );
            sx(end+1:length(sa)+2) = 1;
            sx([true,true,sx(3:end)==sa(3:end)]) = 1;
            if any( sx ~= 1 ),
                args{xt} = repmat( xt, sx );
            end
        end
    end
end

if cz,
    func = p.constant;
elseif ~isempty( p.diagonal ),
    func = p.affine;
    xt = args{1};
    sx = size( xt );
    sd = min(sx(1),sx(2));
    nnzX = nnz( xt );
    if nnzX <= sx(1) * prod(sx(3:end)),
        dX = 1:sx(1)+1:sx(1)*sx(2);
        if mult(1),
            dX = bsxfun(@plus,dX',0:sx(1)*sx(2):prod(sx)-1);
        end
        dX = reshape( cvx_subsref( xt, dX ), [sd,1,sx(3:end)] );
        if nnz( dX ) == nnzX,
            args{1} = dX;
            func = p.diagonal;
        end
    end
else
    func = p.affine;
end

if nm == 1,
    
    y = func( args{:} );
    if oz, y = cvx( y ); end
    
else
    
    if oz, y = cvx( sm, [] );
    else y = zeros( sm ); end
    targs = args;
    for k = 1 : nm,
        for a = find(mult)',
            targs{a} = cvx_subsref( args{a}, ':', ':', k ); 
        end
        yt = func( targs{:} );
        y = cvx_subsasgn( y, k, yt );
    end
    
end

function err = run_test( p )
if ~isfield( p, 'test' ) || isempty( p.test ),
    if max(p.nargs) ~= 1,
        cvx_throw( 'Explicit test generator required for this function.' );
    end
    err = zeros(1,10);
    for k = 1 : 10,
        symm = true;
        switch p.structure,
        case 'eig',
            X = randn(25,25);
            X = X + X';
        case { 'chol', 'psdeig' },
            X = randn(25,25);
            X = X * X';
        otherwise,
            X = randn(25,10);
            symm = false;
        end
        X = X / norm( X, 'fro' );
        y1 = cvx_matrix_op( p, { X } );
        cvx_begin quiet
            cvx_solver sedumi
            cvx_precision best
            if symm,
                variable XX(size(X)) symmetric
            else
                variable XX(size(X))
            end
            z = cvx_matrix_op( p, { XX } );
            if cvx_isconcave( z ),
                maximize( z );
            else
                minimize( z );
            end
            XX == X; %#ok
        cvx_end
        err(1,k) = abs( cvx_optval - y1 );
        fprintf( 'Test %2d: true=%g, computed=%g, err=%g\n', k, y1, cvx_optval, err(1,k) );
        X2(:,:,k) = X; %#ok
        y2(k) = y1; %#ok
    end
    cvx_begin quiet
        cvx_solver sedumi
        cvx_precision best
        if symm,
            variable XX(size(X2)) symmetric
        else
            variable XX(size(X2))
        end
        z = cvx_matrix_op( p, { XX } );
        if cvx_isconcave( z ),
            maximize( sum(z) );
        else
            minimize( sum(z) );
        end
        XX == X2; %#ok
    cvx_end
    err2 = abs( z(:) - y2(:) );
    fprintf( 'Combined test: min %g, med %g, max %g\n', min(err2), median(err2), max(err2) );
    err(2,:) = err2'; 
end

