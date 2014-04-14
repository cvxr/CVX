function y = matrix_op( params, X, varargin )

sx = size( X );
ox = isa( X, 'cvx' );
if params.square(1) && sx(1) ~= sx(2),
    error( 'CVX:ArgError', 'Matrix must be square.' );
elseif ~ox,
    cx = true;
elseif ~cvx_isaffine( X )
    error( 'CVX:ArgError','Input must be affine.' );
else
    cx = cvx_isconstant( X );
end

narg = nargin - 2;
sm = sx( 3 : end );
nm = prod( sm );
xm = nm > 1;
ym = false(1,narg);
for a = 1 : narg,
    ms = length(params.square) > 1;
    arg = varargin{a};
    sa = size( arg );
    if ms && params.square(a+1) && sa(1) ~= sa(2),
        error( 'CVX:ArgError', 'Argument %d must be square.', a+1 );
    end
    if isa( arg, 'cvx' ),
        ox = true;
        cx = cx && cvx_isconstant( arg );
        if ~cvx_isaffine( arg ),
            error( 'CVX:ArgError', 'Input must be affine.' );
        end
    end
    sa = sa( 3 : end );
    na = prod( sa );
    if na > 1,
        ym(a) = true;
        if nm == 1, 
            nm = na;
            sm = sa; 
        elseif ~isequal( sa, sm ),
            error( 'CVX:ArgError', 'Dimensions are not compatible.' );
        end
    end
end

if cx,
    func = params.funcs{1};
elseif isreal( X ),
    func = params.funcs{2};
else
    func = params.funcs{end};
end

if ox && cx,
    X = cvx_constant( X );
    varargin = cellfun( @cvx_constant, varargin, 'UniformOutput', false );
end

if nm == 1,
    
    y = func( X, varargin{:} );
    
elseif cx,
    
    y = zeros(nm,1);
    args = varargin;
    for k = 1 : nm,
        if xm, xt = X(:,:,k); end
        for a = 1 : narg,
            if ym(a), args{a} = varargin{a}(:,:,k); end
        end
        y(k) = func( xt, args{:} );
    end
    
else
    
    y = cvx( sx(3:end), [] );
    xt = X; 
    args = varargin;
    for k = 1 : nm,
        if xm, 
            xt = cvx_subsref( X, ':', ':', k );
        end
        for a = 1 : narg,
            if ym(a), args{a} = cvx_subsref( varargin{a}, ':', ':', k ); end
        end
        yt = func( xt, args{:} );
        y.basis_(1:size(yt.basis_,1),k) = yt.basis_;
    end
    
end

if ox && cx,
    y = cvx( y );
end
