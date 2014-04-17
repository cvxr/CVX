function erase( p )

global cvx___
p = p.index_;
prob = cvx___.problems( p );
if prob.cleared, return; end

nf = length( prob.t_variable ) + 1;
ne = prob.n_equality + 1;
nl = prob.n_linform + 1;
nu = prob.n_uniform + 1;
if nf <= 2,
    cvx___.classes = int8(3);
    cvx___.canslack = false;
    cvx___.readonly = 0;
    cvx___.cones = struct( 'type', {}, 'indices', {} );
    cvx___.exponential = [];
    cvx___.logarithm   = [];
elseif length( cvx___.classes ) >= nf,
    cvx___.classes( nf : end ) = [];
    cvx___.canslack( nf : end ) = [];
    cvx___.readonly( nf : end ) = [];
    if ~isempty( cvx___.cones ),
        tt = true( 1, length( cvx___.cones ) );
        for k = 1 : length( cvx___.cones ),
            cvx___.cones( k ).indices( :, any( cvx___.cones( k ).indices >= nf, 1 ) ) = [];
            if isempty( cvx___.cones( k ).indices ), tt( k ) = false; end
        end
        cvx___.cones = cvx___.cones( 1, tt );
    end
    if ~isempty( cvx___.exponential ),
        cvx___.exponential( nf : end ) = [];
        cvx___.logarithm( nf : end ) = [];
        if ~any( cvx___.exponential ),
            cvx___.exponential = [];
            cvx___.logarithm = [];
        end
    end
end
if nf <= 2 || ne <= 1,
    cvx___.equalities = cvx;
    cvx___.needslack = true( 0, 1 );
elseif cvx___.n_equality >= ne,
    cvx___.equalities( ne: end ) = [];
    cvx___.needslack( ne : end ) = [];
end
if nf <= 2 || nl <= 1,
    cvx___.linforms = cvx;
    cvx___.linrepls = cvx;
elseif length( cvx___.linforms ) >= nl,
    cvx___.linforms( nl : end ) = [];
    cvx___.linrepls( nl : end ) = [];
end
if nf <= 2 || nu <= 1,
    cvx___.uniforms = cvx;
    cvx___.unirepls = cvx;
elseif length( cvx___.uniforms ) >= nu,
    cvx___.uniforms( nu : end ) = [];
    cvx___.unirepls( nu : end ) = [];
end
cvx___.problems( p ).cleared = true;

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
