function cvx_erase( p )

global cvx___
pstr = cvx___.problems( p );
nf = pstr.n_variable + 1;
ne = pstr.n_equality + 1;
if nf <= 2,
    cvx___.classes     = int8(3);
    cvx___.cones       = struct( 'type', {}, 'indices', {}, 'slacks', {} );
    cvx___.exponential = [];
    pstr.n_variable    = 0;
elseif length( cvx___.classes ) >= nf,
    cvx___.classes( nf : end, : ) = [];
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
    end
end
if ~any( cvx___.exponential ),
    cvx___.exponential = zeros(0,1);
    cvx___.logarithm = zeros(0,1);
end
if ne <= 1,
    cvx___.equalities = {};
    cvx___.inequality = false(0,1);
    cvx___.n_equality = 0;
elseif length( cvx___.equalities ) >= ne,
    cvx___.equalities( ne : end ) = [];
    cvx___.inequality( ne : end ) = [];
    cvx___.n_equality = sum( cellfun( @(x)size(x,2), cvx___.equalities ) );
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
