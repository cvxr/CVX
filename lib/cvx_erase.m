function cvx_erase( p )

global cvx___
pstr = cvx___.problems( p );
% n_variables, n_equalities, n_cones
cp = pstr.checkpoint + 1;
nf = cp(1);
if nf <= 2,
    cvx___.classes     = int8(3);
    cvx___.cones       = zeros(0,1);
    cvx___.exponential = zeros(0,1);
elseif length( cvx___.classes ) >= nf,
    cvx___.classes( nf : end, : ) = [];
    if ~isempty( cvx___.exponential ),
        cvx___.exponential( nf : end, : ) = [];
        cvx___.logarithm( nf : end, : ) = [];
        if ~any( cvx___.exponential ),
            cvx___.exponential = zeros(0,1);
            cvx___.logarithm = zeros(0,1);
        end
    end
end
ne = cp(2);
if ne <= 1,
    cvx___.equalities = {};
    cvx___.inequality = false(0,1);
    cvx___.n_equality = 0;
elseif length( cvx___.equalities ) >= ne,
    cvx___.n_equality = cvx___.n_equality - ...
        sum( cellfun( @(x)size(x,2), cvx___.equalities( ne : end ) ) );
    cvx___.equalities( ne : end ) = [];
    cvx___.inequality( ne : end ) = [];
end
nc = cp(3);
if nc <= 1,
    cvx___.cones = [];
else
    cvx___.cones( nc : end ) = [];
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
