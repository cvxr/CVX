function v = cvx_pushexpr( args )

global cvx___
try
    pstr = cvx___.problems(end);
catch
    cvx_throw( 'No CVX model is present.' );
end

name = args(1).name;
if ~isvarname( name ),
    cvx_throw( 'Invalid expression name: %s', name );
elseif isfield( pstr.variables, name ),
    cvx_throw( 'Primal variable name conflict: %s', name );
elseif isfield( pstr.duals,name ),
    cvx_throw( 'Dual variable name conflict: %s', name );
end

%
% Parse the structure
%

xsiz = args(1).args;
if isempty( xsiz ),
    xsiz = [ 1, 1 ];
end

%
% Create the variables
%

v = cvx( xsiz, [] );



