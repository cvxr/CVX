function z = newsubst( prob, src )
error( nargchk( 2, 2, nargin ) );

%
% Check problem
%

if ~isa( prob, 'cvxprob' ),
    error( 'First argument must be a cvxprob object.' );
end
global cvx___
p = index( prob );

%
% Handle trivial cases
%

switch class( src ),
    case { 'double', 'sparse' },
        z = cvx( prob, src );
        return
    case 'cvx',
        oprob = problem( src );
        if oprob == prob,
            z = src;
            return
        elseif cvx_isconstant( src ),
            z = cvx( prob, cvx_constant( src ) );
            return
        end
    otherwise,
        error( [ 'Cannot handle an object of type "', class( src ), '" here.' ] );
end

%
% Insure that the object is from the immediate parent
%

op = index( oprob );
pos = cvx___.problems( p ).stackpos;
opos = cvx___.problems( op ).stackpos;
if opos >= pos,
    error( 'Internal cvx data corruption.\n    Please type ''cvx_clear'' and start over.' );
elseif opos + 1 ~= pos,
    for k = opos + 1 : pos - 1,
        src = newsubst( cvx___.stack{k}, src );
    end
    oprob = problem( src );
end

%
% Add to the substitution table
%

sz = size( src );
bO = cvx_basis( cvx___.problems( p ).substitutions );
bN = cvx_basis( src );
tc = sum( bN ~= 0, 2 ) == ( bN( :, 1 ) ~= 0 );
need_const = any( tc );
if need_const,
    bC = bN( tc, 1 );
    bN( tc, : ) = 0;
end
[ mO, nO ] = size( bO );
[ mN, nN ] = size( bN );
if mO ~= 0,
    nZ = max( nO, nN );
    if nO < nZ, bO( :, nZ ) = 0; end
    if nN < nZ, bN( :, nZ ) = 0; end
    [ bNL, bN ] = cvx_bcompress( [ bO ; bN ], false, mO );
else,
    [ bNL, bN ] = cvx_bcompress( bN );
end
bN = cvx( oprob, size( bN, 1 ), bN );
cvx___.problems( p ).substitutions = bN;
if mO ~= 0, 
    bNL = bNL( mO + 1 : end, : );
    bN  =  bN( mO + 1 : end ); 
end
vars = cvx___.problems( p ).variables;
if ~isempty( bN ),
    dmap = eval( 'vars.map_', '[]' );
    ndx  = length( cvx___.problems( p ).vexity );
    nvar = newvar( prob, '', size( bN ) );
    if isempty( dmap ),
        vars.map_ = nvar;
    else,
        vars.map_ = [ vars.map_ ; nvar ];
    end
    cvx___.problems( p ).vexity( ndx + 1 : end ) = cvx_vexity( bN );
    cvx___.problems( p ).variables = vars;
end

%
% Create the substitution variable
%

z = bNL * cvx_basis( vars.map_ );
if need_const, z( tc, 1 ) = bC; end
z = cvx( prob, sz, z );

%
% Clear the 'complete' flag
%

cvx___.problems( p ).complete = false;

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
