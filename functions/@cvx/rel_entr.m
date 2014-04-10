function z = rel_entr( x, y )

%REL_ENTR   Internal cvx version.

cvx_expert_check( 'rel_entr', x, y );

persistent remap funcs
if isempty( remap ),
    remap = cvx_remap( ...
        { { 'nonnegative' }, { 'real' } }, ...
        { { 'r_affine' }, { 'concave' } } );
    funcs = { @rel_entr_1, @rel_entr_2 };
end

try
    z = binary_op( 'rel_entr', funcs, remap, x, y );
catch exc
    throw( exc );
end

function z = rel_entr_1( x, y )            
% Constant
z = cvx( rel_entr( cvx_constant( x ), cvx_constant( y ) ) );

function z = rel_entr_2( x, y )
z = [];
sz = max( numel(y), numel(x) );
cvx_begin
    epigraph variable z( sz );
    { -z, x, y } == exponential( sz ); %#ok
cvx_end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
