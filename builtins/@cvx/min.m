function z = min( varargin )

%   Disciplined convex/geometric programming information:
%       MAX is convex, log-log-concave, and nondecreasing in its first
%       two arguments. Thus when used in disciplined convex programs, 
%       both arguments must be concave (or affine). In disciplined 
%       geometric programs, both arguments must be log-concave/affine.

if nargin == 2,
    
    %
    % min( X, Y )
    %

    persistent BP %#ok
    if isempty( BP ),
        BP.map = cvx_remap( ...
            { { 'real' } }, ...
            { { 'nonpositive', 'n_nonconst' }, { 'nonnegative', 'p_nonconst' } }, ...
            { { 'nonnegative', 'p_nonconst' }, { 'nonpositive', 'n_nonconst' } }, ...
            { { 'l_concave', 'positive' } }, ...
            { { 'concave' } } );
        BP.constant = 1;
        BP.funcs = { @min_b1, @min_b2, @min_b3, @min_b4, @min_b5 };
        BP.constant = [];
        BP.name = 'min';
    end
    z = cvx_binary_op( BP, varargin{:} );

elseif nargin > 1 && ~isempty( varargin{2} ),

    cvx_throw( 'MIN with two matrices to compare and a working dimension is not supported.' );
        
else
    
    %
    % min( X, [], dim )
    %

    persistent P %#ok
    if isempty( P ),
        P.map      = cvx_remap( { 'real' ; 'l_concave' ; 'concave' } );
        P.funcs    = { @min_r1, @min_r2, @min_r3 };
        P.constant = 1;
        P.zero     = [];
        P.reduce   = true;
        P.reverse  = false;
        P.name     = 'min';
        P.dimarg   = 3;
    end
    z = cvx_reduce_op( P, varargin{:} );
   
end 

function z = min_b1( x, y )
z = builtin( 'min', x, y );

function z = min_b2( x, y ) %#ok
% min( nonnegative, nonpositive ) (pass-through x)
z = x;

function z = min_b3( x, y ) %#ok
% min( nonpositive, nonnegative ) (pass-through y)
z = y;

function z = min_b4( x, y )
z = [];
sz = [max(numel(x),numel(y)),1]; %#ok
% geometric
cvx_begin gp
    hypograph variable z( sz )
    x >= z; %#ok
    y >= z; %#ok
cvx_end

function z = min_b5( x, y )
persistent nneg npos
if isempty( nneg ),
    nneg = cvx_remap( 'nonnegative', 'p_nonconst' );
    npos = cvx_remap( 'nonpositive', 'n_nonconst' );
end
% linear
z = [];
sz = [ max(numel(x),numel(y)), 1 ]; %#ok
vx = cvx_classify( x );
vy = cvx_classify( y );
nn = nneg(vx)&nneg(vy);
np = npos(vx)|npos(vy);
cvx_begin
    hypograph variable z( sz )
    x >= z; %#ok
    y >= z; %#ok
    cvx_setnneg( cvx_fastref( z, nn ) );
    cvx_setnpos( cvx_fastref( z, np ) );
cvx_end

function x =  min_r1( x, y ) %#ok
x = builtin( 'min', x, [], 1 );

function x = min_r2( x, y ) %#ok
nx = x.size_(1);
if nx > 1,
    nv = x.size_(2); %#ok
    cvx_begin gp
        hypograph variable z( 1, nv )
        x >= repmat( z, nx, 1 ); %#ok
    cvx_end
    x = cvx_optval;
end

function x = min_r3( x, y ) %#ok
persistent nneg npos
if isempty( nneg ),
    nneg = cvx_remap( 'nonnegative', 'p_nonconst' );
    npos = cvx_remap( 'nonpositive', 'n_nonconst' );
end
s = size( x );
if s(1) > 1,
    vx = cvx_classify( x );
    nn = all(reshape(nneg(vx),s),1);
    np = any(reshape(npos(vx),s),1);
    cvx_begin
        hypograph variable z( 1, s(2) )
        x >= repmat( z, s(1), 1 ); %#ok
        cvx_setnneg( cvx_fastref( z, nn ) );
        cvx_setnpos( cvx_fastref( z, np ) );
    cvx_end
    x = cvx_optval;
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
