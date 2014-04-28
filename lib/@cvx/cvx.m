function v = cvx( v, b, d )

%CVX   The CVX disciplined convex programming system.
%   CVX is a modeling framework for building, constructing, and solving
%   disciplined convex programs. CVX has online help for many functions and
%   operations, divided into several subsection:
%      cvx/commands   - Top-level commands to control CVX
%      cvx/keywords   - Keywords for declaring variables and objectives
%      cvx/builtins   - Built-in operators and functions supported in CVX models
%      cvx/functions  - Additional functions added by CVX
%      cvx/sets       - Definitions of common convex sets
%      cvx/structures - Matrix structure definitions
%   CVX also provides an extensive user guide in PDF format, which is found
%   in its top directory. This directory can be found by typing
%      cvx_where
%   at the command prompt.

global cvx___
switch nargin,
    case 0,
        s = [ 0, 1 ];
        b = sparse( 1, 0 );
        d = 0;
    case 1,
        if isa( v, 'cvx' ), 
            return; 
        end
        s = size( v );
        d = prod( s );
        b = sparse( reshape( v, 1, d ) );
    case { 2, 3 },
        switch numel( v ),
            case 2,
                s = v;
            case 1,
                s = [ v, 1 ];
            case 0,
                s = [ 0, 1 ];
            otherwise,
                s = v;
                if s(end) == 1,
                    s = s( 1 : max( 2, find( s > 1, 1, 'last' ) ) );
                end
        end
        if isempty( b ),
            b = sparse( 1, prod( s ) );
        else
            b = sparse( b );
        end
        if nargin < 3,
            d = prod( s );
        end
end
id = cvx___.id + 1; 
cvx___.id = id; 
v = class( struct( 'size_', s, 'basis_', b, 'dof_', d, 'dual_', '', 'id_', id ), 'cvx' );
    
% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
