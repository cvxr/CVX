function geometric( varargin )
error( nargchk( 1, Inf, nargin ) );

%GEOMETRIC Declare multiple cvx variables.
%   GEOMETRIC VARIABLE x
%   where x is a valid MATLAB variable name, declares a scalar
%   variable for the current cvx problem. A variable with that
%   nm is added to the problem, and a cvxvar object with that
%   nm is created in the current workspace. An error is
%   generated if a cvx problem isn't in the current workspace.
%
%   GEOMETRIC VARIABLE x(n1,n2,...,nk)
%   declares a vector, matrix, or array geometric variable with
%   dimensions n1, n2, ..., nk, each of which must be positive
%   integers.
%
%   GEOMETRIC VARIABLES x1 x2 x3 ...
%   where x1, x2, x3, etc. are valid variable names declares
%   multiple geometric cvx variables. It is exactly equivalent
%   to issuing a separate GEOMETRIC VARIABLE command for each
%   x1, x2, x3, ...
%
%   Structure modifiers such as "symmetric", "toeplitz", etc.
%   are NOT permitted with geometric variables.
%
%   Examples:
%      geometric variable
%      geometric variables x y z;
%      geometric variables x(100) y z(100,10);
%
%   See also VARIABLE, VARIABLES, DUAL, DUALS.

if ~iscellstr( varargin ),
    error( 'GEOMETRIC must be used in command mode.' );
end

mult = true;
mods = { 'geometric_' };
switch varargin{1},
    case 'epigraph',
        if isequal( varargin{2}, 'variable' ) & nargin > 2,
            varargin(1:2) = [];
            mods{end+1} = 'epigraph_';
        else
            error( 'Syntax: geometric epigraph variable <variable>' );
        end
    case 'hypograph',
        if isequal( varargin{2}, 'variable' ) & nargin > 2,
            varargin(1:2) = [];
            mods{end+1} = 'hypograph_';
        else
            error( 'Syntax: geometric hypograph variable <variable>' );
        end
    case 'variable',
        varargin(1) = [];
    case 'variables',
        varargin(1) = [];
        mult = true;
    otherwise,
        error( sprintf( 'Improper use of the ''geometric'' keyword.\n   Type HELP GEOMETRIC for more info.' ) );
end

if mult,
    for k = 1 : length(varargin),
        evalin( 'caller', sprintf( '%s ', 'variable', varargin{k}, mods{:} ) );
    end
else
    evalin( 'caller', sprintf( '%s ', 'variable', varargin{:}, mods{:} ) );
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
