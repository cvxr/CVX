function sout = cvx_power_warning( flag )

%CVX_POWER_WARNING   Warning message control for power->SOCP conversion.
%   CVX converts power functions like x.^p, for variable x and fixed p, into
%   solvable form using an SOCP transformation. For quadratics x.^2 and square
%   roots x.^(1/2), a single second-order cone is required; for other powers,
%   the number depends on the rational representation of the exponent p.
%
%   CVX_POWER_WARNING(Q) instructs CVX to issue a warning if the resulting
%   transformations requires more than Q second-order cones. The default value
%   is 10, which is not likely to be exceeded for typical choices of P.

if nargin > 0,
    if isempty( flag ),
        ns = 10;
    elseif ~isnumeric( flag ) | ~isreal( flag ) | numel( flag ) > 1 | flag <= 0 | flag ~= floor( flag ),
        error( 'Argument must be a positive integer.' );
    else
        ns = flag;
    end
end
global cvx___
if isempty( cvx___ ), 
    cvx_setpath( 1 ); 
end
cvx_problem = evalin( 'caller', 'cvx_problem', '[]' );
if isa( cvx_problem, 'cvxprob' ),
    s = cvx_problem.rat_growth;
    if nargin > 0,
        cvx___.problems(index(cvx_problem)).rat_growth = ns;
    end
else
    s = cvx___.rat_growth;
    if nargin > 0,
        cvx___.rat_growth = ns;
    end
end
if nargin == 0 | nargout > 0,
    sout = s;
end

% Copyright 2007 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
