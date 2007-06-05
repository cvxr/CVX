function sout = cvx_power_warning( flag )
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
