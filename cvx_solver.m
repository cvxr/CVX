function sout = cvx_precision( flag )
global cvx___
if isempty( cvx___ ), 
    cvx_setpath( 1 ); 
end
if nargin > 0,
    if isempty( flag ),
        flag = 'SDPT3';
    elseif ~ischar( flag ) | size( flag, 1 ) ~= 1,
        error( 'Argument must be a string.' );
    elseif ~isfield( cvx___.path.solvers, lower( flag )  ),
        error( sprintf( 'Unknown, unusable, or missing solver: %s', flag ) );
    end
end
cvx_problem = evalin( 'caller', 'cvx_problem', '[]' );
if isa( cvx_problem, 'cvxprob' ),
    s = cvx_problem.solver;
    if nargin > 0,
        cvx___.problems(index(cvx_problem)).solver = flag;
    end
else
    s = cvx___.solver;
    if nargin > 0,
        cvx___.solver = flag;
    end
end
if nargin == 0 | nargout > 0,
    sout = s;
end
