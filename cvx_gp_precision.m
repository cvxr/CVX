function sout = cvx_gp_precision( flag )
if nargin > 0,
    if ~isa( flag, 'double' ) | length( flag ) ~= 1 | ~isreal( flag ) | flag < 0 | flag > 1,
        error( 'Argument must be a scalar between 0 and 1.' );
    end
end
global cvx___
if isempty( cvx___ ), 
    cvx_setpath( 1 ); 
end
cvx_problem = evalin( 'caller', 'cvx_problem', '[]' );
if isa( cvx_problem, 'cvxprob' ),
    s = cvx_problem.gptol;
    if nargin > 0,
        cvx___.problems(index(cvx_problem)).gptol = flag;
    end
else
    s = cvx___.gptol;
    if nargin > 0,
        cvx___.gptol = flag;
    end
end
if nargin == 0 | nargout > 0,
    sout = s;
end
