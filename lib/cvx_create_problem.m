function p = cvx_create_problem( is_set )
error( nargchk( 0, 1, nargin ) );

global cvx___
cvx_problem = evalin( 'caller', 'cvx_problem', '[]' );
if isa( cvx_problem, 'cvxprob' ) & cvx___.problems( index( cvx_problem ) ).depth == length( dbstack ),
    error( sprintf( 'A cvx problem already exists in this scope.\n(To clear it and start a new one, use the command ''cvx_clear''.)' ) );
end
cvx_problem = cvxprob;
if nargin == 1 & is_set,
    p = index( cvx_problem );
    cvx___.problems( p ).complete  = false;
    cvx___.problems( p ).direction = 'find';
end
assignin( 'caller', 'cvx_problem', cvx_problem );
