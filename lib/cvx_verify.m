function p = cvx_verify
prob = evalin( 'caller', 'cvx_problem', '[]' );
if ~isa( prob, 'cvxprob' ),
    error( 'CVX:NoProblem', 'No CVX model exists in this scope.' );
end
p = validate( prob );
if ~nargout, clear p; end


