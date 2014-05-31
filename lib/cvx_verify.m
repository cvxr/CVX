function p = cvx_verify
prob = evalin( 'caller', 'cvx_problem', '[]' );
if ~isa( prob, 'cvxprob' ),
    cvx_throw( 'No CVX model exists in this scope.' );
end
p = validate( prob );
if ~nargout, clear p; end


