function cvx_cleanup( do_warn )
prob = evalin( 'caller', 'cvx_problem', '[]' );
if ~isa( prob, 'cvxprob' ), return; end
global cvx___
try
    np = length( cvx___.problems );
    [ id, p ] = cvx_id( prob );
    if p <= np,
        p = 0;
    else
        pstr = cvx___.problems( p );
        if pstr.id ~= id
            p = 0;
        elseif do_warn && ( p > np || ...
                length(cvx___.classes) == pstr.checkpoint(1) && ...
            length(cvx___.equalities) == pstr.checkpoint(2) && ...
            length(cvx___.cones) == pstr.checkpoint(3) && ...
            isempty(pstr.objective) )
            warning( 'CVX:Empty', 'A non-empty cvx problem already exists in this scope.\n   It is being overwritten.', 1 ); %#ok
        end
    end
catch
    p = 0;
end
evalin( 'caller', sprintf( 'cvx_pop( %d )', p ) );


