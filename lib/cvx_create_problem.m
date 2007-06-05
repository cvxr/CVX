function p = cvx_create_problem( varargin )
error( nargchk( 0, 2, nargin ) );

global cvx___
cvx_problem = evalin( 'caller', 'cvx_problem', '[]' );
if isa( cvx_problem, 'cvxprob' ) & cvx___.problems( index( cvx_problem ) ).depth == length( dbstack ) - 2,
    error( sprintf( 'A cvx problem already exists in this scope.\n(To clear it and start a new one, use the command ''cvx_clear''.)' ) );
end
cvx_problem = cvxprob;
if nargin > 0,
    p = index( cvx_problem );
    for k = 1 : nargin,
        mode = varargin{k};
        if isempty( mode ),
            continue;
        elseif ~ischar( mode ) | size( mode, 1 ) ~= 1,
            pop( cvx_problem, 'clear' );
            error( 'Arguments must be strings.' );
        end
        switch lower( mode ),
            case 'set',
                cvx___.problems( p ).complete  = false;
                cvx___.problems( p ).direction = 'find';
            case 'sdp',
                if cvx___.problems( p ).gp,
                    pop( cvx_problem, 'clear' );
                    error( 'The GP and SDP modifiers cannot be used together.' );
                end
                cvx___.problems( p ).sdp = true;
            case 'gp',
                if cvx___.problems( p ).sdp,
                    pop( cvx_problem, 'clear' );
                    error( 'The GP and SDP modifiers cannot be used together.' );
                end
                cvx___.problems( p ).gp = true;
            case 'separable',
                cvx___.problems( p ).separable = true;
        end
    end
end
assignin( 'caller', 'cvx_problem', cvx_problem );
