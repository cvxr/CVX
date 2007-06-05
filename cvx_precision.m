function sout = cvx_precision( flag )
if nargin > 0,
    if isempty( flag ),
        ns = [ sqrt(eps), sqrt(sqrt(eps)) ];
    elseif ischar( flag ),
        if size( flag, 1 ) ~= 1,
            error( 'Invalid precision string.' );
        else
            switch flag,
                case 'default',
                    ns = [ eps^0.5,   eps^0.25 ];
                case 'high',
                    ns = [ eps^0.75,  eps^0.375 ];
                case 'best',
                    ns = [ 0,         eps^0.5   ];
                case 'medium',
                    ns = [ eps^0.375, eps^0.25 ];
                case 'low',
                    ns = [ eps^0.25,  eps^0.25 ];
                otherwise,
                    error( [ 'Invalid precision mode: ', flag ] );
            end
        end
    elseif ~isnumeric( flag ) | numel( flag ) > 2,
        error( 'Argument must be a real number or 2-vector, or a string.' );
    elseif ~isreal( flag ) | any( flag < 0 ) | any( flag >= 1 ),
        error( 'Tolerances must be between 0 (inclusive) and 1 (exclusive).' );
    elseif length( flag ) == 2,
        ns = [ min(flag), max(flag) ];
    elseif flag == 0,
        ns = [ 0, eps^0.5 ];
    else
        ns = [ flag, min(sqrt(flag),max(flag,eps^0.25)) ];
    end
end
global cvx___
if isempty( cvx___ ), 
    cvx_setpath( 1 ); 
end
cvx_problem = evalin( 'caller', 'cvx_problem', '[]' );
if isa( cvx_problem, 'cvxprob' ),
    s = cvx_problem.precision;
    if nargin > 0,
        cvx___.problems(index(cvx_problem)).precision = ns;
    end
else
    s = cvx___.precision;
    if nargin > 0,
        cvx___.precision = ns;
    end
end
if nargin == 0 | nargout > 0,
    sout = s;
end
