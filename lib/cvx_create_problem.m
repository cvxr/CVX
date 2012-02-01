function p = cvx_create_problem( varargin )
global cvx___

error( nargchk( 0, 2, nargin ) );
cvx_global
cvx_problem = evalin( 'caller', 'cvx_problem', '[]' );
if isa( cvx_problem, 'cvxprob' ) && cvx___.problems( index( cvx_problem ) ).depth == length( dbstack ) - 2,
    error( 'A cvx problem already exists in this scope.\n(To clear it and start a new one, use the command ''cvx_clear''.)', 1 ); %#ok
end
cvx_problem = cvxprob;
if nargin > 0,
    p = index( cvx_problem );
    for k = 1 : nargin,
        mode = varargin{k};
        if isempty( mode ),
            continue;
        elseif ~ischar( mode ) || size( mode, 1 ) ~= 1,
            cvx_pop( cvx_problem, 'clear' );
            error( 'Arguments must be strings.' );
        end
        switch lower( mode ),
            case 'quiet',
                cvx___.problems( p ).quiet = true;
            case 'set',
                cvx___.problems( p ).complete  = false;
                cvx___.problems( p ).direction = 'find';
            case 'sdp',
                if cvx___.problems( p ).gp,
                    cvx_pop( cvx_problem, 'clear' );
                    error( 'The GP and SDP modifiers cannot be used together.' );
                end
                cvx___.problems( p ).sdp = true;
            case 'gp',
                if cvx___.problems( p ).sdp,
                    cvx_pop( cvx_problem, 'clear' );
                    error( 'The GP and SDP modifiers cannot be used together.' );
                end
                cvx___.problems( p ).gp = true;
            case 'separable',
                cvx___.problems( p ).separable = true;
        end
    end
end
assignin( 'caller', 'cvx_problem', cvx_problem );

% Copyright 2012 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
