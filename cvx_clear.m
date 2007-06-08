function cvx_clear( arg )

% CVX_CLEAR Clears all active cvx data; e.g., in case of an error.

cvx_problem = evalin( 'caller', 'cvx_problem', '[]' );
if ~isa( cvx_problem, 'cvxprob' ),
    global cvx___
    if ~isempty( cvx___ ) & ~isempty( cvx___.problems ),
        ndx = [ cvx___.problems.depth ];
        ndx = min( find( ndx >= length( dbstack ) - 1 ) );
        if isempty( ndx ),
            cvx_problem = [];
        else
            cvx_problem = cvx___.problems( ndx ).self;
            assignin( 'caller', 'cvx_problem', cvx_problem );
        end
    end
end
if ~isempty( cvx_problem ),
    evalin( 'caller', 'pop( cvx_problem, ''clear'' );' );
end
if nargin == 0,
    cvx_clearpath( 1 );
end

% Copyright 2007 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
