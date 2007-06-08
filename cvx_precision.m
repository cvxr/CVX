function sout = cvx_precision( flag )

% CVX_PRECISION    CVX solver precision.
%
% The CVX_PRECISION command controls the precision-related stopping criteria
% for the numerical solver. It defines two precision levels:
%     --- a 'high' precision level, below which a problem is considered
%         accurately solved (returning cvx_status = 'Solved' );
%     --- a 'low' precision level, below which a problem is considered
%         inaccurately solved (returning cvx_status = 'Inaccurate/Solved').
% Problems for which the solver is unable to achieve even the 'low' precision
% level are considered unsolved (returning cvx_status = 'Failed'). These
% tolerance levels apply in appropriate ways to infeasible and unbounded
% problems as well.
%
% CVX_PRECISION(TOL), where TOL is a 2-element vector, sets the low and high
% precisions to MIN(TOL) and MAX(TOL), respectively. If the two values are
% identical, then the low precision level is effectively removed.
%
% CVX_PRECISION(TOL), where TOL is a scalar, sets the low and high precisions
% to min(sqrt(TOL),max(TOL,eps^0.25)) and TOL, respectively. Note that if TOL
% is below eps^0.25, then the low precision level is set equal to the high.
%
% A number of text-based options are provided for convenience:
%     CVX_PRECISION DEFAULT: [eps^0.25,eps^0.5]
%     CVX_PRECISION HIGH   : [eps^0.375,eps^0.75] 
%     CVX_PRECISION MEDIUM : [eps^0.25,eps^0.375]
%     CVX_PRECISION LOW    : [eps^0.25,eps^0.25]
%     CVX_PRECISION BEST   : (see below)
%
% CVX_PRECISION BEST and CVX_PRECISION(0) select a special 'best effort' mode.
% In this case, the solver proceeds as long as it can make progress. If the
% precision achieved is eps^0.5 or better, then the problem is considered
% solved (cvx_status='Solved'); otherwise, the problem is considered unsolved
% (cvx_status='Failed'). This mode can produce higher precision but at cost
% of slower performance.
%
% If CVX_PRECISION(TOL) is called within a model---that is, between the
% statements CVX_BEGIN and CVX_END---then the new precision applies only to that
% particular model. If called outside of a model, then the change applies to
% all subsequent models; that is, it modifies the default tolerance.
%
% On exit, CVX_PRECISION(TOL) returns the *previous* precision, so that it
% can be saved and restored later; for example:
%    otol = cvx_precision(tol);
%    cvx_begin
%       ...
%    cvx_end
%    cvx_precision(otol);
% Of course, this is equivalent to
%    cvx_begin
%        cvx_precision(tol);
%        ...
%    cvx_end
% but the former syntax it may come in handy if you wish to solve several models
% in a row with a different precision.
%
% CVX_PRECISION, with no arguments, returns the current precision value.

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

% Copyright 2007 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
