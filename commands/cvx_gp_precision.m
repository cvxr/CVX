function sout = cvx_gp_precision( flag )

%CVX_GP_PRECISION    Changes the precision of CVX's GP approximation.
%   CVX cannot solve geometric programs exactly, so it constructs semidefinite
%   approximations to them. Specifically, the LOG(SUM(EXP(X))) function is
%   approximated by a function LOGSUMEXP_SDP(X). CVX_GP_PRECISION( TOL )
%   instructs CVX to use an absolute precision of TOL; that is, to insure
%       LOGSUMEXP(X) <= LOGSUMEXP_SDP(X) <= LOGSUMEXP(X) + TOL.
%
%   If CVX_GP_PRECISION(TOL) is called within a model---that is, between the
%   statements CVX_BEGIN and CVX_END---then the new precision applies only to 
%   that particular model. If called outside of a model, then the change applies
%   to all subsequent models; that is, it modifies the default tolerance.
%
%   On exit, CVX_GP_PRECISION(TOL) returns the *previous* precision, so that it
%   can be saved and restored later; for example:
%       otol = cvx_gp_precision(tol);
%       cvx_begin
%           ...
%       cvx_end
%       cvx_gp_precision(otol);
%   Of course, this is equivalent to
%       cvx_begin
%           cvx_gp_precision(tol);
%           ...
%       cvx_end
%   but the former syntax it may come in handy if you wish to solve several 
%   models in a row with a different precision.
%
%   CVX_GP_PRECISION, with no arguments, returns the current precision value.
%
%   For more information, see LOGSUMEXP_SDP.

if nargin > 0,
    if ~isa( flag, 'double' ) | length( flag ) ~= 1 | ~isreal( flag ) | flag <= 0 | flag > 1,
        error( 'Argument must be a scalar between 0 and 1, exclusive.' );
    end
end
cvx_global
if isempty( cvx___.problems ),
    s = cvx___.gptol;
    if nargin > 0,
        cvx___.gptol = flag;
    end
else
    s = cvx___.problems(end).gptol;
    if nargin > 0,
        if ~isequal( s, flag ) & ~isa( evalin( 'caller', 'cvx_problem', '[]' ), 'cvxprob' ),
            warning( 'The global CVX GP precision setting cannot be changed while a model is being constructed.' );
        else
            cvx___.problems(end).gptol = flag;
        end
    end
end
if nargin == 0 | nargout > 0,
    sout = s;
end

% Copyright 2008 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
