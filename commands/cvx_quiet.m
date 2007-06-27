function sout = cvx_quiet( flag )

%CVX_QUIET    CVX output control.
%   CVX_QUIET(TRUE) suppresses all text output from CVX (except for error and
%   warning messages). Specifically, solver progress is not printed.
%
%   CVX_QUIET(FALSE) restores full text output.
%
%   If CVX_QUIET(TF) is called within a model---that is, between the statements
%   CVX_BEGIN and CVX_END---then its effect applies only for the current model.
%   If called outside of a model, the change applies to all subsequent models.
%
%   On exit, CVX_QUIET(TF) returns the *previous* value of the quiet flag, so 
%   that it can be saved and restored later; for example:
%       oflag = cvx_quiet(true);
%       cvx_begin
%           ...
%       cvx_end
%       cvx_quiet(oflag);
%   Of course, this is equivalent to
%       cvx_begin
%       cvx_quiet(true);
%           ...
%       cvx_end
%   but the former syntax is handy if you have a script that solves several 
%   models at once. In general it is good practice to make sure that the
%   CVX_QUIET flag is restored to its previous state upon exit from a script,
%   using either of these techniques.
%
%   CVX_QUIET, with no arguments, returns the current flag value.

if nargin == 1,
    nflag = [];
    if isnumeric(flag) | islogical(flag),
        ns = double(flag) ~= 0;
    elseif ischar(flag) & size(flag,1) == 1,
        switch lower(flag),
            case 'true',
                ns = true;
            case 'false',
                ns = false;
            otherwise,
                error( 'String arugment must be ''true'' or ''false''.' );
        end
    else
        error( 'Argument must be a numeric scalar or a string.' );
    end
end
global cvx___
if isempty( cvx___ ), 
    cvx_setpath( 1 ); 
end
cvx_problem = evalin( 'caller', 'cvx_problem', '[]' );
if isa( cvx_problem, 'cvxprob' ),
    s = cvx_problem.quiet;
    if nargin > 0,
        cvx___.problems(index(cvx_problem)).quiet = ns;
    end
else
    s = cvx___.quiet;
    if nargin > 0,
        cvx___.quiet = ns;
    end
end
if nargin == 0 | nargout > 0,
    sout = s;
end


% Copyright 2007 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.

