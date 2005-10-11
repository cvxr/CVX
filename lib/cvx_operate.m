function [ prob, varargout ] = cvx_operate( prob, varargin )

if isempty( prob ),
    pid = 0;
    for k = 1 : nargin - 1,
        var = varargin{k};
        if ~isnumeric( var ),
        prob_n = problem( var );
        pid_n = id( prob_n );
        if pid_n > pid & ~cvx_isconstant( var ),
            pid = pid_n;
            prob = prob_n;
        end
        end
    end
    if pid == 0,
        prob = cvxprob( 'current' );
    end
end
if isa( prob, 'cvxprob' ),
    for k = 1 : nargin - 1,
        varargout{k} = newsubst( prob, varargin{k} );
    end
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
