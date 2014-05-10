function p = validate( prob, varargin )
p = cvx_validate( prob.index_, prob.id_, varargin{:} );

