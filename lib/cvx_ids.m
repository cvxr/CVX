function y = cvx_ids( varargin )
global cvx___
if cvx___.mversion < 7.1,
    y = zeros( 1, nargin );
    for k = 1 : nargin,
        y(k) = cvx_id( varargin{k} );
    end
else
    y = cellfun( @cvx_id, varargin );
end
