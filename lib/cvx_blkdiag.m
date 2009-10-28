function y = cvx_blkdiag( varargin )

% A reimplemenetation of MATLAB's blkdiag function, which is broken in
% certain versions of MATLAB.

if nargin == 1,
    y = varargin{1};
else
    isYsparse = false;
    for k = 1 : nargin,
        x = varargin{k};
        [ p2(k+1), m2(k+1) ] = size(x);
        if issparse(x), isYsparse = true; end
    end
    p1 = cumsum(p2);
    m1 = cumsum(m2);
    if isYsparse
        y = sparse( varargin{1} );
        for k = 2 : nargin,
            y = [y sparse(p1(k),m2(k+1)); sparse(p2(k+1),m1(k)) varargin{k}];
        end
    else
        y = zeros(p1(end),m1(end));
        for k=1:nargin
            y(p1(k)+1:p1(k+1),m1(k)+1:m1(k+1)) = varargin{k};
        end
    end
end
