function x = cvx_fastasgn( x, ndxs, y )
if isa(y,'cvx'), 
    if ~isa(x,'cvx'), x = cvx(x); end
    y = y.basis_;
else
    y = reshape( y, 1, numel(y) );
end
n1 = size(x.basis_,2);
x.basis_(1:size(y,1),ndxs) = y;
n2 = size(x.basis_,2);
if n1 ~= n2, x.size_ = [n2,1]; end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
