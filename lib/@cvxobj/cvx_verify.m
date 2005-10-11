function s = cvx_verify( varargin )
global cvx___
if nargin == 1,
   x = varargin{1};
   if x.index_ > length( cvx___.problems ) | x.id_ < id( cvx___.problems( x.index_ ).self ), % | cvx___.problems( x.index_ ).locked,
       s = 'Attempt to reference an object from a completed cvx specification.';
       return
   end
else,
    for k = 1 : nargin,
    	s = cvx_verify( varargin{k} );
    	if ~isempty( s ), return; end
    end
end
s = '';

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
