function xi = cvx_invert_structure( x, compact )

%CVX_INVERT_STRUCTURE Compute a right-inverse of a structure mapping.

if nargin == 1,
    
    xi = x'/(x*x');
    [ii,jj,vv] = find(xi);
    [vn,vd] = rat(vv);
    xi = sparse(ii,jj,vn./vd,size(x,2),size(x,1));
    
elseif isreal( x ),
    
    [LL,UU] = lu(x);
    [jj,ii] = find(UU');
    dd = [true;diff(ii)~=0];
    jj = jj(dd);
    [i2,j2,vv] = find( inv(UU(:,jj)) / LL );
    [vn,vd] = rat(vv);
    xi = sparse(jj(i2),j2,vn./vd,size(x,2),size(x,1));
    
else
    
    x = [ real(x), imag(x) ];
    x = x(:,[1:end/2;end/2+1:end]);
    xi = cvx_invert_structure( x, true );
    xi = xi(1:2:end,:) - sqrt(-1) * xi(2:2:end,:);
    
end

% Copyright 2010 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
