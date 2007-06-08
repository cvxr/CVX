function cvx_cd( subdir )

% CVX_CD    Change current working directory to a CVX subdirectory. 
%
% CVX_CD, called with no arguments, changes the current directory to the root
% directory of the CVX package.
%
% CVX_CD <subdir> changes the current directory to the subdirectory <subdir>
% of the CVX package. Useful subdirectories include:
%        doc/          documentation
%        examples/     sample cvx models

oldd = pwd;
cd(cvx_where);
if nargin ~= 0,
	if exist(subdir,'dir'),
	    cd(subdir);
	else
	    cd(oldd);
	    error( sprintf( 'Cannot CD to %s%s (Name is nonexistent or not a directory).', cvx_where, subdir ) );
	end
end

% Copyright 2007 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.

