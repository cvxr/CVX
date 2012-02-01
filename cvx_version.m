function [ cvx_ver2, cvx_bld ] = cvx_version
cvx_ver = 1.22;
cvx_bld = ':::BUILD:::';
cvx_bdate = ':::BDATE:::';
cvx_dbld = ':::DBUILD:::';
cvx_ddate = ':::DDATE:::';
if nargout == 0,
   fprintf( '\n' );
   fprintf( 'CVX version %g\n', cvx_ver );
   fprintf( '    Code: build %s, %s\n', cvx_bld, cvx_bdate );
   fprintf( '    Documentation: build %s, %s\n', cvx_dbld, cvx_ddate );
   if exist( 'OCTAVE_VERSION', 'var' ),
       fprintf( 'GNU Octave %s on %s\n', version, computer );
       fprintf( 'NOTE: Sorry, Octave support is not yet functional.\n' );
   else
       verd = ver('MATLAB');
       fprintf( 'MATLAB version %s %s on %s\n', verd.Version, verd.Release, computer );
   end
   fprintf( '\n' );
else
    cvx_ver2 = cvx_ver;
end

% Copyright 2012 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.

