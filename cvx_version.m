function [ cvx_ver, cvx_bld ] = cvx_version
cvx_ver = 1.21;
cvx_bld = ':::BUILD:::';
if nargout == 0,
   fprintf( 1, 'CVX version %g (build %s)\n', cvx_ver, cvx_bld );
   if exist( 'OCTAVE_VERSION', 'var' ),
       fprintf( 1, 'GNU Octave %s on %s', version, computer );
   else
       verd = ver('MATLAB');
       fprintf( 1, 'MATLAB version %s %s on %s', verd.Version, verd.Release, computer );
   end
end
cvx_bld = str2double(cvx_bld);

