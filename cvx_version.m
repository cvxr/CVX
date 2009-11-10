function [ vv, bb ] = cvx_version
v = '___VERSION___';
b = '___BUILD___';
if nargout == 0,
   fprintf( 1, 'CVX version %s (build %s)\n', v, b );
   if exist( 'OCTAVE_VERSION', 'var' ),
       fprintf( 1, 'GNU Octave %s on %s', version, computer );
   else
       verd = ver('MATLAB');
       fprintf( 1, 'MATLAB version %s %s on %s', verd.Version, verd.Release, computer );
   end
else
   vv = str2double( v );
   bb = str2double( b );
end



