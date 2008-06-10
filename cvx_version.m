function [ vv, bb ] = cvx_version
v = '___VERSION___';
b = '___BUILD___';
if nargout == 0,
   disp( sprintf( 'CVX version %s (build %s)', v, b ) );
   if exist('OCTAVE_VERSION'),
       disp( sprintf( 'GNU Octave %s on %s', version, computer ) );
   else
       verd = ver('MATLAB');
       disp( sprintf( 'MATLAB version %s %s on %s', verd.Version, verd.Release, computer ) );
   end
else
   vv = str2num( v );
   bb = str2num( b );
end


