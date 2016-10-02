function narginchk(imin,imax)
str = sprintf('error(nargchk(%d,%d,nargin))', imin, imax);
evalin('caller', str);

