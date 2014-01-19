function sedumi_binary_error()
%Throw an error indicating the sedumi binaries are not available

if exist('OCTAVE_VERSION','builtin'),
    software = 'Octave';
else
    software = 'Matlab';
end

fs = filesep;
mpath = mfilename('fullpath');
temp = strfind( mpath, fs );
mpath = mpath( 1 : temp(end) - 1 );

error( 'SeDuMi:NoBinaries', ...
    [ 'The SeDuMi binaries for your platform are not installed.\n', ...
      'To rectify, run the following commands from the %s command line:\n\n', ...
      '   cd %s\n', ...
      '   install_sedumi\n\n', ...
      'For more information, read %s%sInstall.txt.' ], ...
    software, mpath, mpath, fs );
