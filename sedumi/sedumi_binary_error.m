function sedumi_binary_error()
%Throw an error indicating the sedumi binaries are not available

err = sprintf('The SeDuMi binaries are not installed.\n%s\n%s',...
              'In Matlab, launch "install_sedumi" in the folder you put the SeDuMi files.',...
              'For more information see the file Install.txt.');

throwAsCaller(MException('SeDuMi:NoBinaries',err));