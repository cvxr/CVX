function cvx_compile(varargin)

% CVX_RECOMPILE    Attempts to recompile the MEX files that ship with CVX.

global cvx___
mexist = cvx_version('-compile');
line = '---------------------------------------------------------------------------';
fprintf('MEX FILE COMPILATION UTILITY (unsupported)\n');
arg = '';
if mexist,
    if any(strcmp(varargin,'-rebuild')),
        arg = '-rebuild';
    else
        cvx_print({line
            'The required CVX MEX files already exist. To recompile, use the command'
            '    cvx_compile -rebuild'
            line});
        return
    end
end

w = cvx___.where;
s = cvx___.fs;
odir = pwd;
mpath = [w,s,'lib'];
cd(mpath);
files = { 'cvx_bcompress_mex.c', 'cvx_classify_mex.c', 'cvx_eliminate_mex.c' };
for k = 1 : length(files),
    tfile = files{k};
    if ~exist(tfile,'file')
        error('ERROR: Source file %s%s%s not found.', mpath, cvx___.fs, tfile);
    end
    fprintf('Attempting to compile %s:', tfile );
    mex('-O',tfile);
    fprintf(' success.\n' );
end
cvx_print({
    'CVX MEX files successfully compiled.'
    'Now please compile the solver MEX files, if necessary, and run CVX_SETUP:'
    '    cd %s%ssedumi'
    '    install_sedumi -nopath %s'
    '    cd %s%ssdpt3'
    '    install_sdpt3 -nopath %s'
    '    cd %s'
    '    cvx_setup'
    line},w,s,arg,w,s,arg,w);
if length(dbstack) <= 1,
    fprintf('\n');
end
    
function cvx_print(fmt,varargin)
if iscell(fmt), fmt = sprintf('%s\\n',fmt{:}); end
fprintf(fmt,varargin{:});

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
