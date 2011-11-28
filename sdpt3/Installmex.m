%%
%% Run this script in Matlab command window 
%%

   function Installmex(recompile)

   if (nargin==0); recompile = 0; end

   curdir = pwd;  
   fprintf(' current directory is:  %s\n',curdir);    
%%
%% generate mex files in Mexfun
%% 
   if strcmp(computer,'PCWIN64') | strcmp(computer,'GLNXA64')
      computer_model = 64; 
   else
      computer_model = 32; 
   end
   matlabversion = sscanf(version,'%f');
   matlabversion = matlabversion(1);
   tmp = version('-release'); 
   matlabrelease = str2num(tmp(1:4));
   fsp = filesep;
%%
%%
%%
   src = [curdir,fsp,'Solver',fsp,'Mexfun']; 
   eval(['cd ','Solver',fsp,'Mexfun']); 
   fprintf ('\n Now compiling the mexFunctions in:\n'); 
   fprintf (' %s\n',src);    

   fname{1} = 'mexProd2'; 
   fname{2} = 'mexProd2nz';
   fname{3} = 'mexinprod';
   fname{4} = 'mexmat';
   fname{5} = 'mexsmat';
   fname{6} = 'mexsvec';
   fname{7} = 'mexschur'; 
   fname{8} = 'mexqops';
   fname{9} = 'mexexpand';
   fname{10} = 'mexskron';
   fname{11} = 'mexnnz';
   fname{12} = 'mexschurfun';
   fname{13} = 'mexMatvec';
   fname{14} = 'mextriang';
   fname{15} = 'mextriangsp';
%%
   if (matlabversion < 7.3) & (matlabrelease <= 2008)
      mexcmd = 'mex  -O  '; 
   else
      mexcmd = 'mex  -O -largeArrayDims  '; 
   end
   ext = mexext; 
   for k = 1:length(fname)
      existmex = exist([fname{k},'.',ext]); 
      if (existmex ~= 3) | (recompile)
         cmd([mexcmd,fname{k},'.c']);  
      end
   end 
   cd .. 
   cd ..
%%***********************************************
   function cmd(s) 
   
   fprintf(' %s\n',s); 
   eval(s); 
%%***********************************************
