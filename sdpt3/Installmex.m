%%
%% Run this script in Matlab command window 
%%

   function Installmex 
   curdir = pwd;  
   fprintf(' current directory is:  %s\n',curdir);    
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
   
   if exist('OCTAVE_VERSION'),
      mexcmd = 'mex ';
   else
      mexcmd = 'mex  -O -largeArrayDims  '; 
   end
   for k = 1:length(fname)
       cmd = [mexcmd,fname{k},'.c'];
       disp( cmd );
       eval( cmd );
   end 
   cd .. 
   cd ..

