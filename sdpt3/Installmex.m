%%
%% Run this script in Matlab command window 
%%

   function Installmex 

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

   if (matlabversion < 7.3) 
      mexcmd = 'mex  -O  '; 
   else
      mexcmd = 'mex  -O -largeArrayDims  '; 
   end
   for k = 1:length(fname)
      cmd([mexcmd,fname{k},'.c']);  
   end 
   cd .. 
   cd ..
%%
%% generate mex files in spchol
%%
   if (matlabversion < 7.3) 
      clear fname
      src = [curdir,fsp,'Linsysolver',fsp,'spchol']; 
      eval(['cd ','Linsysolver',fsp,'spchol']); 
      fprintf ('\n Now compiling the mexFunctions in:\n'); 
      fprintf (' %s\n',src);    

      fname{1} = 'mexordmmd.c ordmmd.c';
      fname{2} = 'mexsymbfct.c  symbfct.c';
      fname{3} = 'choltmpsiz.c';
      fname{4} = 'cholsplit.c';
      fname{5} = 'mexsparchol.c sparchol2.c sdmauxFill.c sdmauxScalarmul.c';
      fname{6} = 'mexfwblkslv.c sdmauxScalarmul.c';
      fname{7} = 'mexbwblkslv.c sdmauxFill.c sdmauxRdot.c';
      mexcmd = 'mex  -O  '; 
      for k = 1:length(fname)
         cmd([mexcmd,fname{k}]);  
      end      
      cd .. 
      cd ..
   end
%%***********************************************
   function cmd(s) 
   
   fprintf(' %s\n',s); 
   eval(s); 
%%***********************************************
