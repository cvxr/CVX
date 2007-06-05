%%
%% generate mex files in spchol
%%
   
   function Installspchol


   fprintf ('\n Now compiling the mexFunctions in spchol:\n'); 

   fname{1} = 'mexordmmd.c ordmmd.c';
   fname{2} = 'mexsymbfct.c  symbfct.c';
   fname{3} = 'choltmpsiz.c';  
   fname{4} = 'cholsplit.c';
   fname{5} = 'mexsparchol.c sparchol2.c sdmauxFill.c sdmauxScalarmul.c';
   fname{6} = 'mexfwblkslv.c sdmauxScalarmul.c';
   fname{7} = 'mexbwblkslv.c sdmauxFill.c sdmauxRdot.c';

   for k = 1:length(fname)
       cmd(['mex  -O ',fname{k}]);  
   end      
   cd .. 
   cd ..
%%***********************************************
   function cmd(s) 
   
   fprintf(' %s\n',s); 
   eval(s); 
%%***********************************************
