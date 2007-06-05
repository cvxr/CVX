%%
%% generate mex files in MA47
%%

   function InstallMA47

   fprintf ('\n Now compiling the mexFunctions in MA47:\n'); 

   fname{1} = 'mexMA47syb.c  depMA47LLB.c  MA47mod.c';
   fname{2} = 'mexMA47fct.c  depMA47LLB.c  MA47mod.c';
   fname{3} = 'mexMA47slv.c  depMA47LLB.c  MA47mod.c';  

   for k = 1:length(fname)
       cmd(['mex  -O ',fname{k}]);  
   end   

%%***********************************************
   function cmd(s) 
   
   fprintf(' %s\n',s); 
   eval(s); 
%%***********************************************


