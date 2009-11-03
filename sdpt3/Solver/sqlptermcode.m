%%*************************************************************************
%% sqlptermcode.m: explains the termination code in sqlp.m
%%               
%%
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%*************************************************************************

   function termcode = sqlptermcode; 

   fprintf('\n   3: norm(X) or norm(Z) diverging'); 
   fprintf('\n   2: dual   problem is suspected to be infeasible') 
   fprintf('\n   1: primal problem is suspected to be infeasible') 
   fprintf('\n   0: max(relative gap,infeasibility) < gaptol'); 
   fprintf('\n  -1: relative gap < infeasibility'); 
   fprintf('\n  -2: lack of progress in predictor or corrector'); 
   fprintf('\n  -3: X or Z not positive definite'); 
   fprintf('\n  -4: difficulty in computing predictor or corrector direction');
   fprintf('\n  -5: progress in relative gap or infeasibility is bad'); 
   fprintf('\n  -6: maximum number of iterations reached'); 
   fprintf('\n  -7: primal infeasibility has deteriorated too much'); 
   fprintf('\n  -8: progress in relative gap has deteriorated'); 
   fprintf('\n  -9: lack of progress in infeasibility'); 
   fprintf('\n')
%%*************************************************************************
