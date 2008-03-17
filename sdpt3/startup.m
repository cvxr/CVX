%%*************************************************************************
%% startup.m: set up search paths, and default parameters
%%            for sqlp.m 
%%
%% SDPT3: version 3.0
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last modified: 2 Feb 01
%%*************************************************************************

   warning off;
   
   path(path,strcat(pwd,'/Solver'));
   path(path,strcat(pwd,'/Solver/Mexfun'));
   path(path,strcat(pwd,'/Linsysolver/spchol'));
   path(path,strcat(pwd,'/Linsysolver/MA47'));
   path(path,strcat(pwd,'/Examples'));
   path(path,strcat(pwd,'/HSDSolver'));

   path(path,strcat(pwd,'/testdir'));
%%
%% specify default parameters for sqlp.m,
%% they are specified in the structure called OPTIONS. 
%% 
   OPTIONS = sqlparameters;  
%%*************************************************************************



