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
   
   SDPT3Home = pwd;
   addpath(genpath(SDPT3Home)); 
%%
%% obsolete:
%% eval(['addpath ',strcat(SDPT3Home,'/')]);
%% eval(['addpath ',strcat(SDPT3Home,'/Solver')]);
%% eval(['addpath ',strcat(SDPT3Home,'/Solver/Mexfun')]);
%% eval(['addpath ',strcat(SDPT3Home,'/HSDSolver')]);
%% eval(['addpath ',strcat(SDPT3Home,'/Examples')]);
%%
%% obsolete: eval(['addpath ',strcat(SDPT3Home,'/Linsysolver/spchol')]);
%%
%% specify default parameters for sqlp.m,
%% they are specified in the structure called OPTIONS. 
%%
%% OPTIONS = sqlparameters;  
%%*************************************************************************



