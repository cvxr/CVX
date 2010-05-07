%%*************************************************************************
%% parameters.m: set OPTIONS structure to specify default
%%               parameters for sqlp.m
%%
%% OPTIONS.vers        : version of direction to use.  
%%                       1 for HKM direction
%%                       2 for NT direction
%%                       0 for the default (which uses the HKM direction if 
%%                       semidefinite blocks exist; and NT direction of SOCP problms)
%% OPTIONS.gam         : step-length parameter,
%% OPTIONS.predcorr    : whether to use Mehrotra predictor-corrector. 
%% OPTIONS.expon       : exponent in decrease of centering parameter sigma. 
%% OPTIONS.gaptol      : tolerance for duality gap as a fraction of the 
%%                       value of the objective functions. 
%% OPTIONS.inftol      : tolerance for stopping due to suspicion of 
%%                       infeasibility.
%% OPTIONS.steptol     : toloerance for stopping due to small steps.
%% OPTIONS.maxit       : maximum number of iteration allowed 
%% OPTIONS.printlevel  : 3, if want to display result in each iteration, 
%%                       2, if want to display only summary,
%%                       1, if want to display warning message,
%%                       0, no display at all.  
%% OPTIONS.stoplevel   : 2, if want to automatically detect termination; 
%%                       1, if want to automatically detect termination, but
%%                          restart automatically with a new iterate
%%                          when the algorithm stagnants because of tiny step-lengths. 
%%                       0, if want the algorithm to continue forever except for 
%%                          successful completion, maximum number of iterations, or 
%%                          numerical failures. Note, do not change this field unless  
%%                          you very sure. 
%% OPTIONS.scale_data  : 1, if want to scale the data before solving the problem, 
%%                          else = 0
%% OPTIONS.rmdepconstr : 1, if want to remove nearly dependent constraints,
%%                          else = 0. 
%% OPTIONS.smallblkdim : block-size threshold determining what method to compute the 
%%                       schur complement matrix corresponding to semidefintie block.
%%                       NOTE: this number should be small, say less than 20. 
%% OPTIONS.parbarrier  : parameter values of the log-barrier terms in the SQLP problem. 
%%                       Default = [],  meaning that the parameter values are all 0.
%% OPTIONS.schurfun    : [], if no user supplied routine to compute the Schur matrix, 
%%                       else, it is a cell array where each cell is either [],
%%                       or contains a string that is the file name where the Schur matrix
%%                       of the associated block data is computed. 
%%                       For example, if the SQLP data has the block structure
%%                       blk{1,1} = '1'; blk{1,2} = 10;
%%                       blk{2,1} = 's'; blk{2,2} = 50;
%%                       and 
%%                       OPTIONS.schurfun{1} = []; 
%%                       OPTIONS.schurfun{2} = 'myownschur', where
%%                             'myownschur' is a function with the calling sequence: 
%%                              function   schur = myownschur(X2,Z2inv,schurfun_par(2,:)); 
%%                       This means that for the first block, the Schur
%%                       matrix is computed by the default method in SDPT3,
%%                       and for the second block, the user supplies the
%%                       routine to compute the associated Schur matrix. 
%% OPTIONS.schurfun_par: [], if no user supplied routine to compute the Schur matrix, 
%%                       else, it is a cell array where the p-th row is either [],  
%%                       or is a cell array containing the parameters needed in 
%%                       the user supplied Schur routine OPTIONS.schurfun{p}.  
%%                       For example, for the block structure described
%%                       above, we may have: 
%%                       OPTIONS.schurfun_par{1} = [];
%%                       OPTIONS.schurfun_par{2,1} = par1;
%%                       OPTIONS.schurfun_par{2,2} = par2;
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%*************************************************************************

   function OPTIONS = sqlparameters; 

   OPTIONS.vers           = 0; 
   OPTIONS.gam            = 0;  
   OPTIONS.predcorr       = 1;
   OPTIONS.expon          = 1; 
   OPTIONS.gaptol         = 1e-8;
   OPTIONS.inftol         = 1e-8; 
   OPTIONS.steptol        = 1e-6; 
   OPTIONS.maxit          = 100; 
   OPTIONS.printlevel     = 3; 
   OPTIONS.stoplevel      = 1; %% do not change this field unless you very sure. 
   OPTIONS.scale_data     = 0; 
   OPTIONS.spdensity      = 0.4; 
   OPTIONS.rmdepconstr    = 0; 
   OPTIONS.smallblkdim    = 50;
   OPTIONS.parbarrier     = []; 
   OPTIONS.schurfun       = [];
   OPTIONS.schurfun_par   = [];
%%*************************************************************************
