%%****************************************************
%% mcpschur: compute schur matrix of HKM or NT direction
%%           for MCP and maxG problems.  
%%
%% [schurmat] = mcpschur(X,Zinv,schurfun_par); 
%%
%%   Ak= -ek*ek';
%%
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 21 May 2004
%%****************************************************

  function [schurmat] = mcpschur(X,Zinv,schurfun_pars); 
  
  schurmat = X .* Zinv;
%%****************************************************
