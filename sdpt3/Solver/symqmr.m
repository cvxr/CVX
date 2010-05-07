%%*************************************************************************
%% symqmr: symmetric QMR with left (symmetric) preconditioner. 
%%         The preconditioner used is based on the analytical
%%         expression of inv(A).  
%%
%% [x,resnrm,solve_ok] = symqmr(A,b,L,tol,maxit) 
%%
%% child function: linsysolvefun.m 
%%
%% A = [mat11 mat12; mat12' mat22].
%% b = rhs vector.
%% if matfct_options = 'chol' or 'spchol' 
%%    L = Cholesky factorization of (1,1) block. 
%%    M = Cholesky factorization of 
%%        Schur complement of A ( = mat12'*inv(mat11)*mat12-mat22).
%% else
%%    L = triangular factors of A.
%%    M = not relevant.
%% end
%% resnrm = norm of qmr-generated residual vector b-Ax. 
%%
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%*************************************************************************

   function  [xx,resnrm,solve_ok] = symqmr(A,b,L,tol,maxit,printlevel) 

   N = length(b); 
   if (nargin < 6); printlevel = 1; end
   if (nargin < 5) | isempty(maxit); maxit = max(30,length(A.mat22)); end;
   if (nargin < 4) | isempty(tol); tol = 1e-10; end; 
   tolb = min(1e-4,tol*norm(b));

   solve_ok = 1; 
   x = zeros(N,1);
   if (norm(x))
      if isstruct(A); Aq = matvec(A,x); else; Aq=A*x; end;
      r = b-Aq;  
   else
      r = b; 
   end
   err = norm(r); resnrm(1) = err; minres = err; xx = x; 
   if (err < 1e-3*tolb); return; end         

   q = precond(A,L,r); 
   tau_old   = norm(q);      
   rho_old   = r'*q; 
   theta_old = 0; 
   d = zeros(N,1); 
   res = r; Ad = zeros(N,1);
%%      
%% main loop
%%
   tiny = 1e-30; 
   for iter = 1:maxit 

       if isstruct(A); Aq = matvec(A,q); else; Aq=A*q; end;     
       sigma = q'*Aq; 
       if (abs(sigma) < tiny)
          solve_ok = 2; 
          if (printlevel); fprintf('*'); end;
          break;
       else
          alpha = rho_old/sigma; 
          r = r - alpha*Aq;
       end
       u = precond(A,L,r); 

       theta = norm(u)/tau_old; c = 1/sqrt(1+theta^2); 
       tau = tau_old*theta*c;
       gam = (c^2*theta_old^2); eta = (c^2*alpha); 
       d = gam*d + eta*q;
       x = x + d; 
%%
       Ad = gam*Ad + eta*Aq;
       res = res - Ad; 
       err = norm(res); resnrm(iter+1) = err; 
       if (err < minres); xx = x; minres = err; end
       if (err < tolb); break; end        
       if (iter > 10) 
          if (err > 0.98*mean(resnrm(iter-10:iter)))
             solve_ok = 0.5; break; 
          end
       end
%% 
       if (abs(rho_old) < tiny)
          solve_ok = 2; 
          if (printlevel); fprintf('*'); end;
          break;
       else
          rho  = r'*u; 
          beta = rho/rho_old; 
          q = u + beta*q; 
       end
       rho_old = rho; 
       tau_old = tau; 
       theta_old = theta; 
   end
   if (iter == maxit); solve_ok = 0.3; end; 
%%
%%*************************************************************************
%% precond: 
%%*************************************************************************

   function Mx = precond(A,L,x)

   m = length(L.perm); m2 = length(x)-m;
   Mx = zeros(length(x),1); 

   for iter = 1:1
      if norm(Mx); r = full(x - matvec(A,Mx)); else; r = full(x); end
      r1 = r(1:m); 
      if (m2 > 0)
         r2 = r(m+[1:m2]);
         w = linsysolvefun(L,r1); 
         z = mexMatvec(A.mat12,w,1) - r2;
         z = L.Mu \ (L.Ml \ (L.Mp*z));
         r1 = r1 - mexMatvec(A.mat12,z); 
      end
      d = linsysolvefun(L,r1);  
      if (m2 > 0)
         d = [d; z];
      end
      Mx = Mx + d;
   end
%%*************************************************************************
%% matvec: matrix-vector multiply.
%% matrix = [A.mat11 A.mat12; A.mat12' A.mat22]
%%*************************************************************************

   function Ax = matvec(A,x);

   m = length(A.mat11); m2 = length(x)-m; 

   if issparse(x); x = full(x); end
   if (m2 > 0)
      x1 = x(1:m); 
   else
      x1 = x; 
   end
   Ax = mexMatvec(A.mat11,x1);
   if (m2 > 0)
      x2 = x(m+[1:m2]);
      Ax = Ax + mexMatvec(A.mat12,x2); 
      Ax2 = mexMatvec(A.mat12,x1,1) + mexMatvec(A.mat22,x2);
      Ax = [full(Ax); full(Ax2)];  
   end
%%*************************************************************************
