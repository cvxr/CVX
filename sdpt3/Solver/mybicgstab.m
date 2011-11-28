%%*************************************************************************
%% mybicgstab
%%  
%% [xx,resnrm,flag] = mybicgstab(A,b,M1,tol,maxit)
%%
%% iterate on  bb - (M1)*AA*x
%%
%% r = b-A*xtrue; 
%%
%%*************************************************************************

  function [xx,resnrm,flag] = mybicgstab(A,b,M1,tol,maxit,printlevel)

  N = length(b); 
  if (nargin < 6); printlevel = 1; end
  if (nargin < 5) | isempty(maxit); maxit = max(30,length(A.mat22)); end;
  if (nargin < 4) | isempty(tol); tol = 1e-10; end; 
  tolb = min(1e-4,tol*norm(b));
  flag = 1;
 
  x = zeros(N,1); 
  if (norm(x)) 
     if isstruct(A); r = b-matvec(A,x); else; r = b-mexMatvec(A,x); end; 
  else
     r =b; 
  end
  err = norm(r); resnrm(1) = err;  minresnrm = err; xx = x;  
  %%if (err < 1e-3*tolb); return; end

  omega  = 1.0;
  r_tld = r;
%%
%%
%%
  breakyes = 0; 
  smtol = 1e-40; 
  for iter = 1:maxit,           

     rho   = (r_tld'*r);                        
     if (abs(rho) < smtol)
        flag = 2;  
        if (printlevel); fprintf('*'); end;
        breakyes = 1; 
        break; 
     end
     if (iter > 1)
        beta  = (rho/rho_1)* (alp/omega);
        p = r + beta*(p - omega*v);
     else
        p = r;
     end
     p_hat = precond(A,M1,p);
     if isstruct(A); v = matvec(A,p_hat); else; v = mexMatvec(A,p_hat); end;
     alp = rho / (r_tld'*v);
     s = r - alp*v;     
     %%
     s_hat = precond(A,M1,s); 
     if isstruct(A); t = matvec(A,s_hat); else; t = mexMatvec(A,s_hat); end;
     omega = (t'*s) / (t'*t);
     x = x + alp*p_hat + omega*s_hat;              
     r = s - omega*t;
     rho_1 = rho;
     %%
     %% check convergence
     %%
     err = norm(r); resnrm(iter+1) = err;     
     if (err < minresnrm); 
        xx = x; minresnrm = err; 
     end
     if (err < tolb)
        break;  
     end
     if (err > 10*minresnrm) 
        if (printlevel); fprintf('^'); end
        breakyes = 2; 
        break;  
     end       
     if (abs(omega) < smtol)
        flag = 2; 
        if (printlevel); fprintf('*'); end
        breakyes = 1; 
        break; 
     end
  end
  if (~breakyes) & (printlevel >=3); fprintf(' '); end
%%
%%*************************************************************************
%%*************************************************************************
%% precond: 
%%*************************************************************************

   function Mx = precond(A,L,x)

   m = L.matdim; m2 = length(x)-m;
   Mx = zeros(length(x),1); 

   for iter = 1
      if norm(Mx); r = x - matvec(A,Mx); else; r = x; end
      if (m2 > 0)
         r1 = full(r(1:m)); 
      else
         r1 = full(r); 
      end
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
%%*************************************************************************
%% matvec: matrix-vector multiply.
%% matrix = [A.mat11, A.mat12; A.mat12', A.mat22]
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
