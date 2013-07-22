function err=dimacserrors(At,b,c,K,x,y)
%Compute the DIMACS error measures for a given solution in SeDuMi format
%Reload the original coefficients
if ~isfield(K,'f')
    K.f=0;
end
s=c-At*sparse(y);       %To make s sparse
cx=sum(c.*x);        %faster than c'*x
by=sum(b.*y);
xs=sum(x.*s);
normb=norm(b,1);
normc=norm(c,1);
err=zeros(1,6);
%     Error measures.
%     Primal infeasibility
err(1)=norm(x'*At-b',2)/(1+normb);
%     Primal cone infeasibility
err(2)=max(0,-min(eigK(full(x(K.f+1:end)),K))/(1+normb));
%     Dual infeasibility
%info.err(3)=0.0; %s is not maintained explicitely
%     Dual cone infeasibility
err(4)=max(0,-min(eigK(full(s(K.f+1:end)),K))/(1+normc));
%     Relative duality gap
err(5)=(cx-by)/(1+abs(cx)+abs(by));
%     Relative complementarity
err(6)=xs/(1+abs(cx)+abs(by));
