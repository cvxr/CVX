% binary hypothesis testing example

P = [0.70  0.10
     0.20  0.10
     0.05  0.70
     0.05  0.10]
[n,m] = size(P);


% tradeoff curve  1-t1'*q1 vs 1-t2'*q2

nopts = 1000;
weights = logspace(-5,5,nopts);
obj = [0;1];
inds = ones(n,1);

% minimize  -t1'*q1 - w*t2'*q2
% s.t.      t1+t2 = 1,  t1,t2 \geq 0

next =2;
for i=1:nopts
   PW = P*diag([1;weights(i)]);
   [maxvals,maxinds] = max(PW');  % max elt in each row 
   if (~isequal(maxinds', inds(:,next-1)))
       inds(:,next) = maxinds';
       T = zeros(m,n);
       for j=1:n
          T(maxinds(1,j),j) = 1;
       end;
       obj(:,next) = 1-diag(T*P);
       next = next+1;
   end;
end;
plot(obj(1,:), obj(2,:),[0 1], [0 1],'--');
grid on
for i=1:size(obj,2)
   text(obj(1,i),obj(2,i),['a', num2str(i)]);
end;
xlabel('x');
ylabel('y');


% min probability of error 
%   max.  w
%   s.t. 1-t'*p1 >= w
%        t'*p2 >= w
%        0 <= t <= 1 

A = [P(:,1)' 1; -P(:,2)' 1; -eye(n) zeros(n,1); eye(n) zeros(n,1)];
b = [1; 0; zeros(n,1); ones(n,1)];
c = [zeros(n,1);-1];
x = linprog(c,A,b);
%[z,x,s] = mpc(A',-c,b);
Tmp = [1-x(1:n) x(1:n)]'
objmp = 1-diag(Tmp*P);
%plot(objmp(1), objmp(2),'o');
text(objmp(1), objmp(2),'b');
xlabel('x'); ylabel('y');
print -deps roc.eps
