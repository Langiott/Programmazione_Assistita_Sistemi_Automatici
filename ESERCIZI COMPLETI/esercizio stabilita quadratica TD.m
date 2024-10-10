gamma=0.4279; %gamma=0.428 NO

A1=[0.8 -0.25 0 1;1 0 0 0;0 0 0.2 0.03; 0 0 1 0]-gamma*[ 0 0 1 0]'*[0.8 -0.5 0 1];
A2=[0.8 -0.25 0 1;1 0 0 0;0 0 0.2 0.03; 0 0 1 0]+gamma*[ 0 0 1 0]'*[0.8 -0.5 0 1];
n=4;


P=sdpvar(n,n);
F=[P>=0];
F=[F, A1'*P*A1-P<=0];
F=[F, A2'*P*A2-P<=0];
diagnostic=optimize(F)
check(F)

P=value(P);

eig(P)
eig(A1'*P*A1-P)
eig(A2'*P*A2-P)
