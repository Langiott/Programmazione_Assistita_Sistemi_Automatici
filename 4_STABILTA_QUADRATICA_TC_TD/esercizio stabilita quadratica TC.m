
%STABILITA' QUADRATICA
mu=0.7526; %mu=0.753 NO;
A1=[-4 4;-5 0]-mu*[-2 2 ;-1 4];
A2=[-4 4;-5 0]+mu*[-2 2 ;-1 4];

n=2;

P=sdpvar(n,n);
F=[P>=0];
F=[F, A1'*P+P*A1<=0];
F=[F, A2'*P+P*A2<=0];
diagnostic=optimize(F)
check(F)

P=value(P);
eig(P)
eig(A1'*P+P*A1)
eig(A2'*P+P*A2)
