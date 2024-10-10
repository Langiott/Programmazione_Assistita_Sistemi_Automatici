%T.D.
A1=[0.8 0.3;0.1 -1.5];
A2=[0.8 0.3;0.1 -1];
n=2;
B=[0;1];
m=1; %dim di u
C=[1 1];
q=1; %dim di y
%x(k+1)=A*x(k)+B*u(k)

%SOF  u=Ky=K*C*x(k) 
W=sdpvar(n);  %simmetrica
M=sdpvar(q,q,'full');
N=sdpvar(m,q,'full');

F=[W>=0];
F=[F,[-W (W*A1'+C'*N'*B');A1*W+B*N*C -W]<=0];
F=[F,[-W (W*A2'+C'*N'*B');A2*W+B*N*C -W]<=0];
F=[F,M*C-C*W==0];
sol=optimize(F);
check(F)
W=value(W)
N=value(N)
M=value(M)
K=N*inv(M)

%alpha_i<=1
alfa1=0,3;
alfa2=0.51;
Aalfa=alfa1*A1+alfa2*A2;
abs(eig(Aalfa+B*K*C))



