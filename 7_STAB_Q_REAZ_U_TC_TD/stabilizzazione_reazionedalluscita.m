%t.d. x(k+1)=A(alpha)x(k)+Bu(k)
%       y(k)=Cx(k)
%A(rho)=[0.25 1 0;0 0.1 0;0 0 0.6+rho];
B=[1;0;1];
C=[1 0 2];
% rho \in [0 0.5]


A1=[0.25 1 0;0 0.1 0;0 0 0.6];
A2=[0.25 1 0;0 0.1 0;0 0 0.6+0.5];

Anom=(A1+A2)/2;

n=3;
m=1;
q=1;
%STEP 1: progetto L : (A(alpha)+L*C) quadr. stabile per ogni alpha
S=sdpvar(n);
Z=sdpvar(n,q,'full');

F=[S>=0];
F=[F,[-S S*A1-Z*C;A1'*S-C'*Z' -S]<=0];
F=[F, [-S S*A2-Z*C;A2'*S-C'*Z' -S]<=0];
diagnostic=optimize(F)
check(F)

S=value(S);
Z=value(Z);
L=inv(S)*Z;
eig(A1-L*C)
eig(A2-L*C)

alpha1=0.3;
alpha2=0.7;
Aalpha=(alpha1*A1+alpha2*A2);
abs(eig(Aalpha-L*C))
%
%STEP 2: progetto di K
Ahat_1=[Anom L*C;A1-Anom A1-L*C];
Ahat_2=[Anom L*C;A2-Anom A2-L*C];
Bhat=[B;zeros(n,m)];


clear F
Q1=sdpvar(n);
Q2=sdpvar(n);
Y1=sdpvar(m,n,'full');
Q=blkdiag(Q1,Q2);
Y=[Y1 zeros(m,n)];

F=[Q1>=0];
F=[F,Q2>=0];
F=[F,[-Q Q*Ahat_1'-Y'*Bhat'; Ahat_1*Q-Bhat*Y -Q]<=0];
F=[F,[-Q Q*Ahat_2'-Y'*Bhat'; Ahat_2*Q-Bhat*Y -Q]<=0];
diagnostic=optimize(F)
check(F)
Q=value(Q);
Q1=value(Q1);
Q2=value(Q2);
Y=value(Y);
Y1=value(Y1);

K=Y1*inv(Q1)
Khat=Y*inv(Q)

alpha1=0.6;
alpha2=0.4;
Aalpha=alpha1*A1+alpha2*A2;

Af=[Anom-B*K L*C;Aalpha-Anom Aalpha-L*C];
abs(eig(Af))



alpha1=0.9;
alpha2=0.1;
Aalpha=alpha1*A1+alpha2*A2;

Af=[Anom-B*K L*C;Aalpha-Anom Aalpha-L*C];
abs(eig(Af))


alpha1=0.3;
alpha2=0.7;
Aalpha=alpha1*A1+alpha2*A2;

Af=[Anom-B*K L*C;Aalpha-Anom Aalpha-L*C];
abs(eig(Af))

