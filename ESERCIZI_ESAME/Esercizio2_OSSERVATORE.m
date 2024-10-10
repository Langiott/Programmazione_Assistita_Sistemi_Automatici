    %T.D. 
%%
A1=[0.8 0.3;0.1 -1.5];
A2=[0.8 0.3;0.1 -1];
n=2;
B=[0;1];
m=1;
C=[1 1];
q=1;
%step 1: progetto L : (A(alpha)+L*C) sia quad. stabile
%(A(alpha)+L*C)
S=sdpvar(n);
Z=sdpvar(n,q,'full');

F=[S>=0];
F=[F,[-S S*A1-Z*C;A1'*S-C'*Z' -S]<=0];
F=[F,[-S S*A2-Z*C;A2'*S-C'*Z' -S]<=0];
sol=optimize(F)
check(F)
S=value(S);
Z=value(Z);
L=inv(S)*Z

alfa1=rand(1);
alfa2=1-alfa1;
Aalfa=alfa1*A1+alfa2*A2;
abs(eig(Aalfa-L*C))
%step 2: fissato L : Khat: Ahat(alfa)-Bhat*Khat sia quad. stabile

Anom=(A1+A2)/2;   
Ahat1=[Anom L*C;A1-Anom A1-L*C];
Ahat2=[Anom L*C;A2-Anom A2-L*C];

Bhat=[B;zeros(n,m)];

%Khat=Y*inv(Q)=[K 0]  Y=[Y1 0]
Q1=sdpvar(n);
Q2=sdpvar(n);
Q=blkdiag(Q1,Q2);

Y1=sdpvar(m,n,'full');
Y=[Y1 zeros(m,n)];

F1=[Q>=0];
F1=[F1,[-Q Q*Ahat1'-Y'*Bhat';Ahat1*Q-Bhat*Y -Q]<=0];
F1=[F1,[-Q Q*Ahat2'-Y'*Bhat';Ahat2*Q-Bhat*Y -Q]<=0];
sol=optimize(F1)
check(F1)
Q=value(Q);
Y=value(Y);
Khat=Y*inv(Q)

eig(Ahat1-Bhat*Khat)
eig(Ahat2-Bhat*Khat)

alfa1=rand(1);
alfa2=1-alfa1;
abs(eig((alfa1*Ahat1+alfa2*Ahat2)-Bhat*Khat))