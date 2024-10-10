%xdot(t)=A(p)*x(t)+B(p)*u(t)
%y(t)=C*x(t)  C quadrata e invertibile x(t)=inv(C)*y(t)  
%nell'ipotesi di stato accessibile   u(t)=K*x(t): (A(alfa)+B(alfa)*K) sia
%quadraticamente stabile sul dominio di incertezza

%A=[p 1.2;0.4 1];
%B=[p; 1];
%C=[1 1.2;0 1];

%p \in [1 1.5]
n=2;  %dim x
m=1;  %dim u
q=2;  %dim y

%progettare u(t) tale che il sistema a ciclo chiuso sia quadraticamente
%stabile


 %(Ai,Bi) i=1,2
A1=[1 1.2;0.4 1];
B1=[1;1];

A2=[1.5 1.2;0.4 1];
B2=[1.5;1];

%V(x(t))=x^T(t)*Px(t)  P=P^T>0
%Q=inv(P)   KQ=Y   
%  Q,Y     Ai*Q+Bi*Y+Q*Ai'+Y'*Bi'<0  i=1,..,N, Q=Q^T>0
Q=sdpvar(n); 
Y=sdpvar(m,n,'full');  %u(t) m componenti
F=[Q>=0];
F=[F,A1*Q+B1*Y+Q*A1'+Y'*B1'<=0];
F=[F,A2*Q+B2*Y+Q*A2'+Y'*B2'<=0];
sol=optimize(F)
check(F)

Q=value(Q);
Y=value(Y);
K=Y*inv(Q)
alfa1=rand(1);
alfa2=1-alfa1;

Aalfa=alfa1*A1+alfa2*A2;
Balfa=alfa1*B1+alfa2*B2;

eig(Aalfa+Balfa*K)  %Re[autovalore]<0


%progettare u(t) tale che il sistema a ciclo chiuso sia D stabile
%Re[z]<=-1=-alpha  D={z \in C: f_D(z)=z+zconiugato+2*alpha<0}  
%D regione LMI  caratterizzata da 2 matrici: L=L^T =2*alpha, M=1

A1=[1 1.2;0.4 1];
B1=[1;1];

A2=[1.5 1.2;0.4 1];
B2=[1.5;1];
alpha=1;
Q=sdpvar(n);
Y=sdpvar(m,n,'full');
F=[Q>=0];
F=[F,2*alpha+A1*Q+B1*Y+Q*A1'+Y'*B1'*Q<=0];
F=[F,2*alpha+A2*Q+B2*Y+Q*A2'+Y'*B2'+Q<=0];
sol=optimize(F)
check(F)

Q=value(Q);
Y=value(Y);
K=Y*inv(Q)
alfa1=rand(1);
alfa2=1-alfa1;

Aalfa=alfa1*A1+alfa2*A2;
Balfa=alfa1*B1+alfa2*B2;

eig(Aalfa+Balfa*K)
