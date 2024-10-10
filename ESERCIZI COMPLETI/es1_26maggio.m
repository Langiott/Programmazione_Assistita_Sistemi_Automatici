%x(k+1)=A(p,q)*x(k)  
%A=[0.75 0.3;p q] p \in [-0.3 0.3]  q \in [-0.1 0.1]

%verificare se la matrice A e' quadraticamente stabile
n=2;
A1=[0.7 0.3;-0.3 -0.1];
A2=[0.7 0.3;-0.3 0.1];
A3=[0.7 0.3;0.3 -0.1];
A4=[0.7 0.3;0.3 0.1];

%A(alfa)=alfa1*A1+alfa2*A2+alfa3*A3+alfa4*A4;
%alfa1+alfa2+alfa3+alfa4=1  alfa_i>=0, i=1,..,N
%V(x(k))=x^T(k)*P*x(k) P=P^T>0  approccio quadratico P e' unica e
%indipendente da alfa (incertezza)  P=inv(Q)

N=4;
r=1;
Q=sdpvar(n);
F=[Q>=0];
F=[F,[-r*Q A1*Q;Q*A1' -r*Q]<=0];
F=[F,[-r*Q A2*Q;Q*A2' -r*Q]<=0];
F=[F,[-r*Q A3*Q;Q*A3' -r*Q]<=0];
F=[F,[-r*Q A4*Q;Q*A4' -r*Q]<=0];
sol=optimize(F)
check(F)

alfa1=0.2;
alfa2=0.6;
alfa3=0.1;
alfa4=0.1;
Aalfa=alfa1*A1+alfa2*A2+alfa3*A3+alfa4*A4;
abs(eig(Aalfa))

%verificare se la matrice A e'D stabile (r=0.85)
N=4;
r=0.85;
Q=sdpvar(n);
F=[Q>=0];
F=[F,[-r*Q A1*Q;Q*A1' -r*Q]<=0];
F=[F,[-r*Q A2*Q;Q*A2' -r*Q]<=0];
F=[F,[-r*Q A3*Q;Q*A3' -r*Q]<=0];
F=[F,[-r*Q A4*Q;Q*A4' -r*Q]<=0];
sol=optimize(F)
check(F)

alfa1=0.2;
alfa2=0.6;
alfa3=0.1;
alfa4=0.1;
Aalfa=alfa1*A1+alfa2*A2+alfa3*A3+alfa4*A4;
abs(eig(Aalfa))