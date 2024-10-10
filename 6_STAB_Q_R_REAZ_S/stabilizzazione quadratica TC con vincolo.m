%------------
%Esercizio 1: stabilizzazione quadratica con reaziona stato
%T.C.
%-----------
%A(p)=[1 0.5;1 p]; p \in [-mu mu]
%A(p) come A(alpha)  alpha \in simplesso di ordine 2
%A(alpha)=alpha1*A1+alpha2*A2;
%B nota ; B=alpha1*B+alpha2*B  alpha1+alpha2=1

mu=0.4;
A1=[1 0.5;1 -mu];
A2=[1 0.5;1 mu];
B=[1;1];
n=2;
m=1;

Q=sdpvar(n);  %P=Q^-1
Y=sdpvar(m,n,'full');  %K=YP
F=[Q>=0];
F=[F,A1*Q+B*Y+Q*A1'+Y'*B'<=0];
F=[F,A2*Q+B*Y+Q*A2'+Y'*B'<=0];
diagnostic=optimize(F);
check(F)
Q=value(Q);
Y=value(Y);
P=inv(Q)  %V=x^TPx
K=Y*P
eig(Q)
eig(A1*Q+B*Y+Q*A1'+Y'*B')
eig(A2*Q+B*Y+Q*A2'+Y'*B')

%alpha \in simplesso
alfa1=0.3;
alfa2=0.7;
Aalfa=alfa1*A1+alfa2*A2;
eig(Aalfa+B*K)

sys_u=ss(Aalfa+B*K,zeros(n,1),-K,zeros(m,1));
x0=[1;1];
t=0:0.1:20;
figure;
initial(sys_u,x0,t)


alfa1=0.5;
alfa2=0.5;
Aalfa=alfa1*A1+alfa2*A2;
eig(Aalfa+B*K)

alfa1=0.7;
alfa2=0.3;
Aalfa=alfa1*A1+alfa2*A2;
eig(Aalfa+B*K)

sys_u_2=ss(Aalfa+B*K,zeros(n,1),-K,zeros(m,1));
figure;
initial(sys_u_2,x0,t)

%--------------------------------------------
%Esercizio 1 con sforzo di controllo ||u||<=2
%--------------------------------------------
Q=sdpvar(n);
Y=sdpvar(m,n,'full');
F=[Q>=0];
F=[F,A1*Q+B*Y+Q*A1'+Y'*B'<=0];
F=[F,A2*Q+B*Y+Q*A2'+Y'*B'<=0];
F=[F, [1 x0';x0 Q]>=0];
F=[F, [Q Y';Y 2^2]>=0];
diagnostic=optimize(F);
check(F)
Q=value(Q);
Y=value(Y);
P=inv(Q);
K=Y*P
eig(Q)

alfa1=0.3;
alfa2=0.7;
Aalfa=alfa1*A1+alfa2*A2;
eig(Aalfa+B*K)
sys_u=ss(Aalfa+B*K,zeros(n,1),-K,zeros(m,1));
x0=[1;1];
t=0:0.1:20;
figure;
initial(sys_u,x0,t)


alfa1=0.6;
alfa2=0.4;
Aalfa=alfa1*A1+alfa2*A2;
eig(Aalfa+B*K)
sys_u=ss(Aalfa+B*K,zeros(n,1),-K,zeros(m,1));
x0=[1;1];
t=0:0.1:20;
figure;
initial(sys_u,x0,t)

%------------%
%Esercizio 2  stabilizzazione %
%------------%

%A(gamma)  B(gamma)  gamma \in [-mu mu]
%A(alpha)=alfa1*A1+alfa2*A2
%B(alpha)=alfa1*B1+alfa2*B2
%unico politopo di vertici (Ai,Bi) i=1,2
mu=0.3;
A1=[1 0.5;1 -mu];
B1=[-mu;1];

A2=[1 0.5;1 mu];
B2=[mu;1];

n=2;
m=1;

Q=sdpvar(n);
Y=sdpvar(m,n,'full');
F=[Q>=0];
F=[F,A1*Q+B1*Y+Q*A1'+Y'*B1'<=0]; %(A1,B1)
F=[F,A2*Q+B2*Y+Q*A2'+Y'*B2'<=0]; %(A2,B2)
diagnostic=optimize(F);
check(F)
Q=value(Q)
Y=value(Y);
P=inv(Q);
K=Y*P
eig(Q)
alfa1=0.3;
alfa2=0.7;

Aalfa=alfa1*A1+alfa2*A2;
Balfa=alfa1*B1+alfa2*B2;
eig(Aalfa+Balfa*K)


alfa1=0.1;
alfa2=0.9;

Aalfa=alfa1*A1+alfa2*A2;
Balfa=alfa1*B1+alfa2*B2;
eig(Aalfa+Balfa*K)

alfa1=0;
alfa2=1;

Aalfa=alfa1*A1+alfa2*A2;
Balfa=alfa1*B1+alfa2*B2;
eig(Aalfa+Balfa*K)






%----------------------------------
%ESERCIZIO 3  (2 parametri incerti)
%-----------------------------------

%A(gamma)  B(rho)  gamma \in [-mu mu], rho \in [-mu mu]
%A(alpha)=alfa1*A1+alfa2*A2 alpha \in simplesso_2
%B(beta)=beta1*B1+beta2*B2  beta \in simplesso_2
mu=0.2;
A1=[1 0.5;1 -mu];
A2=[1 0.5;1 mu];
B1=[-mu;1];
B2=[mu;1];
n=2;
m=1;

Q=sdpvar(n);
Y=sdpvar(m,n,'full');
F=[Q>=0];  %(Ai,Bj) i=1,2  j=1,2   N=2 M=2  N*M
F=[F,A1*Q+B1*Y+Q*A1'+Y'*B1'<=0]; %(A1,B1)
F=[F,A2*Q+B1*Y+Q*A2'+Y'*B1'<=0]; %(A2,B1)
F=[F,A1*Q+B2*Y+Q*A1'+Y'*B2'<=0]; %(A1,B2)
F=[F,A2*Q+B2*Y+Q*A2'+Y'*B2'<=0]; %(A2,B2)
diagnostic=optimize(F);
check(F)
Q=value(Q)
Y=value(Y);
P=inv(Q);
K=Y*P
eig(Q)

alfa1=0.3;
alfa2=0.7;
beta1=0.2;
beta2=0.8;
Aalfa=alfa1*A1+alfa2*A2;
Bbeta=beta1*B1+beta2*B2;
eig(Aalfa+Bbeta*K)



alfa1=0.2;
alfa2=0.8;
beta1=0.7;
beta2=0.3;
Aalfa=alfa1*A1+alfa2*A2;
Bbeta=beta1*B1+beta2*B2;
eig(Aalfa+Bbeta*K)


% sys_u=ss(Aalfa+Bbeta*K,zeros(n,1),-K,zeros(m,1));
% x0=[1;1];
% t=0:0.1:5;
% figure;
% initial(sys_u,x0,t)

alfa1=1;
alfa2=0;
beta1=0;
beta2=1;

Aalfa=alfa1*A1+alfa2*A2;
Bbeta=beta1*B1+beta2*B2;
eig(Aalfa+Bbeta*K)
sys_u=ss(Aalfa+Bbeta*K,zeros(n,1),-K,zeros(m,1));
x0=[1;1];
t=0:0.1:5;
figure;
initial(sys_u,x0,t)

alfa1=0.5;
alfa2=0.5;
beta1=0.2;
beta2=0.8;
Aalfa=alfa1*A1+alfa2*A2;
Bbeta=beta1*B1+beta2*B2;
eig(Aalfa+Bbeta*K)
% sys_u=ss(Aalfa+Bbeta*K,zeros(n,1),-K,zeros(m,1));
% x0=[1;1];
% t=0:0.1:5;
% figure;
% initial(sys_u,x0,t)

%----------------------------------
%ESERCIZIO 3  (2 parametri incerti)
%- vincolo sullo sforzo ||u||<=30

mu=0.2;
A1=[1 0.5;1 -mu];
A2=[1 0.5;1 mu];
B1=[-mu;1];
B2=[mu;1];
n=2;
m=1;
x0=[1;1];
%u=Kx
Q=sdpvar(n);
Y=sdpvar(m,n,'full');
F=[Q>=0];
F=[F,A1*Q+B1*Y+Q*A1'+Y'*B1'<=0];
F=[F,A2*Q+B1*Y+Q*A2'+Y'*B1'<=0];
F=[F,A1*Q+B2*Y+Q*A1'+Y'*B2'<=0];
F=[F,A2*Q+B2*Y+Q*A2'+Y'*B2'<=0];
F=[F, [1 x0';x0 Q]>=0];
F=[F, [Q Y';Y 30^2]>=0];
diagnostic=optimize(F);
check(F)
Q=value(Q)
Y=value(Y);
P=inv(Q);
K=Y*P
eig(Q)

alfa1=0.3;
alfa2=0.7;
beta1=0.2;
beta2=0.8;
Aalfa=alfa1*A1+alfa2*A2;
Bbeta=beta1*B1+beta2*B2;
eig(Aalfa+Bbeta*K)
sys_u=ss(Aalfa+Bbeta*K,zeros(n,1),-K,zeros(m,1));
x0=[1;1];
t=0:0.1:5;
figure;
initial(sys_u,x0,t)

alfa1=1;
alfa2=0;
beta1=0;
beta2=1;

Aalfa=alfa1*A1+alfa2*A2;
Bbeta=beta1*B1+beta2*B2;
eig(Aalfa+Bbeta*K)
sys_u=ss(Aalfa+Bbeta*K,zeros(n,1),-K,zeros(m,1));
x0=[1;1];
t=0:0.1:5;
figure;
initial(sys_u,x0,t)

alfa1=0.5;
alfa2=0.5;
beta1=0.4;
beta2=0.6;

Aalfa=alfa1*A1+alfa2*A2;
Bbeta=beta1*B1+beta2*B2;
eig(Aalfa+Bbeta*K)
sys_u=ss(Aalfa+Bbeta*K,zeros(n,1),-K,zeros(m,1));
x0=[1;1];
t=0:0.1:5;
figure;
initial(sys_u,x0,t)

