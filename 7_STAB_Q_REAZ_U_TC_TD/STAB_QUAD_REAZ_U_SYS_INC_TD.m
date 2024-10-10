   %stabilizzazione con reazione statica dall'uscita T.D.

%A(p)=[1.5 0.5 ;p 0.2] p \in [-1 1];
n=2; %dim stato
B=[1;0];
m=1; %num ingressi
C=[1 0];  
q=1; %num uscita

%A(alpha)=alpha1*A1*alpha2*A2;
%B(alpha)=alpha1*B1*alpha2*B2;
A1=[1.5 0.5;-1 0.2];
B1=B;

A2=[1.5 0.5;1 0.2];
B2=B;

W=sdpvar(n); %W=W'
M=sdpvar(q,q,'full');
N=sdpvar(m,q,'full');

S=[W>=0];
S=[S,[-W W*A1'+C'*N'*B1'; A1*W+B1*N*C -W]<=0];
S=[S,[-W W*A2'+C'*N'*B2'; A2*W+B2*N*C -W]<=0];
S=[S,M*C-C*W==0];

diagnostic=optimize(S)
check(S)

W=value(W); %matrice di Lyapunov
N=value(N);
M=value(M);
K=N*inv(M);  %u=Ky

alpha1=0.2;
alpha2=0.8;

Aalpha=alpha1*A1+alpha2*A2;
Balpha=alpha1*B1+alpha2*B2;

eig(Aalpha+Balpha*K*C)


alpha1=0.9;
alpha2=0.1;

Aalpha=alpha1*A1+alpha2*A2;
Balpha=alpha1*B1+alpha2*B2;

eig(Aalpha+Balpha*K*C)
