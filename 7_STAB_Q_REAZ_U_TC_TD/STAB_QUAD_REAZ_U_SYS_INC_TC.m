%Stabilizzazione con reazione statica dall'uscita (SOF)
%(T.C.)

%A(p)=[1 0.5 ;p -0.2]
%A(alpha)=alpha1*A1+alpha2*A2
%p \in [0 0.5]
n=2; %dim stato
m=1; %num ingressi
q=1; %num uscite
B=[1;0];
C=[1 0];
A1=[1 0.5; 0 -0.2];
B1=B;
A2=[1 0.5;0.5 -0.2];
B2=B;

W=sdpvar(n); %W=W'
M=sdpvar(q,q,'full');
N=sdpvar(m,q,'full');
S=[W>=0]; %LMI
S=[S, A1*W+B1*N*C+W*A1'+C'*N'*B1'<=0]; %LMI
S=[S, A2*W+B2*N*C+W*A2'+C'*N'*B2'<=0]; %LMI
S=[S, M*C-C*W==0]; %LME
diagnostic=optimize(S)
check(S)
W=value(W); %matrice di Lyapunov  V(x)=x'Wx
M=value(M);
N=value(N);
K=N*inv(M);  

alpha1=0.2;
alpha2=0.8;
Aalpha=alpha1*A1+alpha2*A2;
Balpha=alpha1*B1+alpha2*B2;

eig(Aalpha+Balpha*K*C)


alpha1=0.7;
alpha2=0.3;
Aalpha=alpha1*A1+alpha2*A2;
Balpha=alpha1*B1+alpha2*B2;

eig(Aalpha+Balpha*K*C)



