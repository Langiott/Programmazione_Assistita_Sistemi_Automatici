% (T.C) stabilita' robusta P(alpha) mu=1.6666;
%\dot x=A(delta)x+Bu e' robustamente stabilizzabile con reazione dallo
%stato u=Kx nel dominio di incertezza delta \in [-mu mu]?
mu=1.8;
A1=[-4 4;-5 0]-mu*[-2 2;-1 4];
A2=[-4 4;-5 0]+mu*[-2 2;-1 4];
B=[1;0];
n=2;
m=1;
P1=sdpvar(n);
P2=sdpvar(n);
Y=sdpvar(m,n,'full');
G=sdpvar(n,n,'full');
epsilon=1; %gridding (0,10] step 1
S=[P1>=0];
S=[S,P2>=0];
S=[S,[A1*G+B*Y+G'*A1'+Y'*B' P1-G'+epsilon*(A1*G+B*Y); P1-G+epsilon*(G'*A1'+Y'*B') -epsilon*(G+G')]<=0];
S=[S,[A2*G+B*Y+G'*A2'+Y'*B' P2-G'+epsilon*(A2*G+B*Y); P2-G+epsilon*(G'*A2'+Y'*B') -epsilon*(G+G')]<=0];
diagnostic=optimize(S);
check(S)
P1=value(P1)
P2=value(P2)
G=value(G)
Y=value(Y)
K=Y*inv(G)
alfa1=0.3;
alfa2=0.7;
Aalfa=alfa1*A1+alfa2*A2;
eig(Aalfa+B*K)
%%%%%%%%%%%%%%%%%%%%
%(T.D) stabilita' robusta P(alpha) mu=0.4619;
%x(k+1)=A(delta)x(k)+B(delta)u(k) e' robustamente stabilizzabile con reazione dallo
%stato u=Kx nel dominio di incertezza delta \in [-mu mu]?
mu=0.6;
A1=[0.8 -0.25 0 1;1 0 0 0;0 0 0.2 0.03;0 0 1 0]-mu*[0;0;1;0]*[0.8 -0.5 0 1];
B1=[0;0;-mu;0];
A2=[0.8 -0.25 0 1;1 0 0 0;0 0 0.2 0.03;0 0 1 0]+mu*[0;0;1;0]*[0.8 -0.5 0 1];
B2=[0;0;mu;0];
n=4;
m=1;
P1=sdpvar(n);
P2=sdpvar(n);
Y=sdpvar(m,n,'full');
G=sdpvar(n,n,'full');
S=[P1>=0];
S=[S,P2>=0];
S=[S,[-P1 A1*G+B1*Y; G'*A1'+Y'*B1' P1-(G+G')]<=0];
S=[S,[-P2 A2*G+B2*Y; G'*A2'+Y'*B2' P2-(G+G')]<=0];
diagnostic=optimize(S);
check(S)
P1=value(P1)
P2=value(P2)
G=value(G)
Y=value(Y)
K=Y*inv(G)
alfa1=0.3;
alfa2=0.7;
Aalfa=alfa1*A1+alfa2*A2;
Balfa=alfa1*B1+alfa2*B2;
abs(eig(Aalfa+Balfa*K))
%%%%%%%%%%%%%%%%%