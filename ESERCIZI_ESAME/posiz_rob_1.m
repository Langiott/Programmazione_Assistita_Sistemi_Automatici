%posizionamento robusto nella regione D di stabilita  Re[z]<-4  alfa_1=4

%A(p)=[1 0.6;0.8 p]  B=[0;1]  p \in [-1 1]

%(A(alfa),B(alfa))=\sum_i alfa_i (A_i,B_i)
A1=[1 0.6;0.8 -1];
B1=[0;1];

A2=[1 0.6;0.8 1];
B2=[0;1];

n=2;
m=1;

alfa_1=4;

X=sdpvar(n);
Y=sdpvar(m,n,'full');
F=[X>=0];
F=[F, 2*alfa_1*X+A1*X+B1*Y+X*A1'+Y'*B1'<=0];
F=[F, 2*alfa_1*X+A2*X+B2*Y+X*A2'+Y'*B2'<=0];

diagnostic=optimize(F)
check(F)

X=value(X);
Y=value(Y);

K=Y*inv(X);

alfa1=0.4;
alfa2=0.6;
Aalfa=alfa1*A1+alfa2*A2;
Balfa=alfa1*B1+alfa2*B2;

eig(Aalfa+Balfa*K)