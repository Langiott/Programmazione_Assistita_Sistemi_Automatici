
%A(p)=[-p -1;1 0]; p\[1.5 2]
%B= [1,0];
%A(alfa)=alfa1*A1+alfa2*A2
%(Ai,B) vertici i_esimo

B=[1;0];
A1=[-1.5 -1;1 0];
A2=[-2 -1;1 0];

n=2;
m=1;

%settore conico di apertura theta=0.4 rad
% 2*pi : 360° = tetha_rad : tetha_gradi

%tetha_gradi = (tetha*360)/2*pi;

theta_gradi=45;
thet = (theta_gradi*360)/2*pi;
X=sdpvar(n);
Y=sdpvar(m,n,'full');
       
F=[X>=0];
F=[F,[sin(thet)*(A1*X+B*Y+X*A1'+Y'*B') cos(thet)*(A1*X+B*Y-X*A1'-Y'*B');-cos(thet)*(A1*X+B*Y-X*A1'-Y'*B') sin(thet)*(A1*X+B*Y+X*A1'+Y'*B') ]<=0];
F=[F,[sin(thet)*(A2*X+B*Y+X*A2'+Y'*B') cos(thet)*(A2*X+B*Y-X*A2'-Y'*B');-cos(thet)*(A2*X+B*Y-X*A2'-Y'*B') sin(thet)*(A2*X+B*Y+X*A2'+Y'*B') ]<=0];

diagnostic=optimize(F)
check(F)
X=value(X);
Y=value(Y);
K=Y*inv(X);

alpha1=0.3;
alpha2=1-alpha1;
Aalpha=alpha1*A1+alpha2*A2;
eig(Aalpha+B*K)
alpha1=1;
alpha2=1-alpha1;
Aalpha=alpha1*A1+alpha2*A2;
eig(Aalpha+B*K)

alpha1=0.5;
alpha2=1-alpha1;
Aalpha=alpha1*A1+alpha2*A2;
eig(Aalpha+B*K)
%AMMETTE SOLUZIONE!!!!

