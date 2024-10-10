A=[-1 1 0.5;0 0.2 0; 1 0 -1];
B=[1 0;0 1;0 0];
C=[0 1 0];
P=[1;0;1];
Q=1;
n=3; %dim stato
m=2; %num ingressi
p=1; %num uscite
d=1; %dim disturbo
%w(t)=cos(2t)   %yd(t)=2t+4
%costruzione esosistema 
S=[0 -2 0 0;2 0 0 0;0 0 0 1;0 0 0 0];
r=4; %dim esosistema
sys_w=ss(S,zeros(r,1),eye(r),zeros(r,1));
t=0:0.1:10;
w0=[1;0;4;2];
figure;
initial(sys_w,w0,t)
wtilde=initial(sys_w,w0,t);
% rappresentazione processo soggetto disturbo generalizzato wtilde
Ptilde=[P zeros(n,1) zeros(n,1) zeros(n,1)];
Qtilde=[-Q zeros(p,1) eye(p) zeros(p,1)];
Ctilde=-C;
%verifica H1
eig(S) %ok
%verifica H2 (A,B) stabilizzabile
Mr=ctrb(A,B)
rho=rank(Mr) % rho=n=3 (A,B) raggiungibile
%risolvo il problema come se fossi nel caso di IC
K=place(A,B,[-1 -2 -3]);
K=-K;
eig(A+B*K)
%determino L=Gamma-K*Pigreco risolvendo le equazioni del regolatore
%nel caso di IC
Rtilde=[Ptilde;Qtilde];
W=[A B;Ctilde zeros(p,m)];
J=[eye(n) zeros(n,m);zeros(p,n) zeros(p,m)];
X=sdpvar(n+m,r,'full');
F=[W*X+Rtilde-J*X*S == 0];
diagnostic=optimize(F)
X=value(X);
Pigreco=X(1:n,:);
Gamma=X(n+1:n+m,:);
L=Gamma-K*Pigreco;
%estensione al caso Reazione dall'errore 
%progetto osservatore asintotico dello stato del processo esteso
% [x;wtilde]
%verifica H3 (Ae,Ce) rilevabile
Ae=[A Ptilde;zeros(r,n) S];
Ce=[Ctilde Qtilde];
Mo=obsv(Ae,Ce)
rho=rank(Mo)  %rho=5<n+r=7
Be=[B;zeros(r,m)];
%decomposizione rispetto l'osservabilita'
[a b c T l]=obsvf(Ae,Be,Ce)
Ano=a(1:2,1:2);
eig(Ano) %sottospazio non osservabile e' asint. stabile
%(Ae,Ce) e' rilevabile H3 e' soddisfatta
% alloco autovalori al sottospazio osservabile
Ao=a(3:7,3:7);
Co=c(:,3:7);
Lot=place(Ao',Co',[-2 -2.5 -3 -3.5 -4]);
Lo=Lot';
Gehat=[zeros(2,p);Lo];
Ge=inv(T)*Gehat;
eig(Ae-Ge*Ce)
%determino terna (F,G,H) del controllore dinamico basato sull'osservatore
G=Ge;
G0=G(1:n,:);
G1=G(n+1:n+r,:);
F=[A+B*K-G0*Ctilde Ptilde+B*L-G0*Qtilde;-G1*Ctilde S-G1*Qtilde];
H=[K Gamma-K*Pigreco];
%
%rappresentazione in spazio di stato del sistema a ciclo chiuso
%(Af,Bf,Cf,Df)
Af=[A B*H;G*Ctilde F];
eig(Af)
Bf=[Ptilde;G*Qtilde];
Cf=[Ctilde zeros(p,n+r)];
Df=Qtilde;
sys_e=ss(Af,Bf,Cf,Df);
e=lsim(sys_e,wtilde,t,zeros(2*n+r,1)); %risposta forzata
figure;
plot(t,e),grid
%
Cy=[C zeros(p,n+r)];
Qy=[Q  zeros(p,1) zeros(p,1) zeros(p,1)];
sys_y=ss(Af,Bf,Cy,Qy);
y=lsim(sys_y,wtilde,t,zeros(2*n+r,1)); %risposta forzata

figure;
plot(t,y,'r',t,wtilde(:,3),'b'),grid








