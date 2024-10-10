A=[0 1 0;-1 0 2;0 1 -2];
B=[1 0;0 1; 0 0];
C=[1 0 0;0 1 0];
P=[1;0;0];
Q=[0;1];
n=3; %dim stato
m=2; %num ingressi
p=2; %num uscite controllate
d=1; %dim disturbo
%w(t)=1.5   yd(t)=[3sen(0.5t); t]

%costruzione esosistema
S=[0 0 0 0 0;0 0 -0.5 0 0;0 0.5 0 0 0;0 0 0 0 1;0 0 0 0 0];
r=5;
sys_w=ss(S,zeros(r,1),eye(r),zeros(r,1));
t=0:0.01:50;
w0=[1.5;3;0;0;1];
figure;
initial(sys_w,w0,t);
wtilde=initial(sys_w,w0,t);

%definizione processo soggetto al disturbo generalizzato
Ptilde=[P zeros(n,1) zeros(n,1) zeros(n,1) zeros(n,1)];
Qtilde=[-Q zeros(p,1) eye(p) zeros(p,1)];
Ctilde=-C;
%%%%%%%%%
%risolvo il problema nel caso di IC
%H1
eig(S) %S antistabile
%H2( A,B) stabilizzabile
Mr=ctrb(A,B);
rho=rank(Mr) % rho=n (A,B) raggiungibile
%progetto di K   : A+B*K
K=place(A,B,[-1 -2 -3]);
K=-K;
eig(A+B*K)
%progetto L=Gamma-K*Pigreco
%LME equzioni del regolatore nel caso di IC
Rtilde=[Ptilde;Qtilde];
W=[A B;Ctilde zeros(p,m)];
J=[eye(n) zeros(n,m); zeros(p,n) zeros(p,m)];
%Yalmip
X=sdpvar(n+m,r,'full');
F=[W*X+Rtilde-J*X*S == 0];
diagnostic=optimize(F);
chek(F);
X=value(X);
Pigreco=X(1:n,:);
Gamma=X(n+1:n+m,:);
L=Gamma-K*Pigreco;
%%%%%%%%%
%ipotesi IC rappresentazione del sistema a c.c.
sys_e=ss(A+B*K,B*L+Ptilde,Ctilde,Qtilde);
e=lsim(sys_e,wtilde,t,zeros(n,1));
figure;
subplot(211),plot(t,e(:,1)),grid
subplot(212),plot(t,e(:,2)),grid

%estensione al caso di reazione dall'errore
%progetto dell'osservatore asintotico dello stato
Ae=[A Ptilde;zeros(r,n) S];
Ce=[Ctilde Qtilde];
Mo=obsv(Ae,Ce);
rho=rank(Mo) % rho=n+r (Ae,Ce) osservabile
%
Get=place(Ae',Ce',[-2.5 -3 -3.5 -4 -4.5 -5 -5.5 -6]);
Ge=Get';
eig(Ae-Ge*Ce)

%
G0=Ge(1:n,:);
G1=Ge(n+1:n+r,:);
%costruzione (F,G,H)  del controllore basato sull'osservatore
F=[A+B*K-G0*Ctilde Ptilde-G0*Qtilde+B*L; -G1*Ctilde S-G1*Qtilde];
G=Ge;
H=[K L]; %L=Gamma-K*Pigreco
%rappresentazione in spazio di stao del sistema a ciclo chiuso nel caso di
%reazione dall'errore



Af=[A B*H; G*Ctilde F];
Bf=[Ptilde;G*Qtilde];
Cf=[Ctilde zeros(p,n+r)];
Df=Qtilde;
sys_e=ss(Af,Bf,Cf,Df);
e=lsim(sys_e, wtilde,t,zeros(2*n+r,1));
figure;
subplot(211),plot(t,e(:,1)),grid
subplot(212),plot(t,e(:,2)),grid
%
%rappresentazione sistema a c.c. dove
%l'errore di uscita e' sostituito con l' uscita controllata
Cy=[C zeros(p,n+r)];
Qy=[Q zeros(p,1) zeros(p,1) zeros(p,1) zeros(p,1)];
sys_y=ss(Af,Bf,Cy,Qy);
y=lsim(sys_y, wtilde,t,zeros(2*n+r,1));
figure;
subplot(211),plot(t,y(:,1),'b',t,wtilde(:,3),'r'),grid
subplot(212),plot(t,y(:,2),'b',t,wtilde(:,4),'r'),grid



