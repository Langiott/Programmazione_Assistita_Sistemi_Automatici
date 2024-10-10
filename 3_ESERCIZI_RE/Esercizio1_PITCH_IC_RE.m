%rappresentazione in spazio di stato (A,B,C) aeroplano

n=4; %stato
m=2; %ingressi
p=1; %uscite
A=[-0.01580 0.02633 -9.810 0;-0.1571 -1.03 0 120.5;0 0 0 1; 0.0005274  -0.01652 0 -1.466];
B=[0.0006056 0;0 -9.496;0 0;0 -5.565];
C=[0 0 1 0];
d=1; %dim. disturbo w=0.1sen(2t);
P=zeros(n,d);
Q=1;

%ydes=30;


%--COSTRUZIONE ESOSISTEMA---%
%S=diag(S1,S2) con  S1=[0 -\omega;0 \omega], wtilde_1(0)=[0.1;0] e S2=0 wtilde_2(0)=30
S=[0 -2 0;2 0 0;0 0 0];
r=3;
sys_w=ss(S,zeros(r,1),eye(r),zeros(r,1));
t=0:0.1:50;
w0=[0.1;0;30];
figure;
initial(sys_w,w0,t);
%wtilde=[0.1cos(2t);0.1sen(2t) ;30]
wtilde=initial(sys_w,w0,t);

%----le matrici (Ptilde Ctilde Qtilde) del processo soggetto al disturbo generalizzato
Ctilde=-C;
Ptilde=[zeros(n,1) P zeros(n,1)]; %Ptilde*wtilde=P*w
Qtilde=[zeros(p,1) -Q eye(p)];  %Qtilde*wtilde=ydes-Q*w

%verifica H1
eig(S) %ok
%verifica H2
Mr=ctrb(A,B);
rank(Mr) %ok

%risolviamo il problema nel caso di informazione completa

%- progetto di K
K=place(A,B,[-1 -1.2 -1.3 -1.4]);
K=-K;
eig(A+B*K)

%-progetto di L
%le equazioni del regolatore date nel teorema IC
Rtilde=[Ptilde;Qtilde];
W=[A B;Ctilde zeros(p,m)];
J=[eye(n),zeros(n,m);zeros(p,n),zeros(p,m)];

X=sdpvar(n+m,r,'full');
F=[W*X+Rtilde-J*X*S == 0];
diagnostic=optimize(F)
check(F)
X=value(X);

 
Pigreco=X(1:n,:);
Gamma=X(n+1:n+m,:);
L=Gamma-K*Pigreco;


sys_e=ss(A+B*K,B*L+Ptilde,Ctilde,Qtilde);
e=lsim(sys_e,wtilde,t,zeros(n,1));
figure;
plot(t,e(:,1)),grid




%simulazione del sistema di controllo a ciclo chiuso la cui uscita e'
%la variabile y(t)
Qy=[zeros(p,1) Q zeros(p,1)];
sys_y=ss(A+B*K,B*(Gamma-K*Pigreco)+Ptilde,C,Qy); %Qy
y=lsim(sys_y,wtilde,t,zeros(n,1));
figure
plot(t,y(:,1),'b',t,wtilde(:,3),'r'), grid, title('pitch')

%setensione al caso di reazione dall'errore di uscita 

%verifica H3
Ce=[Ctilde,Qtilde];
Ae=[A,Ptilde;zeros(r,n),S];
nr=rank(obsv(Ae,Ce))  %nr=n+r=7 OK

Get=place(Ae',Ce',[-1.5 -1.6 -1.7 -1.8 -1.9 -2 -2.1]);
Ge=Get';
G0=Ge(1:n,:);
G1=Ge(n+1:n+r,:);

% La terna (F,G,H) del controllore dinamico con reazione  dall'errore basato
% sull' osservatore
F=[A-G0*Ctilde+B*K,Ptilde-G0*Qtilde+B*L;-G1*Ctilde,S-G1*Qtilde];
G=[G0;G1];
H=[K L];
%
eig([A,B*H;G*Ctilde,F])
%Simulazione del sistema a ciclo chiuso nel caso di reazione dall'errore
sys_e=ss([A,B*H;G*Ctilde,F],[Ptilde;G*Qtilde],[Ctilde,zeros(p,n+r)],Qtilde);
e=lsim(sys_e,wtilde,t,zeros(2*n+r,1));
figure;
plot(t,e(:,1)),grid



%Simulazione del sistema a ciclo chiuso nel caso di reazione dall'errore
%dove in uscita abbiamo la la variabile y(t)
%
%y(t)=Cx+Qw   y(t)=f(x,\xi,wtilde)   
%wtilde=[0.1cos(2t); 0.1sen(2t); 30]
Cy=[C,zeros(p,n+r)]; %Cy*x_cc=C*x
Qy=[zeros(p,1),Q,zeros(p,1)]; %Qy*wtilde=Q*w

sys_y=ss([A,B*H;G*Ctilde,F],[Ptilde;G*Qtilde],Cy,Qy);
y=lsim(sys_y,wtilde,t,zeros(2*n+r,1));

figure
plot(t,y(:,1),'b',t,wtilde(:,3),'r'), grid, title('pitch')

