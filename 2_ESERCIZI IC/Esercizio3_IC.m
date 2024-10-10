
n=3; m=2; p=2;
A=[-1 0 1;-2 1 1;0 2 1];
B=[1 0;0 0;0 1];
C=[1 0 0;0 1 0];
P=[1 0 0]';
Q=[1 0]';
%w=sen(t) 
d=1; %size di w
%esosistema
S=[0,-1,0,0,0;1,0,0,0,0;0,0,0,0,0;0,0,0,0,1;0,0,0,0,0];
r=5;
sys_w=ss(S,zeros(r,1),eye(r),zeros(r,1));
t=0:0.01:10;
w0=[1,0,1,0,2]';
figure;
initial(sys_w,w0,t);
wtilde=initial(sys_w,w0,t);
%definizione delle matrici del processo soggetto al disturbo generalizzato
Ctilde=-C;
Qtilde=[zeros(p,1),-Q,eye(p),zeros(p,1)];
Ptilde=[zeros(n,1),P,zeros(n,1),zeros(n,1),zeros(n,1)];
%verifica H1 H2  u=K*x+L*wtilde  L=Gamma-K Pigreco
eig(S)

rank(ctrb(A,B))
%progetto di K   tale che (A+B*K) sia asint. stabile
K=place(A,B,[-1 -2 -3]);
K=-K;

%risoluzione equazioni del regolatore
Rtilde=[Ptilde;Qtilde];
W=[A,B;Ctilde,zeros(p,m)];
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
subplot(211),plot(t,e(:,1)),grid
subplot(212),plot(t,e(:,2)),grid
%simulazione del sistema di controllo a ciclo chiuso la cui uscita è
%la variabile y(t)
sys_y=ss(A+B*K,B*L+Ptilde,C,[zeros(p,1),Q,zeros(p,1),zeros(p,1),zeros(p,1)]); %Qy
y=lsim(sys_y,wtilde,t,zeros(n,1));
figure
subplot(211),
plot(t,y(:,1),'b',t,wtilde(:,3),'r'), grid, title('gradino ')
subplot(212),
plot(t,y(:,2),'b',t,wtilde(:,4),'r'), grid, title('rampa 2t')
