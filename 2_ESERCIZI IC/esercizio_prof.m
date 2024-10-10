A=[-0.01 0;0 -0.02];
B=[1 1;-0.25 0.75];
C=[0.01 0;0 1];
n=2;
m=2; %ingressi
p=2; %uscite

S=[0 0;0 0];
r=2;
sys_w=ss(S,zeros(r,1),eye(r),zeros(r,1));
t=0:0.01:100;
w0=[0.05;1.5];
figure;
initial(sys_w,w0,t);
wtilde=initial(sys_w,w0,t);

Ctilde=-C;

Ptilde=zeros(n,p);
Qtilde=eye(p);
%verifica H1
eig(S)
 %ok
%verifica H2
Mr=ctrb(A,B);
rank(Mr) %ok

%risolviamo il problema nel caso di informazione completa
K=place(A,B,[-0.1 -0.2]);
K=-K;
eig(A+B*K)

%risolvere le equazioni del regolatore date nel teorema IC
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


sys_e=ss(A+B*K,B*(Gamma-K*Pigreco)+Ptilde,Ctilde,Qtilde);
e=lsim(sys_e,wtilde,t,zeros(n,1));
figure;
plot(t,e(:,1)),grid


figure;
plot(t,e(:,2)),grid

%simulazione del sistema di controllo a ciclo chiuso la cui uscita Ë
%la variabile y(t)
sys_y=ss(A+B*K,B*(Gamma-K*Pigreco)+Ptilde,C, zeros(p)); %Qy
y=lsim(sys_y,wtilde,t,zeros(n,1));
figure
plot(t,y(:,1),'b',t,wtilde(:,1),'r'), grid, title('Portata ')
figure
plot(t,y(:,2),'b',t,wtilde(:,2),'r'), grid, title('Concentrazione')

