
% SCOPO: 
% Progettare una legge di controllo del sistema conoscendo lo 
%stato e il distrurbo , siamo nel caso di INFORMAZIONE COMPLETA 

%SISTEMA:
% Il sistema ha in ingresso il  disturbo w_tilde(t), e in uscita 
% troviamo l'errore generalizzato e(t)
%
%    dx(t) =  A*x(t) + B*u(t) + P_tilde* w_tilde(t)
%    e(t) = C_tilde*x(t) + Q_tilde * w_tilde(t)= y(t)-yd(t)
%
%LEGGENDA:
%    x--->n   e--->p        A--->n,n    
%    u--->m   w_tilde--->r  B--->m,1

%IPOTESI:
%H1: S è antistabile 
%H2: A+BK sia stabilizzabile 

%TEOREMA:
%Dato un sistema a ciclo chiuso e sotto ipotesi H1 e H2 , il problema di 
%regolazione asintotica ha soluzione:
% <---> Esiste GAMMA & PI :
%       | PI*S= A*PI+B*GAMMA+P_tilde
%       | 0= C_tilde*PI+Q_tilde
%E la possibile soluzione del controllore sarà:
%      u = K*x + (GAMMA-K*PI)*W_tilde
%Ho costruito insomma:
% K---> Spettro(A+BK)c C-
% L---> è una combinazione tra PI & GAMMA---> L = GAMMA-KPI

%ECOSISTEMA:
% | D(w_tilde) = S*w_tilde
% | w_tilde(0) = w_tilde_0
% | (P*w) = (P_tilde * w_tilde)
%      w_tilde=[w,yd]

%SEMPLIFICAZIONI: 

%COSTRUISCO MATRICI DEL SISTEMA
A= [-0.01 0;0 -0.02];
n=2;
B=[1 1; -0.3 0.75];
B=-B;
m=2; 
C=[0.01 0;0 1]; %matrice invertibile 
p=2;
P= [0 0;0.015 0.005];
d=2;
Q=zeros(p,d);

%COSTRUZIONE ECOSISTEMA S

%INIZAILIZZAZIONE 
S1=0;
r1=1;
S2=0;
r2=1;
omega=1;
S3=[0 -omega;omega 0];
r3=2;


S=blkdiag(S1,S2,S3);
r=r1+r2+r3;
sys_w=ss(S,zeros(r,1),eye(r),zeros(r,1));

%RAPPRESENTAZIONE GRAFICA 
w0=[0.5;1.5;0.01;0];
t= 0:0.1:10;
figure;
initial(sys_w,w0,t);grid
w_tilde= initial(sys_w,w0,t);

%RAPPRESENTAZIONE SPAZIO DI STATO DEL PROCESSO SOGGETTO A DISTURBO
%GENERALIZZATO: devo costruire matrici (P_tilde,Q-tilde,C_tilde)

C_tilde = -C;
P_tilde = [zeros(n,p) P];
Q_tilde = [eye(p,p),-zeros(p,d)];


%VERIFICO IPOTESI:

eig(S);
M_r = ctrb(A,B);
n_r = rank(M_r);

%PROGETTO K :(1° componente della legge di controllo) 

K = place (A,B,[-0.1 -0.2]);
K = -K;
eig(A+B*K)

%PROGETTO L (2° componente della legge di controllo) 
%devo torvare le matrici PI e TAU tali per cui soddisfo tali condizioni 
%1°condizione: PI = A*PI + B*TAU
%2°condizione: 0 = C_tide*PI + Q-tilde
%Per semplicita unisco queste due equazioni in una sola 
%    J*X*S- W*X+R_tilde
%     X= [PI;TAU]
%     J= [I 0; 0 0]
%     S= ECOSISTEMA 
%     W= [A B;C_tilde 0]
%     R_tilde= [P_tilde; Q_tilde]

%RISOLUZIONE CON YAMILP

R_tilde= [P_tilde; Q_tilde];
W= [A B;C_tilde zeros(p,m)];
J= [eye(n) zeros(n,m);zeros(p,n) zeros(p,m)];

X = sdpvar(n+m,r,'full');

F =[(W*X+R_tilde)-(J*X*S)==0];
dignostic = optimize(F);
check(F)
X=value(X);


PI  = X(1:n,:);
TAU = X(n+1:n+m,:);

L = TAU -K*PI;

sys_e = ss(A+B*K,B*(TAU-K*PI)+P_tilde,C_tilde,Q_tilde);
e= lsim(sys_e,w_tilde,t,zeros(n,1));

%RAPPRESENTAZIONE GRAFICA 

figure
plot(t,e(:,1));grid

figure
plot(t,e(:,2));grid


%simulazione del sistema di controllo a ciclo chiuso la cui uscita Ë
%la variabile y(t)
sys_y=ss(A+B*K,B*(TAU-K*PI)+P_tilde,C, zeros(p,r)); %Qy
y=lsim(sys_y,wtilde,t,zeros(n,1));
figure
plot(t,y(:,1),'b',t,wtilde(:,1),'r'), grid, title('Portata')
figure
plot(t,y(:,2),'b',t,wtilde(:,2),'r'), grid, title('Concentrazione')
















