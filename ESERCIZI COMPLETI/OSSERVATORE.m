%SCOPO + SIATEMA ORIGINARIO:
%PROGETTO UN OSSERVATORE ASINTOTICO DELLO STATO PER UN PROCESSO
%INCERTO POLITOPICO (td):
%     dotx = A*x+B*u
%     u = k*x_s & y=C*x
%     dot x_s = A_nom*x + B*K*X + L*(y-y_s)
%     y=C+x & y_s = C*x_s
%     A(alpha)= Sum(alpha_i*A_i)(alpha appartiene al simplesso di ordine N)

%CONDIZIONI:
%-->(A,B) RAGGIUNGIBILE
%-->(A,C) OSSERVABILE 

%SISTEMA [x,x-x_s]:
%  dx   = (A(alpha)-B*K)  *  (x_s) +       B*K      *(x-x_s)
%  dx-dx_s = (A(alpha)-A_norm)* (x_s) + (A(alpha)-LC)  *(x-x_s)
%  A_norm = Sum(A1+A2+....+An)/n

%SISTEMA [x,x_s]:
%  dx_s    = (A_norm-B*K)  *  (x_s) +       L*C      *(x-x_s)
%  dx-dx_s = (A(alpha)-A_norm)* (x_s) + (A(alpha)-LC)  *(x-x_s)
%
%  A_norm = Sum(A1+A2+....+An)/n
%CON ESTENSIONE PILITOPICA LE MATRICE NON DIPENDONO DA ALPHA:
% A(alpha)---> A_i i=1,2,...,l

%INIZIALIZZAZIONE VARIABILI DEL SISTEMA:
A1= [0.8 0.3; 0.1 -1.5];
A2 =[0.8 0.3;0.1 1];
n=2;
B=[0;1];
m=1;
C=[1 1];
q=1;
A_norm= (A1+A2)/n;

%PROGETTO L: (A(alpha)*LC) sia quadraticamente stabile per ogni alpha
%La disequazione di Lyapunov da rispettare sarà:
% <--> S>0    (A(alpha)+LC)'*S*(A(alpha)+LC)-S<0   
% <--> -S-(A(alpha)+LC)'*S*(inv(S))*S*(A(alpha)+LC)<0  CON S=inv(P)
% <--> D_i = | -S           S*A_i-Z*C; |   
%            |A_i'*S-Z'*C'    -S       |<0
%    CON A(alpha) = sum(alpha_i*A_i) (per ogni alpha appart al simplesso-n)
%        alpha_1+alpha_2+......+alpha_n=1  alpha_i>0 i=1,..,n
%========>      L= INV(S)*Z  & S>=0
%
%oss: NON VALE IL PRNCIPIO DI SEPARAZIONE!!!

S=sdpvar(n);
Z=sdpvar(n,q);

F=[S>=0];
F=[F,[-S S*A1-Z*C;A1'*S-C'*Z' -S]<=0];
F=[F,[-S S*A2-Z*C;A2'*S-C'*Z' -S]<=0];
sol=optimize(F)
check(F)
S=value(S);
Z=value(Z);
L=inv(S)*Z;

%VERIFICA |z|<1:
alpha_1=rand(1);
alpha_2=1-alpha_1;
A_alpha= alpha_1*A1+alpha_2*A2;

abs(eig((A_alpha+L*C)))


%PROGETTO K:(A_hat+B_hat*K_hat) sia quadraticamente stabile 
% |A_norm-B*K              L*C|=|A_norm                  L*C|+|B|*|K 0|
% |A_alpha-A_norm  A_alpha-L*C| |A_alpha-A_norm  A_alpha-L*C| |0|
%
% = A_hat+B_hat*K_hat
%
%la disequazione di Lyapunov da rispetattare :
%   (A_alpha+B*K)*P*(A_alpha+B*K)-P<0
% <--> ESISTE Q=inv(P) & P=P'>0
% <-->       | -Q           A_i*Q-B*Y; |
%            |Q*A_i'-Y'*B'    -Q       |<0
% ===>        K=Q*INV(Y) & Q>=0
%          Q=|Q1 0;Q2 0|  & Y=|Y1 0|


A_norm= (A1+A2)/n;
A1_hat= [A_norm L*C; A1-A_norm A1-L*C];
A2_hat= [A_norm L*C; A1-A_norm A1-L*C];
B_hat= [B;zeros(n,m)];

Q1=sdpvar(n);
Q2=sdpvar(n);
Y1=sdpvar(m,n,'full');
Q=blkdiag(Q1,Q2);
Y=[Y1 zeros(m,n)];

%clear F
F1=[F1,Q>=0];
F1=[F1,[-Q (Q*A1_hat'-Y'B_hat');(A1_hat*Q-B*Y) -Q]<=0];
F1=[F1,[-Q (A2_hat*Q-B*Y);(Q*A2_hat'-Y'B_hat') -Q]<=0];
soluzione_1=optimize(F1);
check(F1)

Y=value(Y);
Q=value(Q);

K_hat=Y*inv(Q);

%VERIFICA |Z|<1

eig(A1_hat-B_hat*K_hat);
eig(A2_hat-B_hat*K_hat);

% sum(alpha_i)=1
a1= rand(1);
a2= 1-a1;

A_alpha_hat = a1*A1_hat+a2*A2_hat;

Af=[A_norm-B_hat*K_hat L*C;A_alpha-A_norm A_alpha-L*C];
abs(eig(Af))
%abs(eig(A_alpha_hat+B_hat*K_hat));




