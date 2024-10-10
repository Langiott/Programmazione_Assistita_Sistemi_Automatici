%posizionamento robusto con sistema con 
%matrice di incerta politopica con 2 parametri incerti

%Definizione matrici 
A1=[-1 0.1;0.5  0.1];
A2=[-1 0.1;1  0.1];
A3=[-1 0.3;0.5  0.1];
A4=[-1 0.3;1  0.1];

B=[1;0];
C=[1];

x_0= [0.2;0.1];

%matrice C é invertibile quindi stato accessibile

n=2;%dim. stato
m=1;%dim. uscita 

%regione desiderata
alfa=3;

%vincolo 
mu=5;

%definisco le lmi rispetto alla mia regione D
X=sdpvar(n);
Y=sdpavar(m,n,'full');

F=[F, 2*alfa*X+A1*X+B1*Y+X*A1'+Y'*B1'<=0];
F=[F, 2*alfa*X+A2*X+B2*Y+X*A2'+Y'*B2'<=0];
F=[F, 2*alfa*X+A3*X+B1*Y+X*A3'+Y'*B1'<=0];
F=[F, 2*alfa*X+A4*X+B2*Y+X*A4'+Y'*B2'<=0];

%Aggiunta vincolo di controllo 
Q=X;
F=[F,X>=0];
F=[F,[1 x_0';x_0 Q] >=0]; %ERRORE!!!!!!!!Q=X ,MI HA BOCCIATO PER QUESTO
F=[F,[Q Y';Y (mu*mu)*eye(n)] >=0]

diagnostic = optimize(F);
check(F);

X=value(X);
Y=value(Y);

%Posso progettare K:

K=Y*inv(X);



%Costruzione matrice politopica espressa  rispetto ai vertici 
% sum(alpha_i)<1 perche sono parametri\\ convessi 

alpha_1=0,2;
alpha_2=0,2;
alpha_3=0,3;
alpha_4=0,3;

A_alpha=alpha_1*A1+alpha_2*A2+alpha_3*A3+alpha_4*A4;
B_alpha=alpha_1*B+alpha_2*B+alpha_3*B+alpha_4*B;


eig(A_alfa+B_alfa*K)

%Rappresentazione sistema a ciclo chiuso
sys_u=ss(A_alfa+B_alfa*K,zeros(n,1),K,zeros(m,1));
x_0=[0.2;0.1];
t=0:0.1:10;
figure;
initial(sys_u,x0,t);





