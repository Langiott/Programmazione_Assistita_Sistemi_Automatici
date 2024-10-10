%ES: POSIZIONAMENTO ROBUSTO.Devo progettare una legge di 
% controllo u(t) tale che il sistema a ciclo chiuso sia
% D-Stabile per il domini di incertezza tale che:
%             re[autovalori]<-alpha
%osservazioni:
%->alpha>0 & alpha =1
%-> regione in questione sarà: D={f_D = 2*alpha*X+A1*X+B1*Y+X*A1'+Y'*B1'}
%->la funzione di lyapnouv sarà: V(X(t)) = x'*P*x P=X^-1 & K=Y*INV(X)


%inizializzazioni matrici del sistema

A1 =[1 1.2;0.4 1];
B1 = [1;1];

A2 = [1.5 1.2; 0.4 1];
B2 =[1.5; 1];

C = [1 1.2 ; 0 1];

x0=[0.2;0.1];

alpha= 1;
mu= 2;

%Definisco le LMI's che mi definiscono una regione di D_stabilità

X=sdpvar(n);
Y=sdpvar(m,n,'full');
F=[X>=0];
%vincolo robusto
F=[F, 2*alpha*X+A1*X+B1*Y+X*A1'+Y'*B1'<=0];
F=[F, 2*alpha*X+A2*X+B2*Y+X*A2'+Y'*B2'<=0];
%mantengo lo sforzo di controllo |u(t)|<=2
F=[F, [1 x0';x0 X]>=0];%Q=X
F=[F, [X Y';Y mu^2*eye(m)]>=0];
diagnostic=optimize(F)
check(F)

X=value(X);
Y=value(Y);

K=Y*inv(X);

%verifica delle LMI's siano ben costruite 
% sum [alpha_i] = 1

alpha_1=0.4;
alpha_2=0.6;
A_alpha=alpha_1*A1+alpha_2*A2;
B_alpha=alpha_1*B1+alpha_2*B2;
%Balpha= alpha_1*B+alpha_2+B come nel compito, B non dipende da alpha  

%Verifica modulo degli autovalori sia a aprte reale negativo
eig(Aalfa+Balfa*K)

%Rappresentazione sistema a ciclo chiuso
sys_u=ss(A_alpha+B_alpha*K,zeros(n,1),K,zeros(m,1));
x_0=[0.2;0.1];
t=0:0.1:10;
figure;
initial(sys_u,x0,t);


