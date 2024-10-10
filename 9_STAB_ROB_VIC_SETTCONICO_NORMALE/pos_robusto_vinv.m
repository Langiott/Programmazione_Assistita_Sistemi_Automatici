%posizionamento robusto dei poli a ciclo chiuso  nella regione D desiderata 

%A(p)=[-1 0.6;0.1 -p];  p\in [0.5 1]
%B=[0;1]
%regione D desiderata  Re[z]<=-2
A1=[-1 0.6;0.1 -0.5];
B1=[0;1];

A2=[-1 0.6;0.1 -1];
B2=[0;1];

n=2;
m=1;

alfa_1=2;

X=sdpvar(n);
Y=sdpvar(m,n,'full');
F=[X>=0];
F=[F, 2*alfa_1*X+A1*X+B1*Y+X*A1'+Y'*B1'<=0];
F=[F, 2*alfa_1*X+A2*X+B2*Y+X*A2'+Y'*B2'<=0];

diagnostic=optimize(F)
check(F)

X=value(X);
Y=value(Y);

K=Y*inv(X);


%simulazione
alfa1=rand(1);
alfa2=1-alfa1;
Aalfa=alfa1*A1+alfa2*A2;
Balfa=alfa1*B1+alfa2*B2;

eig(Aalfa+Balfa*K)
%%
% VADO A RAPPRESENTARE LO SFORZO DI CONTROLLO DEL SISTEMA , IN SPAZIO DI STATO CONIL 
% COMANDO ss, perche voglio vedere cosa succede a u, diventa quindi la nostra uscita.
% con condizione iniziale X(0)!=0

% dx=(A(alpha )+B(alpha)*K)
% U= K*x

sys_u=ss(Aalfa+Balfa*K,zeros(n,1),K,zeros(m,1));
x0=[0.2;0.1];
t=0:0.1:5;
figure;
initial(sys_u,x0,t);
%vedo l evoluzione libera del sistema, l uscita è proprio lo sfozo di controllo u=K*x
%Il sistema è quadraticamete stabile , quindi anche lo sforzo di controllo
%si assesta intorno all'origine (VEDI HRAFICO!!)


%%%%%%%%% vincolo sullo sforzo di controllo  |u(t)|<=2=mu (upper bound)

% sto imponendo oltre ad  Re[z]<-4 , voglio imporre un altra condizione |u(t)|<=2
% abbiamo la norma perche u  è uno scalare 
% NON VA BENE !!!!! non rispetto la seconda condizioni , quindi devo
% modificare le condizioni di stabilita aggiungendo tale condizione. Sono 
% condizioni delle sforzo di controllo 
mu=2;

X=sdpvar(n);
Y=sdpvar(m,n,'full');
F=[X>=0];
%vincolo robuto
F=[F, 2*alfa_1*X+A1*X+B1*Y+X*A1'+Y'*B1'<=0];
F=[F, 2*alfa_1*X+A2*X+B2*Y+X*A2'+Y'*B2'<=0];
%mantengo lo sforzo di controllo |u(t)|<=2
F=[F, [1 x0';x0 X]>=0];%Q=X
F=[F, [X Y';Y mu^2*eye(m)]>=0];
diagnostic=optimize(F)
check(F)

X=value(X);
Y=value(Y);

K=Y*inv(X);

alfa1=rand(1);
alfa2=1-alfa2;
Aalfa=alfa1*A1+alfa2*A2;
Balfa=alfa1*B1+alfa2*B2;

eig(Aalfa+Balfa*K)

sys_u=ss(Aalfa+Balfa*K,zeros(n,1),K,zeros(m,1));
x0=[0.2;0.1];
t=0:0.1:5;
figure;
initial(sys_u,x0,t);
%vado a vederel'evoluzione libera 
%%Ho un transitorio piu lento perche devo imporre piu condizioni 


 
