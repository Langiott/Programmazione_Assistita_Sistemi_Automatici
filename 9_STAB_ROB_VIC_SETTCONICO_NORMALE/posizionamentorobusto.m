%posizionamento robusto nella regione D di stabilita  Re[z]<-4  alfa_1=4

%A(p)=[1 0.6;0.8 p]  B=[0;1]  p \in [-1 1]

%(A(alfa),B(alfa))=\sum_i alfa_i (A_i,B_i)
A1=[1 0.6;0.8 -1];
B1=[0;1];

A2=[1 0.6;0.8 1];
B2=[0;1];
%B=B1=B2 
%LA SOMMATORIA DI ALPHA CON I UGUALE A 1, PROPRIETA DEL SIMPLESSO
n=2;
m=1;

alfa_1=4;

%IMPONGO LE CONDIZIONE DA RISPETTARE 
X=sdpvar(n);
Y=sdpvar(m,n,'full');
%LE LMI, PRIMA VONCOLI DI X E POI LE CONDIZIONI AI VERTICI 
F=[X>=0];
F=[F, 2*alfa_1*X+A1*X+B1*Y+X*A1'+Y'*B1'<=0];
F=[F, 2*alfa_1*X+A2*X+B2*Y+X*A2'+Y'*B2'<=0];

diagnostic=optimize(F)
check(F)

%Le LMI danno risultati positivi 

X=value(X);
Y=value(Y);

%stablizza degli autovaori della matrrice trasposta , che cprriswponde
%anche a stabilizzare la matrice non trasposta 
K=Y*inv(X);
%valori negativi 

alfa1=0.4;
alfa2=0.6;
Aalfa=alfa1*A1+alfa2*A2;
Balfa=alfa1*B1+alfa2*B2;

eig(Aalfa+Balfa*K)
%vado a vedere gli autovaòori e noto che sono negativi e  Re[z]<-4

alpha1 = rand(1);
alpha2 =1- alpha1;
%il risultato dovrebbe verificare sempre che autovalori abbiabno  Re[z]<-4


%%
%%%%%%%%%%%%FINE %%%%%%%%%%%%%%
%con ilcomando place vado ilo a mettere gli autovaolri 
% con le lmi posso sia mettere igli autovalorie posso imporre anche delle
% condizioni delle variabili interne alle equaizioni LMI 

%posizionamento robusto nella regione D di stabilita |z|<r 

%A(p)=[-1.2 p;1.1 0.9]  B(p)=[1;p]  p \in [0 0.5]

%(A(alfa),B(alfa))=\sum_i alfa_i (A_i,B_i)
% IL numero di parametri indica il numero di vertici...2 vertivi , 1
% parametro. nel caso generale N=2^n_p
% A1=[-1.2 0;1.1 0.9];
% B1=[1;0];
% 
% A2=[-1.2 0.5;1.1 0.9];
% B2=[1;0.5];
% 
% n=2;
% m=1;
% 
% r=0.5; %|z|<r(RAGGIO R)  
% 
% X=sdpvar(n);
% Y=sdpvar(m,n,'full');
% F=[X>=0];
% %matrice costruita 
% F=[F,[-r*X A1*X+B1*Y;X*A1'+Y'*B1' -r*X]<=0];
% F=[F,[-r*X A2*X+B2*Y;X*A2'+Y'*B2' -r*X]<=0];
% 
% diagnostic=optimize(F)
% check(F)
% % dato che sono condizioni NECCESARIE E SUFFICENTE  fissando r 
% % non ammette soluzione allora l imìntero problema non ammette soluzioni
% % devo mdificare il valore di r , o intervengo a livello harware,
% % aggiungendo un sensore, cioe aggiungo vettore lin.ind. a B
% 
% %AMMETTE SOLUZIONE!!
% 
% X=value(X);
% Y=value(Y);
% 
% K=Y*inv(X);%posiziona autovalori nella regione desiderata 
% 
% alfa1=0.3;
% alfa2=0.7;
% Aalfa=alfa1*A1+alfa2*A2;
% Balfa=alfa1*B1+alfa2*B2;

abs(eig(Aalfa+Balfa*K))

alfa1=0.6;
alfa2=0.4;
Aalfa=alfa1*A1+alfa2*A2;
Balfa=alfa1*B1+alfa2*B2;

abs(eig(Aalfa+Balfa*K))




