% Per assegnare degli autovalori a ciclo chiuso con stato NON accessibile,
% devo prima partire a ciclo aperto ed inizializzare le variabili:
% n-->dim. stato         A:(3,3)
% m-->dim.ingresso       B:(1,3)
% q-->dim.uscite         C:(3,1)
% 
% SISTEMA A CICLO APERTO
% dx(t) = A*x(t)+B*u(t)
% y(t) = C*x(t)}%

A=[0.1 0 0;1 -1 0;0 0.2 -0.5];
n=3; 
B=[1;0;0]; 
m=1; 
C=[1 0 0];  
q=1; 

% SISTEMA A CICLO CONSIDERANDO STIMA s_x(t) ED ERRORE e(t)
% dx(t) = (A-B*K)*x(t)+B*e(t)+B*e(t)
% de(t)=(A-B*L)*e(t)
% e(t) = x(t) - s_x(t)

% tramite calcoli mi sono riduco ad avere questa forma ,cosi da avere le 
% informazioni per trivare le matrici a ciclo chiuso (Af,Bf,Cf).
% Il principio di separazione consiste nel modificare la matrice di
% trasferimento (A-B*K) scomponendola ma gli autovalori sono gli stessi,cioè:
%        spettro(A) = spettro(A-B*K) + spettro(A-L*C)
% Ora posso costruire le matrici (Ac,Bc,Cc) a ciclo chiuso 
% Af =[(A-B*K) B*K,0 (A-L*C) ] ha questa forma:
% -posso trovare autovalori(A-B*K) se (A,B) raggiungibile
% -posso trovare autovalori(A-L*C) SE (A,C) osservabile 


%Proprieta strutturali richieste: ragg(A,B) E oss(A,C)

Mr=ctrb(A,B);
rank(Mr)  %ragg
Mo=obsv(A,C);
rho=rank(Mo) %non oss

% Decompongo secondo l'osservabilità

[a b c T l]=obsvf(A,B,C);
a_no=a(1:2,1:2);
eig(a_no)   %-1 -0.5
%(a_o,c_o)  sottospazio osservabile
a_o=a(3,3);
c_o=c(1,3);

% Nella costruzione di L e k posixioni gli autovalori nel semipiano 
% sinistro con comando place , ne troppo vicino a zero ne troppo lontano
%Costruisco matrice L,con sistema modificato ma ora è osservabile

Lot=place(a_o',c_o',-2);
Lo=Lot';
L_f=[zeros(n-rho,q);Lo];
L=inv(T)*L_f;
eig(A-L*C)  %-1 -0.5 -2

%Costruisco matrice L, sistema gia raggiungibile 

K=place(A,B,[-0.4 -0.6 -0.8]);
eig(A-B*K)


%Costruzione del sistema a ciclo chiuso 
Ac=[A -B*K;L*C A-B*K-L*C];
Bc=[B;B];
Cc=[C zeros(q,n)];

% Con eig (Ac) vedo se gli autovaolori sono a parte reale
% negativa ciò sistema è stabile 

eig(Ac)
