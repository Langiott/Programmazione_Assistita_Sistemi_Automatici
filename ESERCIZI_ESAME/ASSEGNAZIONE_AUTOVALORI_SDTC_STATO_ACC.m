%%
% L'obiettivo è di progettare una legge di controllo in retroazione.
% Con la tecnica di assegnazione degli autovalori, assegno dei autovalori a 
% ciclo chiuso in modo da avere un trasitorio smorzato e rapido o per
% stabilizzare un sistema instabile (con lo studio del luogo delle radici
% dobbiamo disporre gli autovalori nel semipiano sinistro )

%%
% Per assegnare degli autovalori a ciclo chiuso con stato accessibile,
% devo prima partire a ciclo aperto ed inizializzare le variabili:
% n-->dim. stato         A:(3,3)
% m-->dim.ingresso       B:(1,3)
% q-->dim.uscite         C:(3,1)
% 
% SISTEMA A CICLO APERTO
% dx(t) = A*x(t)+B*u(t)
% y(t) = C*x(t)

A=[0.1 0 0;1 -1 0;0 0.2 -0.5];
n=3;
B=[1;0;0]; 
m=1; 
C=[1 0 0]; 
q=1; 

%%
% Studio il sistema a ciclo chiuso con legge di controllo u(t). Stimo lo
% stato x(t),la stima la chiamiamo s_x(k),l'errore lo indichiamo con
% e(t),con lo stato accessibile.L'errore converge a zero asintoticamente
% cioè nel transitorio l'errore converge a zero e(t)-->0 mentre la stima
% converge al valore esatto, s_x(t)-->x(t).


% SISTEMA A CICLO CHIUSO & LEGGE DI CONTROLLO 
% dx(t) = (A-B*K)*x(t)+B*r(t)
% y(t) = C*x(t) 
% u(k)=-K*s_x(k)+r(k) 
% e(k)=x(k)-x(k)

%%
% Devo verificare le condizioni necessarie e sufficienti di (A,B),cioè la
% raggiunggibilità e (A,C) osservabilita, determinando se le matrici
% abbiano rango massimo (det(Af)!=0)se si verifica possiamo definire che le 
% matrici siano ammissibili pertanto il sistema ammette soluzione

% Se il sistema di partenza non è ragg e/o oss ,aggiungo o scarto un
% sensore, cioè intervengo a livello hardware aggiungendo appunto un
% sensore; a livello software modifichiamo la matrice C aggiungendo un
% vettore lin. ind.

Mr=ctrb(A,B);
rank(Mr) 
Mo=obsv(A,C);
rho=rank(Mo);

V=[0 0 1];
C=[1 0 0;V];
q=2; 
Mo=obsv(A,C);
rank(Mo);

%%
% Il principio di separazione consiste nel modificare la matrice di
% trasferimento (A-B*K) scomponendola ma gli autovalori sono gli stessi,cioè:
%        spettro(A) = spettro(A-B*K) + spettro(A-L*C)
% Cosi posso costruire le matrici (Ac,Bc,Cc) a ciclo chiuso 


K=place(A,B,[-0.5 -0.8 -1]);
eig(A-B*K);
Lt=place(A',C',[-1.2 -1.3 -1.5]);
L=Lt';
eig(A-L*C);
Ac=[A -B*K;L*C A-B*K-L*C];
Bc=[B;B];
Cc=[C zeros(q,n)];
eig(Ac);

