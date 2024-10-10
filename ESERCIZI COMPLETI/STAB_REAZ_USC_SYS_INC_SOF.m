%SCOPO:
%STABILIZZAZIONE CON REAZIONE STATICA DALL'USCITA A TEMPO DISCRETO 
%DI UN SISTEMA INCERTO, ANCHE DETTO SOF. sono nel caso non ho nesuna 
%informazione sullo stato.

%SISTEMA:
%  | x(k+1) = A(p)*x+B*u     x--->dim n
%  | y(k) = C*x              y--->dim q
%  | u= K*y=K*C*x            u--->dim m 
%  | A(p)---->A_alpha = SUM(alpha_i*A_i) 

%FUNZIONE DI LYAPNOV DA RISPETTARE TD:
%Si costruisce tai LMI's per avere una stabilizzazione mediante un
%controllore statico dell uscita, la lmi sarà:
%   Esiste W=W'<=0 : (A+BKC)'*W*(A+BKC)-W<0 
%
% <---> (K=N*INV(M)) & (M*C=C*W) & (W>0)
% ---->(A+BNM^-1M^-1CW)'*(W)*(A+BNCW)-W<0
% ----> W*(A+BNCW)'*(-W^-1)*(A+BNCW)<0
% ----> | -W    W(A+BNCW)'|          oss:complemento si schur
%       |(A+BNCW)W     -W |<0
% ----> |-W     WA_i'+C'N'B'|
%       |A_iW+BNC        -W |<0      oss:estensione incertezza politopica
%
% u = N*M^-1*y; --> legge di controllo che stabilizza il sistema 
%oss: Ho definito che (A+BKC)' sia asintoticamente stabile e per ipotesi
%     di stazionarieta anche (A+BKC) lo è

%FUNZIONE DI LYAPUNOV DA RISPETTARE TC
%Il sistema stazionario:
%   |dx = A*x+B*u
%   |y = C*x
%   |u = k*y
%E' stablizzabile con un controllore u=ky:
%<---> Esiste W=W'>0 : (A+BK)'W+W(A+BK)<0
%<---> K = N*INV(M) & M*C=C*W (C=INV(M)*C*W)& W>0
%      (A+BK)'W+W(A+BK)<0
%<---> (A+BN(M^-1)(M^-1)CW)'W+W(A+BN(M^-1)(M^-1)CW)<0
%<---> WA'+ C'N'B'+ AW + BNC <0 
% caratterizzando su matlab scriverò, con stensione incertezza politopica:
% S=[S,W*A1'+ C'*N'*B'+ A1*W + B*N*C <0 ];
% S=[S,W*A2'+ C'*N'*B'+ A2*W + B*N*C <0 ];


%Inizializzazione variabile di sistema 
A1=[1.5 0.5;-1 0.2];
A2=[1.5 0.5; 1 0.2];
n=2;
B==B1==B2==[0;1];
m=1; 
C=[1 1];
q=1; 

%Costruisco le lmi's:
W=sdpvar(n);
M=sdpvar(q,q,'full');
N=sdpvar(m,q);
S=[W>=0];
S=[S,M*C-C*W==0];
S=[S,[-W (W*A1'+C'*N'*B');(A1*W+B*N*C)  W]<=0];
S=[S,[-W (W*A2'+C'*N'*B');(A2*W+B*N*C)  W]<=0];

soluzione= optimize(S);
check(S);
M=value(M);
N=value(N);
W=value(W);

K=N*inv(M);

%Verifica 

alpha_1= rand(1);
alpha_2= 1-alpha_1;
A_alpha = alpha_1*A1+alpha_2*A2;
B_alpha = alpha_1*B1+alpha_2*B2;

eig (A_alpha+B_alpha*K*C);







