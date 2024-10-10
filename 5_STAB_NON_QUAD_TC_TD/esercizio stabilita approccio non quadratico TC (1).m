
%STABILITA' Approccio non quadratico
mu=1.6666; %mu=1.67 NO soluzione
A1=[-4 4;-5 0]-mu*[-2 2 ;-1 4];
A2=[-4 4;-5 0]+mu*[-2 2 ;-1 4];

n=2;

P1=sdpvar(n);
P2=sdpvar(n);
F=sdpvar(n,n,'full');
G=sdpvar(n,n,'full');

F1=[P1>=0];
F1=[F1, P2>=0];
F1=[F1, [A1'*F'+F*A1  P1-F+A1'*G; P1-F'+G'*A1 -(G+G')]<=0];
F1=[F1, [A2'*F'+F*A2  P2-F+A2'*G; P2-F'+G'*A2 -(G+G')]<=0];
diagnostic=optimize(F1)
check(F1)

P1=value(P1);
P2=value(P2);
F=value(F);
G=value(G);
eig(P1)
eig(P2)
eig([A1'*F'+F*A1  P1-F+A1'*G; P1-F'+G'*A1 -(G+G')])
eig([A2'*F'+F*A2  P2-F+A2'*G; P2-F'+G'*A2 -(G+G')])
