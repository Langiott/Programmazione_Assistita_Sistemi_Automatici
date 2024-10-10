gamma=0.4619; %gamma=0.47 NO

A1=[0.8 -0.25 0 1;1 0 0 0;0 0 0.2 0.03; 0 0 1 0]-gamma*[ 0 0 1 0]'*[0.8 -0.5 0 1];
A2=[0.8 -0.25 0 1;1 0 0 0;0 0 0.2 0.03; 0 0 1 0]+gamma*[ 0 0 1 0]'*[0.8 -0.5 0 1];
n=4;



P1=sdpvar(n);
P2=sdpvar(n);
F=sdpvar(n,n,'full');
G=sdpvar(n,n,'full');

F1=[P1>=0];
F1=[F1, P2>=0];
F1=[F1, [A1*F+F'*A1'-P1 -F'+A1*G;-F+G'*A1'  P1-(G+G')]<=0]
F1=[F1, [A2*F+F'*A2'-P2 -F'+A2*G;-F+G'*A2'  P2-(G+G')]<=0]
diagnostic=optimize(F1)
check(F1)

P1=value(P1);
P2=value(P2);
F=value(F);
G=value(G);

eig(P1)
eig(P2)

eig([A1*F+F'*A1'-P1 -F'+A1*G;-F+G'*A1'  P1-(G+G')])
eig([A2*F+F'*A2'-P2 -F'+A2*G;-F+G'*A2'  P2-(G+G')])
