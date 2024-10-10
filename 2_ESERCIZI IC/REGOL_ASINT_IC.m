A= [-0.1 0;0 -0.2];
n=2; %dim stato
B=[1 1; -0.25 0.5];
m=2; %numero ingressi
C=[0.01 0;0 1];
p=2;
P= [0 0;0.015 0.005];
d=2; %dim w
Q=zeros(p,d);

%COSTRUZIONE ECOSISTEMA S

%INIZAILIZZAZIONE 
S1=0;
r1=1;
S2=0;
r2=1;
omega=1;
S3=[0 -omega;omega 0];
r3=2;


S=blkdiag(S1,S2,S3);
r=r1+r2+r3;
sys_w=ss(S,zeros(r,1),eye(r),zeros(r,1));

%RAPPRESENTAZIONE GRAFICA 
w0=[0.5;1.5;0.01;0];
t= 0:0.1:10;
figure;
initial(sys_w,w0,t);grid
wtilde= initial(sys_w,w0,t);




