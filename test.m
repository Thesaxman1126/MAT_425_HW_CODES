clear all; clc; close all;
%eps y''+2y'+exp(y)=0, y(0)=y(1)=0
%Approx solution y=ln(2/(x+1))-ln2exp(-2x/eps)
Mesh = [100, 200, 500]; 

for N = Mesh
x=linspace(0,1,N+1);
tol=1e-6;
eps=.02;
dx=1/N;
D2=-2*eye(N-1);
for i=1:N-2
    D2(i,i+1)=1;
    D2(i+1,i)=1;
end
D2=D2/dx^2;
D1=zeros(N-1);
for i=1:N-2
    D1(i,i+1)=1;
    D1(i+1,i)=-1;
end
D1=D1/2/dx;
A=eps*D2+2*D1;

y=zeros(N-1,1);
F=Fy(A,y);
J=Jaco(A,y);
yprev=y;
y=y-J\F;
while norm(yprev-y)>tol
    yprev=y;
    F=Fy(A,y);
    J=Jaco(A,y);
    y=y-J\F;
end
yy=[0;y;0];
hold on
plot(x,yy);

end
plot(x,log(2./(x+1))-log(2)*exp(-2*x/eps),'k')
legend(["N = 100", "N = 200", "N = 500", "Asypmtotic Aprrox"])

function F=Fy(A,y)
F=A*y+exp(y);
end

function J=Jaco(A,y)
J=A+diag(exp(y));
end




