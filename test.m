clear all

tic

% constants
N=29; 
x=linspace(0,pi,N+1);
x=x';
y=linspace(0,1,N+1); 
dx=pi/(N);
dy=1/(N);
[X,Y] = meshgrid(x,y);
b=cos(pi*y/2);
U=0*Y;
U(N+1,:)=b;
f=sin(X).*cos(pi*Y/2);
error=1;
k=0;

% method
while error>1e-6
    Uinterior=U(2:N,1:N);
    for i=2:N
        U(i,1)=(dy^2*U(i-1,1)+dy^2*U(i+1,1)+2*dx^2*U(i,2)-dx^2*dy^2*f(i,1))/(2*(dx^2+dy^2));
        for j=2:N
            U(i,j)=(dy^2*U(i-1,j)+dy^2*U(i+1,j)+dx^2*U(i,j-1)+dx^2*U(i,j+1)-dx^2*dy^2*f(i,j))/(2*(dx^2+dy^2));
        end
    end
    error=norm(Uinterior-U(2:N,1:N));
    k=k+1;
end
U=U';

k

% plots
figure
hold on
contourf(X,Y,U); shading interp
%contourf(X,Y,U);
hold off

figure
hold on
Z = (sinh(pi*X/2).*cos(pi*Y/2))/sinh(pi^2/2) - (sin(X).*cos(pi*Y/2))/(1+pi^2/4);
contourf(X,Y,Z); shading interp
hold off

figure
hold on
contourf(X,Y,Z-U); shading interp
hold off

toc