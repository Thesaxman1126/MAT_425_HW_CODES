clear all
%% Problem 3

% Initial Coditions
N = 9;

% Construct A
A=4*eye(N); %Diagonal for the matrix
for i=1:N-1
    A(i,i+1)=1;
    A(i+1,i)=1;
end 
A(N,1) = 1;
A(1,N) = 1;
b = zeros(N,1);
for i=1:N
    b(i) = sin(2*pi*i/10); % d_i = sin(2 i pi/10)
end
b0 = b;
x0 = zeros(N,1); % Initial guess
tic
omega0_5 = SOL(A, b, x0, 0.5, 1e-10)
toc
error0_5 = A\b0-omega0_5
tic
omega1 = SOL(A, b, x0, 1, 1e-10)
toc
error0_1 = A\b0-omega1
tic
omega1_5 = SOL(A, b, x0, 1.5, 1e-10)
toc
error1_5 = A\b0-omega1_5
function xn = SOL(A, b, x0, omega, tol)
% Extract info from A%
D1=diag(A); 
D=diag(D1); % creates D
N=size(A,1);
L=ones(N,N);
U=L; 
for i=1:N % Create U and L
    L(i,i:N)=0;
    U(i,1:i)=0;
end
L=L.*A; U=U.*A; % Define the splitting matrices

x=x0;
i=0;
rhs=(b-U*x);
xn = zeros(N,1);
xn(1)=rhs(1)/D1(1);
for j=2:N
    xn(j)=1/D1(j)*(rhs(j)-L(j,1:j-1)*xn(1:j-1));
end
while norm(xn-x)>tol
    x=xn;
    rhs=(b-U*x);
    xn(1)=rhs(1)/D1(1);
    for j=2:N
        xn(j)=omega/D1(j)*(rhs(j)-L(j,1:j-1)*xn(1:j-1));
    end
    i=i+1;
    if i>10000
        break;
    end
    
end
disp("Num iterations "+i)
end