clear 
clc
%% Problem 2
% Initial conditions 
x = [0,1,2,3,4];
y = [0,1,4,1,0];
h = 1;
n = length(x)-2; % determine size of A

% SET UP OF TRI-DIAGONAL MATRIX
% Main diag of tri-diag
a_main = 4*ones(1,n);
a_main(1) = 5;
a_main(n) = 5;

% Sub/Super diag of tri-diag
a_s = ones(1,n-1);

% Make tridiag matrix
Am = diag(a_main);
Asup = diag(a_s,1);
Asub = diag(a_s,-1);
A = Am+Asup+Asub;

% SOLUTION VECTOR 
b = zeros(n,1);
for i=1:n
    b(i) = y(i)-2*y(i+1)+y(i+2); 
end
b = (6/(h^2)) * b;

% Solve A*x = b
fxx = linsolve(A,b);

% Store f'' values along with the parabolic runout conditions
f_xx = zeros(1,n+2);
for i = 2:4
   f_xx(i) = fxx(i-1);
end
f_xx(1) = fxx(1);
f_xx(5) = fxx(n);

% Constructing the splines
xx = min(x):0.01:max(x);
f = zeros(n+1,401);
for i = 1:4
    f(i,:) = (f_xx(i)/6)*((x(i+1)-xx).^3 - (x(i+1)-xx)) + (f_xx(i+1)/6)*((xx-x(i)).^3 - (xx-x(i)))+y(i)*(x(i+1)-xx) + y(i+1)*(xx-x(i));
end

% Graphing the splines
hold on
i = 1;
k=101;
while i <= n+1
    plot(xx(k-100:k),f(i,k-100:k))
    k = k+100;
    i = i+1;
end
plot(x, y, "ko", 'MarkerFaceColor', 'k')
title("Cubic Spline interpolation with Quadratic Run-Out Conditions")
legend("f0","f1","f2","f3","Raw Data")