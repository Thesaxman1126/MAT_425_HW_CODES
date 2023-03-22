clear;
%% Problem 2
% initial conditions 
x1 = 3;
x2 = 5;
tol = eps;
iteration = 0;
% Secant method
while abs(x2-x1) > tol
   iteration = iteration + 1; 
   x3 = x2 -((x2-x1)*poly(x2))/(poly(x2) - poly(x1));
   x1 = x2;
   x2 = x3;
end
disp("Number of iterations")
disp(iteration)
disp("Approximate zero")
disp(x2)

%% Functions
function f = poly(x)
    f = x^3-2*x+2;
end
%% Analysis
% This took 40 iterations while bisection took about 56 iterations