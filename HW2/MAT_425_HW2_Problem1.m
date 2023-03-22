clear;
%% Problem 1
% initial conditions
x1 = 3;
x2 = 5;
xm = (x1+x2)/2;
tol = eps;%2^(-48);

% Braket corrector
while poly(x1)*poly(x2) > 0
    x1 = 2*(x1-xm) + xm;
    x2 = 2*(x2-xm) + xm;
end

% Bisection method
iterations = 0;
while x2-x1 > tol
    iterations = iterations+1;
    xm = (x1+x2)/2;
    fm = poly(xm);
    f1 = poly(x1);
    f2 = poly(x2);
    if f1 * fm < 0
        x2 = xm;
        %f2 = poly(x2);
    elseif f1 * fm >= 0
        x1 = xm;
        %f1 = poly(x1);
    end
end
disp("Iterations taken")
disp(iterations)
disp("Interval size")
disp(x2-x1)
disp("Approximate zero")
disp(xm)
%% Functions
function f = poly(x)
    f = x^3-2*x+2;
end