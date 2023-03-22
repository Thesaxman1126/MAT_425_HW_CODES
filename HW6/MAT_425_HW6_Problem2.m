clear all
%% PROBLEM 2 
% Setup
dt = 0.01;
t = 0:dt:20;
N = size(t);
max_iter = N(2);

% E.E. Method
x_EE = zeros(N);
x_EE(1) = 1;
y_EE = zeros(N);
for i = 1:max_iter-1
    x_EE(i+1) = x_EE(i) + dt*y_EE(i);
    y_EE(i+1) = y_EE(i) - dt*x_EE(i);
end

% Symplectic Method
x_SM = zeros(N);
x_SM(1) = 1;
y_SM = zeros(N);
for i = 1:max_iter-1
    x_SM(i+1) = x_SM(i) + dt*y_SM(i);
    y_SM(i+1) = y_SM(i) - dt*x_SM(i+1);
end

% Exact soln
x_exact = cos(t);

% plotting
figure(1)
hold on 
title("Numerical Solns to x''+x=0, x(0)=1 x'(0)=0")
plot(t, x_EE,"b-")
plot(t, x_SM,"r-")
plot(t,x_exact,"k-.")
xlabel("t")
ylabel("x")
legend("EE","SM","EXACT")
hold off
