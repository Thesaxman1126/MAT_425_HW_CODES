clc; clear all; close all;
%% Problem 2
Nx = 100; % x-grid size
Nt = 19900; % t-grid size
dx = 1/Nx; % x-grid spacing
dt = 1/Nt; % t-grid spacing
a = dt/(dx^2); 
x = 0:dx:1; % x-grid including ends
xx = x(2:Nx); % x-drid without ends
D2 = -2*eye(Nx-1) + diag(ones(1,Nx-2),1) + diag(ones(1,Nx-2), -1); % 
D = eye(Nx-1) + a*D2; % Differentiation matrix 
old = 0; 
max_iter = Nt;
iter = 0; %iteration counter
u = sin(2*pi*xx);
while norm(u-old) > eps
    soln = u*D;
    hold on
    if mod(iter,1000) == 0
        plot(x,[0,soln,0])
    end
    iter = iter+1;
    old = u;
    u = soln;
    check = norm(u-old);
    if iter > max_iter
        disp("DNC")
        disp(check)
        break
    end
end