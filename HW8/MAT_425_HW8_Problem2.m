clear all; close all; clc;

%% Problem 2

% constants
N = 100; dx = 1/N; % Number of gridpoints and x-grid spacing
x1 = 1; x2 = 2; y1 = 1; y2 = 0; % Boundary values
x = 1:dx:2; % spacial array for x
tol = eps; % Tolerance
max_iter = 1e3; % Max number of iterations

% D2 Matrix
D2 = (-2*eye(N-1) + diag(ones(1,N-2),-1)+diag(ones(1,N-2),1)) / (dx^2);

% D1 Matrix
D1 = (diag(ones(1,N-2)./x(2:N-1), 1)-diag(ones(1,N-2)./x(1:N-2), -1)) / dx;

% A for BVP method
A = D2 + D1;

tic
[soln, F_Final, Jacobian_Final] = solve_BVP_Jacobian(A, y1, y2, N, dx, tol, max_iter);
toc
figure()
hold on
plot(x, soln,"r")
title("y''+(2/x)y'-e^y = 0, y(1) = 1, y(2) = 0")
xlabel("x")
ylabel("y(x)")
hold off
F_Final;
Jacobian_Final;
%% Functions
function [F] = F_y(A, y)
    F = A*y - exp(y);
end

function [J] = Jacobian(A, y)
    J = A - diag(exp(y));
end

function [soln, F_final, J_final] = solve_BVP_Jacobian(A, y1, y2, N, dx, tol, max_iter)
    y_in = zeros(N-1, 1); % set up interior points
    y_prev = y_in;
    F = F_y(A, y_in);
    J = Jacobian(A, y_in);
    E = zeros(N-1,1);
    E(1,1) = 2*y1/dx - y1/dx^2;
    y_in = y_in + J\(E-F);
    iter = 1;
    while norm(y_prev - y_in) > tol
        y_prev = y_in;
        F=F_y(A, y_in);
        J=Jacobian(A, y_in);
        y_in=y_in + J\(E-F);
        iter = iter + 1;
        
        if iter > max_iter
            disp("Did not converge within tolerance given")
            disp("norm of before kill sequence")
            norm(y_prev - y_in)
            break
        end
    end
    iter
    soln = [y1; y_in; y2];
    F_final = F;
    J_final = J;
    norm(y_prev - y_in) 
end