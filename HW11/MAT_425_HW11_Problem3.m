clear all; clc; close all;
%% Problem 3
N = 100;
x = linspace(0, pi, N);
y = linspace(0, 1, N);
[X,Y] = meshgrid(x,y);
L = (sinh(pi*X/2).*cos(pi*Y/2)/sinh(pi^2/2));
phi = (sin(X).*cos(pi*Y/2));
test1 = (sin(X).*cos(pi*Y/2));
% plotting exact solution
u_exact = L - phi/(1+pi^2/4);
figure
contourf(X,Y,u_exact,20)

figure
mesh(u_exact)
%% Functions
function [soln] = Poisson_Jacobi(u, phi, Delta)
    
end