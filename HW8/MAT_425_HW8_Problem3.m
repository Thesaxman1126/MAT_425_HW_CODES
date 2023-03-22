clear all; clc; close all;
%% Problem 3 
tic
% CONSTANTS %
dt = 0.00001; % timestep
t = 0:dt:1; % t-Mesh
N = length(t);% how many times we do RK4
y0 = 0; %y(0) = 0
zu = 1.25; % undershooting condition
zo = 1.5; % overshooting condition

[soln,fa,fb] = solve_BVP_Shooting(y0, zu, zo, 1, dt, eps);
toc

figure()
hold on
title("Solving BVP: y'' = y^2-1, y(0) = 0, y(1) = 1")
plot(t, soln, 'k')
plot(t, fa, 'r')
plot(t, fb, 'b')
xlabel("x")
ylabel("y(x)")
legend(["Final", "Undershoot", "Overshoot"])
hold off
%% functions
% y'' = y^2-1
% let y' = f = z
% and z' = g = y^2-1
function dydx = eval_f(zn) % evaluates f from system above
    dydx = zn;
end

function dzdx = eval_g(yn) % evaluates g from system above
    dzdx = (yn^2)-1;
end


function rk4g = solve_rk4g(yn,dt) % solves for k values related to g's RK4 
    k1 = eval_g(yn);
    k2 = eval_g(yn+(dt/2)*k1);
    k3 = eval_g(yn+(dt/2)*k2);
    k4 = eval_g(yn+(dt)*k3);
    rk4g = (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
end

function rk4f = solve_rk4f(zn,dt) % solves for k values related to f's RK4 
    k1 = eval_f(zn);
    k2 = eval_f(zn+(dt/2)*k1);
    k3 = eval_f(zn+(dt/2)*k2);
    k4 = eval_f(zn+(dt)*k3);
    rk4f = (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
end

function [soln] = RK4_Solver(y0, z0, dt, max_iter) % Solves the system above using RK4
    % allocate memory for arrays arrays
    yrk4 = zeros(1,max_iter);
    zrk4 = zeros(1,max_iter);
    yrk4(1) = y0; % populate initial condition
    zrk4(1) = z0; % populate initial condition
    
    for i = 2:max_iter % solve RK4
       yrk4(i) = yrk4(i-1) + solve_rk4f(zrk4(i-1),dt); 
       zrk4(i) = zrk4(i-1) + solve_rk4g(yrk4(i-1),dt);
    end
    
    soln = yrk4;
end
% finds correct endpoint (nice)
function [BVP_sol, fa, fb] = solve_BVP_Shooting(y0, za, zb, target, dt, tol) 
    t = 0:dt:1;
    N = length(t);
    max_iter = 500;
    iter = 0;
    fa = RK4_Solver(y0, za, dt, N);
    fb = RK4_Solver(y0, zb, dt, N);
    while zb - za > tol
        iter = iter+1;
        zm = (za+zb)/2;
        fm = RK4_Solver(y0, zm, dt, N);

        if fm(N) - target >0
            zb = zm;
        elseif fm(N) - target < 0
            za = zm;
        end
    
        if iter > max_iter
            BVP_sol = fm;
            disp("DNC after")
            iter
            break
        end
    end
    BVP_sol = fm;
    disp("Converged at")
    iter
    disp("final IVP guess")
    zm
end