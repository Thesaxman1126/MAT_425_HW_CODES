clear all 
%% Problem 1
tic % start timer
% Constants 
dt = 1000; 
t = 0:dt:200000;
y = zeros(1,length(t));
y(1) = 0.00001;
y(2) = y(1);

% Solve y' = y^2 - y^3 using BDF2
for i = 3:length(t)
    y(i) = NR_Solve(y(i-1), y(i-1), y(i-2), dt, eps, 500000,i);
end
% to impliment: exact solution %
% Plotting
figure(1)
hold on
title("Solution y' = y^2-y^3, y(0) = 1e-5 Using BDF2")
plot(t,y,'r')
ylim([-0.5 1.5])
xlabel('t')
ylabel('y(t)')
hold off
toc % end timer


%% NR-Method for y_{n+1}
function yp1 = NR_Solve(y0, yn, ym1, dt, tol, max_iter,j)
    y_prev = y0;
    [g, dg] = G(y_prev, yn, ym1, dt);
    yp1 = y_prev - g/dg;
    iter = 1;
    
    while abs(yp1 - y_prev) > tol
        y_prev=yp1;
        [g, dg] = G(y_prev, yn, ym1, dt);
        yp1 = y_prev - g/dg;
        iter = iter+1;
        
        if iter > max_iter
            value = j
            disp("Did not converge")   
            abs(yp1 - y_prev)
            break 
        end
    end
end
% the function and its derivative which we're finding the root for
function [g, dg] = G(x, yn, ym1, dt)
    g = x - (4/3)*yn + (1/3)*ym1 - dt*(2/3)*(x^2 - x^3);
    dg = 1 - dt*(2/3)*(2*x - 3*x^2);
end

