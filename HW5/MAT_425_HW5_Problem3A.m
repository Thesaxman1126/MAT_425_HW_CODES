clear all
%% PROBLEM 3A ANALYSIS PLEASE READ!
% What goes wrong with this choice of Delta-t is that when doing the
% Crank-Nicholason method you'd get the growth factor to be zero meaning
% you'd get the solution to be y_0 = 1 y_n = 0 which is not particulaly
% useful. We also see that the implicit Euler method over estimates the 
% exact solution as it has a growth factor of 1/3.
%% PROBLEM 3 CODE
% INITIAL CONDITIONS %
dt = 2; % Time-step
a = -1; % sets up y' = -y
max_iter = 10; % Max time-steps
y = zeros(1,max_iter+1); % initialize y array
y(1) = 1; % y initial condition
t = 0:2:20; % t values for plotting
% IMPLICIT EULER % 
for i = 1:max_iter-1
   y(i+1) = y(i)/(1-a*dt); 
end
hold on
figure(1)
title("Numerical solns to y'=-y using dt = 2")
plot(t,y,'b')

% CRANK-NICOLSON %
y = zeros(1,max_iter+1); % initialize y array
y(1) = 1; % y initial condition
for i = 1:max_iter-1
   y(i+1) = ((1+a*dt/2)/(1-a*dt/2))*y(i); 
end
plot(t,y,'r')
% EXACT %
tt = 0:0.1:20;
plot(tt,exp(-tt),'m')
xlabel("t")
ylabel("y")
legend("IMPLICIT EULER","CRANK-NICOLSON", "EXACT")
hold off