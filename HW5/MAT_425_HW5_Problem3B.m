clear all
%% PROBLEM 3B ANALYSIS PLEASE READ!
% Here we see that the CN method has a non-zero growth factor of -1/41
% which shows that it under estimates the exact solution of the
% differential equation. Here we see that the IE method has a growth factor
% of 10/31 which we see still over estimates the exact solution of the
% differential but not by as much.
%% PROBLEM 3 CODE
% INITIAL CONDITIONS %
dt = 2.1; % Time-step
a = -1; % sets up y' = -y
max_iter = 10; % Max time-steps
y = zeros(1,max_iter+1); % initialize y array
y(1) = 1; % y initial condition
t = 0:2.1:21; % t values for plotting
% IMPLICIT EULER % 
for i = 1:max_iter-1
   y(i+1) = y(i)/(1-a*dt); 
end
hold on
figure(1)
title("Numerical solns to y'=-y using dt = 2.1")
plot(t,y,'b')
% CRANK-NICOLSON %
y = zeros(1,max_iter+1); % initialize y array
y(1) = 1; % y initial condition
for i = 1:max_iter-1
   y(i+1) = ((1+a*dt/2)/(1-a*dt/2))*y(i); 
end
plot(t,y,'r')
% EXACT %
tt = 0:0.1:21.1;
plot(tt,exp(-tt),'m')
xlabel("t")
ylabel("y")
legend("IMPLICIT EULER","CRANK-NICOLSON", "EXACT")
hold off