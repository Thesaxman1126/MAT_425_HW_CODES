%This code computes the Lorenz equations comparing RK4, RK4/5 and ode45
clear all; clc; close all;

t_interval=[0 40]; %t_interval
y0=[-12.7121529 -33.3189869 57.7069384]';

% r = 28
t_interval=[0 40]; %t_interval
y0=[-12.7121529 -33.3189869 57.7069384]';
[t1,y] = ode45(@(t,y) dydt(y,t,28), t_interval, y0);
y=y';
yode45=y';
figure()

plot3(y(1,:),y(2,:),y(3,:))%,'r')

% r = 100
t_interval=[0 40]; %t_interval
y0=[-12.7121529 -33.3189869 57.7069384]';
[t1,y] = ode45(@(t,y) dydt(y,t,100), t_interval, y0);
y=y';
yode45=y';
figure()
plot3(y(1,:),y(2,:),y(3,:))%,'r')
%scatter3(y(1,1),y(2,1),y(3,1),'g','filled')
%scatter3(y(1,end),y(2,end),y(3,end),'k','filled')

% r = 160
t_interval=[0 40]; %t_interval
y0=[-12.7121529 -33.3189869 57.7069384]';
[t1,y] = ode45(@(t,y) dydt(y,t,160), t_interval, y0);
y=y';
yode45=y';
figure()

plot3(y(1,:),y(2,:),y(3,:))%,'r')


% r = 350
t_interval=[0 40]; %t_interval
y0=[-12.7121529 -33.3189869 57.7069384]';
[t1,y] = ode45(@(t,y) dydt(y,t,350), t_interval, y0);
y=y';
yode45=y';
figure()

plot3(y(1,:),y(2,:),y(3,:))%,'r')

% r = 5
t_interval=[0 40]; %t_interval
y0=[-12.7121529 -33.3189869 57.7069384]';
[t1,y] = ode45(@(t,y) dydt(y,t,5), t_interval, y0);
y=y';
yode45=y';
figure()

plot3(y(1,:),y(2,:),y(3,:))%,'r')
function k=dydt(y,t,r)
P=10; b=8/3;
A=[-P P 0;r -1 -y(1);y(2) 0 -b];
k=A*y;
end