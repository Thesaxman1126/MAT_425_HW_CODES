%This code computes the Lorenz equations comparing RK4, RK4/5 and ode45
clear all; clc; close all;

t_interval=[0 40]; %t_interval
y0=[-12.7121529 -33.3189869 57.7069384]';
[t1,y] = ode45(@(t,y) dydt(y,t), t_interval, y0);
y=y';
yode45=y';
plot3(y(1,:),y(2,:),y(3,:),'r')
figure(1)
hold on;
scatter3(y(1,end),y(2,end),y(3,end),'k','filled')
%
t_interval=[0 40]; %t_interval
y0=[-12.7121530 -33.3189869 57.7069384]';
[t1,y] = ode45(@(t,y) dydt(y,t), t_interval, y0);
y=y';
yode45=y';
plot3(y(1,:),y(2,:),y(3,:),'b')
scatter3(y(1,end),y(2,end),y(3,end),'m','filled')




function k=dydt(y,t)
P=10; r=28; b=8/3;
A=[-P P 0;r -1 -y(1);y(2) 0 -b];
k=A*y;
end