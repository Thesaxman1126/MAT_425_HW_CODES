clear
%% Problem 2
tic
% MATRIX CONSTRUCTION %
n = 9;
% A-TriDiag %
a_main = 4*ones(1,n);
a_s = ones(1,n-1);
Am = diag(a_main);
Am(n) = 1;
Asup = diag(a_s,1);
Asub = diag(a_s,-1);
A = Am+Asup+Asub;
A(1,n) = 1;
A(n,1) = 1;

% b-soln array %
b = zeros(n,1);
for i=1:n
    b(i) = sin(2*pi*i / 10); 
end
toc