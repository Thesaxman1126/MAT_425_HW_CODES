clear all

%% ANALYSIS PLEASE READ
%
% After some experimental testing I 
% believe this uses somewhere around O(19*n)
% operations. IDK How correct this is.
%
% I have a possible other way to solve this but don't
% have time to test it. 
% My idea was to use the TDMA on the n-2xn-2 sub array using sinc that is 
% tri-diagonal to get an array that looks like
%
% | b1 c1 00 ... 00 a1   | | x1   |    | b1    |
% | a2 01 00 00 ... 00   | | x2   |    | b2'   |
% | 00 00 01 00 ... 00   | | x3   | =  | b3'   |
% | :  :     \      :    | | :    | =  | :     |
% | 00 00 ... 00 01 cn-1 | | xn-1 |    | bn-1' |
% | cn 00 ... 00 an bn   | | xn   |    | bn    |
%
% From here we know x3 -> xn-2. And thus there is a much simpler system 
% that we can solve using either matrix inverses or GJ-Elimination.
% If my calculations were correct it should take about O(n+5) operations 
% given that TDMA uses O(n) operations the sub system takes O(n-2) using
% TDMA and to solve the rest using GJ-Elimination we get an extra 7
% operations giveing O(n+5) operations
% 
% The smaller system in question would be a 4x4 system 
% | b1 c1 00 a1  || x1   |   | b1    |
% | a2 01 00 00  || x2   | = | b2'   |
% | 00 00 01 cn-1|| xn-1 | = | bn-1' |
% | cn 00 an bn  || xn   |   | bn    |
%
%
%% Problem 2
tic

% Initial conditions %
N = 9; % size
% Set-up A
a = ones(N,1); 
b = 4*a;
c = a;

% d vector in Ax = d
d = zeros(N,1);
for i=1:N
    d(i) = sin(2*pi*i/10); % d_i = sin(2 i pi/10)
end
d0 = d;

% initialize cost
cost = 0;

% Solving Ax = d
[x, cost] = TDMA_Mod(a,b,c,d,N,cost);
disp("solution vector")
disp(x)

% Check soln with matlab builtins and disp cost
A = 4*eye(N);
A(1,N) = 1;
A(N,1) = 1;
for i=1:N-1
    A(i,i+1) = 1;
    A(i+1,i) = 1;
end
disp("Time to perform TDMA_Mod")
toc


tic
A\d0;
disp("Time to perform so solve Ax = d")
toc 
disp("Error between matlab builtin and TDMA_Mod")
disp(A\d0-x)
disp("Total operations in TDMA_Mod")
disp(cost)


%% Problem 1
% MODIFIED TDMA %
function [sol, comp_cost] = TDMA_Mod(a,b,c,d,N,cost)
% Set-up for column tracker     
    NCol = zeros(N,1);
    NCol(1) = a(1);
    NCol(N-1) = c(N-1);
    NCol(N) = b(N);
% Set-up for row tracker
    NRow = zeros(N,1);
    NRow(1) = c(N);
    NRow(N-1) = a(N);
    NRow(N) = b(N);
% Forward TDMA
    for i=1:N-2
        NCol(i+1) = NCol(i+1) - NCol(i) * a(i+1) / b(i);
        cost = cost + 3; % 1 add/sub 2 mul/div 
        
        a(i+1) = -a(i+1)/b(i);
        cost = cost + 2; % 2 mu;/div 
        
        b(i+1) = b(i+1)+a(i+1)*c(i);
        cost = cost + 2; % 1 add/sub 1 mul/div
        
        d(i+1) = d(i+1)+a(i+1)*d(i);
        cost = cost + 2; % 1 add/sub 1 mul/div
        
        NRow(i+1) = NRow(i+1) - c(i) * NRow(i) / b(i);
        cost = cost + 3; % 1 add/sub 2 mul/div
        
        b(N) = b(N) - NCol(i) * NRow(i) / b(i);
        cost = cost + 3; % 1 add/sub 2 mul/div
        
        d(N) = d(N) - d(i) * NRow(i) / b(i);
        cost = cost + 3; % 1 add/sub 2 mul/div
    end
    b(N) = b(N) - NCol(N-1) * NRow(N-1) / b(N-1);
    cost = cost + 3; % 1 add/sub 2 mul/div
    
    d(N) = d(N) - d(N-1) * NRow(N-1) / b(N-1);
    cost = cost + 3; % 1 add/sub 2 mul/div
    
% Backwards TDMA
    d(N) = d(N)/b(N);
    cost = cost + 1; % 1 mul/div
    
    d(N-1) = (d(N-1)-NCol(N-1)*d(N))/b(N-1);
    cost = cost + 3; % 1 add/sub 2 mul/div
    
    for i=N-2:-1:1
        d(i) = (d(i)-c(i)*d(i+1)-NCol(i)*d(N))/b(i);
        cost = cost + 5; % 2 add/sub 3 mul/div
    end
    
    sol = d;
    comp_cost = cost;
end