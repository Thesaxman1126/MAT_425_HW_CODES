clear
%% Problem 3
%% Analysis on lowest degree polynomial construction
% Lagrange: 2
% Vandermond: 2
%% Initial Conditions/Data
x = [1,0,-1];
y = [0,-1,0];
xx = min(x)-.1:0.1:max(x)+.1; % used to make smooth curve
%% Lagrange interpolation method
% initialize coeff and final polynomial
a_i = 0;
pL = 0;
% solve for coeff
for i=1:length(x)
    L = 1;
    for j=1:length(x)
        if j~=i
            L = conv(L,poly(x(j))/(x(i)-x(j)));
        end
    end
    term = L*y(i);
    a_i = a_i + term;
end
% construct final polynomial
for m = 1:length(x)
    pL = a_i(length(x)-m+1)*xx.^(m-1) + pL;
end
% plot polynomial
figure(1)
plot(xx,pL,'b', x, y, "ko", 'MarkerFaceColor', 'k')
title('Lagrangian Interpolation Polynomial')
%% Vandermonde interpolation method
% solve Vandermonde Coeff
xx = min(x)-.1:0.1:max(x)+.1;
V = fliplr(vander(x));
bV = linsolve(V,y');
% Construct Polynomial
k = 0;
pV = 0;
while k < length(x)
    pV = pV + bV(k+1)*xx.^(k);
    k = k+1;
end
% Plot
figure(2)
plot(xx, pV, x, y, "ko", 'MarkerFaceColor', 'k')
title('Vandermonde Interpolation Polynomial')
%% Both on same plot
figure(3)
plot(xx, pL, "r", xx, pV, "b--", x, y, "ko", 'MarkerFaceColor', 'k');
title('Lagrangian and Vandermonde Interpolation Polynomials')
legend('Lagrangian','Vandermonde',"Raw Data")