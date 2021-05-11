% Beam problem
% Main Progam
% Written by Nandha Kumar S on 11th May 2021

% The description of the problem and details, refer text.
% Required functions - beamstiff.m, beamload_vec.m

clear all 
clc

E = 210e9;
l1= 2; l2=2; l3=1; l4=1;
I=0.5*10^(-6);
q=2000;
M=3000;
Fp=10000;

% Element 1
x = [0, 2];
k1 = beamstiff(E,I,x);


% Element 2
x = [2, 4];
k2 = beamstiff(E,I,x);
f2=  beamload_vec(q,x);

% Element 3
x = [4, 5];
k3 = beamstiff(E,I,x);


% Element 4
x = [5, 6];
k4 = beamstiff(E,I,x);

%Assembly
K = zeros(10,10);
F = zeros(10,1);
uvec = zeros(10,1);

K(1:4, 1:4) = k1;
K(3:6, 3:6) = K(3:6, 3:6) + k2;
K(5:8, 5:8) = K(5:8, 5:8) + k3;
K(7:10, 7:10) = K(7:10, 7:10) + k4;

F(3:6) = f2; % Distributed load
F(3) = F(3) + Fp; % Global dof of, Point load associated with node 2
F(8) = F(8) + M;  % Concentrated moment at node 4

%Imposing BC
Kreduce = K([2:8,10], [2:8,10]);
Freduce = F([2:8,10]);

% Solution

ureduce = Kreduce \ Freduce;

uvec([2:8,10]) = ureduce;

% Reaction forces

Fr = K * uvec;

% FE solution Interpolation
xn = [0,2,4,5,6];
xnume = [];
unume = [];

for iel =1:4
    x_n = xn(iel:iel+1);
    u_n = uvec((iel-1)*2+1:(iel+1)*2);
    le = x_n(2) - x_n(1);
    xi = [-1:0.2:1]';
    Nx= [ (1- xi)/2, (1+xi)/2 ];
    N1=  (2 - 3*xi + xi.^3)./4;
    N2=  (1 - xi - xi.^2 +xi.^3)/4;
    N3=  (2 + 3*xi - xi.^3)./4;
    N4=  (-1 - xi + xi.^2 +xi.^3)/4;
    Nu=  [N1, le*N2/2, N3, le*N4/2];
    xnume = [ xnume; Nx*x_n'];
    unume = [ unume; Nu*u_n];
end


% Analytical solution
l =6; a=2; b=4; c=5;

xa_1 = [0:0.01:1]*a;
ya_1 = (28.1589*xa_1 - 8.17*xa_1.^3/6) /(E*I);

xa_2 = a + [0:0.01:1]*(b-a);
ya_2 = (28.1589*xa_2 - 8.17*xa_2.^3/6 + 5*(xa_2 -a).^3/3 + (xa_2 -a).^4/12)/(E*I);

xa_3 = b + [0:0.01:1]*(c-b);
ya_3 = (28.1589*xa_3 - 8.17*xa_3.^3/6 + 5*(xa_3 -a).^3/3 + (xa_3 -a).^4/12 - (xa_3 -b).^4/12 )/(E*I);

xa_4 = c + [0:0.01:1]*(l-c);
ya_4 = (28.1589*xa_4 - 8.17*xa_4.^3/6 + 5*(xa_4 -a).^3/3 + (xa_4 -a).^4/12 - (xa_4 -b).^4/12 - 3*(xa_4-c).^2/2)/(E*I);

x_ana = [xa_1, xa_2, xa_3, xa_4];
y_ana = [ya_1, ya_2, ya_3, ya_4]*1000;

%Plotting
plot(x_ana, y_ana, 'r-', xnume, unume, 'bo-');
legend('Analytical', 'FE solution');
xlabel('x axis');
ylabel('Vertical displacement');
title('Beam Problem solution plot')
grid on 
hold on