% Beam problem Example 2
% Main Progam
% Written by Nandha Kumar S on 12th May 2021

% The description of the problem and details, refer text.
% Required functions - beamstiff.m, beamload_vec2.m

clear all 
clc

E = 200e9;
L=3;
I= 30*10^(-6);
q0=20000;
M0=15000;
F0=50000;

% Element 1
x = [0, 1];
k1 = beamstiff(E,I,x);
f1 = beamload_vec2(q0, L, x);

% Element 2
x = [1, 2];
k2 = beamstiff(E,I,x);
f2 = beamload_vec2(q0, L, x);

% Element 3
x = [2, 3];
k3 = beamstiff(E,I,x);
f3 = beamload_vec2(q0, L, x);

K = zeros(8,8);
F = zeros(8,1);
uvec = zeros(8,1);

% Assembling the system
K(1:4, 1:4) = k1;
K(3:6, 3:6) = K(3:6, 3:6) + k2;
K(5:8, 5:8) = K(5:8, 5:8) + k3;

F(1:4) = f1;
F(3:6) = F(3:6) + f2;
F(5:8) = F(5:8) + f3;
F(7) = F(7) + F0;
F(8) = F(8) + M0;

% Imposing BC
K_new = K(3:8, 3:8);
F_new = F(3:8);

%Solving
u_new = K_new \ F_new;

%Reaction force calculation 
uvec(3:8) = u_new;
Fr = K * uvec;

% FE solution Interpolation
xn = [0,1,2,3];
xnume = [];
unume = [];

for iel =1:3
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

xana = L*[0:0.01:1];
xbar = xana /L;
fac1 = q0*L^4/120/E/I;
fac2 = F0*L^3/6/E/I;
fac3 = M0*L^2/2/E/I;

uana = fac1 * (10*xbar.^2 - 10*xbar.^3 + 5*xbar.^4 -xbar.^5) + ...
    fac2*(-xbar.^3 + 3*xbar.^2) + fac3*xbar.^2;

% Plotting

plot(xana, uana, 'r-', xnume, unume, 'bo')
legend('Analytical', 'FE solution');
xlabel('x axis');
ylabel('Vertical displacement');
title('Beam Problem (2) solution plot')
grid on 
hold on