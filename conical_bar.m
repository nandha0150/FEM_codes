%Conical bar program
clear all
clc

rho=2700;
b=9.81;
d1 = 20/1000;
d2 = 50/1000;
E=70e9;
P1 = 20;
P2= 30;
L =0.8;


xi1= -0.774597;  w1= 5/9;
xi2=  0.774597;  w2= 5/9;
xi3=  0;         w3= 8/9;

%Element 1

xvec = [0; 0.1; 0.2];
xvec1 = xvec;
K1 = stiff_bar(xi1,xvec,d1,d2,L,E)*w1 + stiff_bar(xi2,xvec,d1,d2,L,E)*w2 + ...
    stiff_bar(xi3,xvec,d1,d2,L,E)*w3;
F1 = loadvec_bar(xi1, xvec, d1, d2, L, rho,b)*w1 + loadvec_bar(xi2, xvec, d1, d2, L, rho,b)*w2+...
    loadvec_bar(xi3, xvec, d1, d2, L, rho,b)*w3;


%Element 2

xvec = [ 0.2; 0.35; 0.5];
xvec2 = xvec;
K2 = stiff_bar(xi1,xvec,d1,d2,L,E)*w1 + stiff_bar(xi2,xvec,d1,d2,L,E)*w2 + ...
    stiff_bar(xi3,xvec,d1,d2,L,E)*w3;
F2 = loadvec_bar(xi1, xvec, d1, d2, L, rho,b)*w1 + loadvec_bar(xi2, xvec, d1, d2, L, rho,b)*w2+...
    loadvec_bar(xi3, xvec, d1, d2, L, rho,b)*w3;

%Element 3

xvec = [0.5; 0.65; 0.8];
xvec3 = xvec;
K3 = stiff_bar(xi1,xvec,d1,d2,L,E)*w1 + stiff_bar(xi2,xvec,d1,d2,L,E)*w2 + ...
    stiff_bar(xi3,xvec,d1,d2,L,E)*w3;
F3 = loadvec_bar(xi1, xvec, d1, d2, L, rho,b)*w1 + loadvec_bar(xi2, xvec, d1, d2, L, rho,b)*w2+...
    loadvec_bar(xi3, xvec, d1, d2, L, rho,b)*w3;

%Assembly of matrices

K = zeros(7,7);
K(1:3, 1:3) = K1;
K(3:5, 3:5) = K(3:5, 3:5) + K2;
K(5:7, 5:7) = K(5:7, 5:7) + K3;

F = zeros(7,1);
F(1:3) = F1;
F(3:5) = F(3:5) + F2;
F(5:7) = F(5:7) + F3;


%External Load

F(3) = F(3) + 20;
F(5) = F(5) - 30;


%Imposing BC

K_reduce = K(2:7,2:7);
F_reduce = F(2:7);

%Solution

uvec = K_reduce\F_reduce;

unode = [0; uvec];
xnode= [0; 0.1; 0.2; 0.35; 0.5; 0.65; 0.8];

%Calculation of reactions

F_reac = K*unode;
% Post processing for FEM results

xi = [-1:0.25:1]';

N = [ -xi.*(1- xi)/2, (1-xi.^2), xi.*(1+xi)/2 ];

xn1 = N*xnode(1:3);
un1 = N*unode(1:3);
xn2 = N*xnode(3:5);
un2 = N*unode(3:5);
xn3 = N*xnode(5:7);
un3 = N*unode(5:7);

xn = [xn1; xn2; xn3];
un = [un1; un2; un3];

%%
% Analytical solution

delx = 0.0025;
tht = (d2- d1)/L;
g = 9.81;

% Deflection in portion 1 due to point load

P1 = -10;
x = 0:delx:0.2;
up1 = 4*P1*x/pi/E/d1./(d1 +tht*x);
x1= x;
uA = up1(end);

% Deflection in portion 2 due to point load

P2 = -30;
x = 0.2:delx:0.5;
xbar= x-0.2;
d1bar = d1 + tht*0.2; % new top surface
up2 = 4*P2*xbar/pi/E/d1bar./(d1bar +tht*xbar);

x2 = x;
up2 = uA + up2;
uB = up2(end); %deflection of point B

% third element 
% No point loads, 
% so each point is deflected by
% a constant amount. 

x = 0.5:delx:0.8;
x3 = x;
up3 = uB*ones(1,size(x,2));

%Combining the point loads
xana = [x1, x2, x3];
up = [up1, up2, up3];

%Analytical deflection due to Body force

dx = d1 + tht*xana;
ug = (rho*g*d2^3/3/E/tht/d1)*xana./dx -rho*g*d1*xana/3/E/tht - rho*g*xana.^2/6/E;

%Total Analytical deflection

uana = up + ug;

%Visualisation 
h = figure(1)
plot(xana, uana, 'b-', xn, un, 'ro', 'linewidth',2,'MarkerEdgeColor',...
    'k', 'MarkerFaceColor', 'r','MarkerSize',8);
hold on
set(gcf, 'Position', get(0, 'Screensize'));
set(gca, 'Fontsize', 12, 'Fontweight', 'demi');
set(gcf, 'defaultTextInterpreter', 'latex');
xlabel('x', 'fontsize', 18);
ylabel('u', 'fontsize', 18);
legend('Analytical', 'FEM')
grid on 
hold on



