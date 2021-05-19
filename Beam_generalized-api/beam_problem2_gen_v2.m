% Generalized function for solving Beam Problem 
% Written by Nandha Kumar S on 19th May 2021
% krrish0150@gmail.com

% Required functions are 
    % input_file_beamgen
    % stiffload_beamgen
    % 
% Calling input data 
input_file_beamgen;

% Calculation of Stiffness Matrix and Load vector

[K,F] = stiffload_beamgen(nele, ngauss,coord,connect,xivec,wvec,E,Ie, q0,L);

% Applying the point load and moment

F = point_ld_moment(F, P_load, P_moment);

% Imposing BC

Kg = K;
Fg = F;

[K, F] = impose_bc(K, F, BC_data);

% Solving the system
ureduce = K \ F;

un = [0;0; ureduce];
Fr = Kg* un;

% Post processing

xi = [-1:0.2:1]';
[xnume, unume] = postprocessing_beam_gen(nele, coord, connect, un, xi);

% Analytical solution 

xana = L*[0:0.01:1];
xbar = xana /L;
fac1 = q0*L^4/120/E/Ie;
fac2 = P_load(2)*L^3/6/E/Ie;
fac3 = P_moment(2)*L^2/2/E/Ie;

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
