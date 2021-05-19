function f = beam_eleload_gen(xi,q0,L,x)
% This program is generalized calculation of beam element
% load vector
% Using Gauss Quadrature
% Written by Nandha Kumar S on 17th May 2021

le = x(2) - x(1);

fac1 = q0*le/(2*L);
Nx= [ (1- xi)/2, (1+xi)/2 ];

N1 =  (2 - 3*xi + xi.^3)./4;
N2 =  (1 - xi - xi.^2 +xi.^3)/4;
N3 =  (2 + 3*xi - xi.^3)./4;
N4 =  (-1 - xi + xi.^2 +xi.^3)/4;
N  =  [N1, le*N2/2, N3, le*N4/2];

xe = Nx*x';
q = q0*(1- xe/L);
f = N'*q*le/2;

end

