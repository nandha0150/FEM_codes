function ke = beam_elestiff_gen(xi, E, Ie,x)

% This program evaluates the integrand(Stiffness Matrix)
% at particular gauss point
% Written by Nandha Kumar S on 17th May 2021

le = x(2) - x(1);
B1 = 3*xi/2;
B2 = le * (3*xi - 1)/4;
B3 = -3*xi/2;
B4 = le* (3*xi + 1)/4;

B = (4/le^2)*[B1, B2, B3, B4];

ke = (E*Ie*le/2)*(B'*B);

end
