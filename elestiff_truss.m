function Ke = elestiff_truss(E,A,x)
    
% This program calculcates the elemental stiffness matrix for single truss
% element.
% Written by Nandha Kumar S 10th May 2021

lx = (x(3) - x(1));
ly = (x(4) - x(2));

L = sqrt(lx^2 + ly^2);
l = lx/L; m = ly/L;

fac = A*E /L;

Ke = fac*[ l^2, l*m, -l^2, -l*m;
           l*m, m^2, -l*m, -m^2;
           -l^2, -l*m, l^2, l*m;
           -l*m, -m^2, l*m, m^2];
end
