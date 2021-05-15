function ke = framestiff(E, A, Ie, x)

% This program calculcates the elemental stiffness matrix for frame
% element.
% Written by Nandha Kumar S 15th May 2021

lx = x(3)  - x(1);
ly = x(4)  - x(2);
le = sqrt (lx^2 + ly^2);

l = lx/le; m = ly/le;
lm = l*m;

fac1 = E*Ie/(le^3);
Ke_beam = [12,   6*le,  -12,   6*le;
      6*le, 4*le^2, -6*le, 2*le^2;
      -12,  -6*le,  12,  -6*le;
      6*le,  2*le^2,  -6*le, 4*le^2]*fac1;

Ke_bar = (E*A/le)*[1 -1; -1 1];

% Overall elemental stiffness Matrix for the frame element
kep = zeros(6,6);
kep(1:3:4, 1:3:4) = Ke_bar;
kep([2:3, 5:6]) = Ke_beam;

Q = zeros(6,6);
Q([1:2, 1:2]) = [l, m; m -1];
Q([4:5, 4:5]) = [l, m; m -1];
Q(3,3) = 1;
Q(6,6) = 1;

ke = Q' * kep * Q;
end
