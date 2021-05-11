function Ke = beamstiff(E, Ie, x)

le = x(2) - x(1);
fac = E*Ie/(le^3);

Ke = [12,   6*le,  -12,   6*le;
      6*le, 4*le^2, -6*le, 2*le^2;
      -12,  -6*le,  12,  -6*le;
      6*le,  2*le^2,  -6*le, 4*le^2]*fac;

end
