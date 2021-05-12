function fe = beamload_vec2(q0, L, x)

le = x(2) - x(1);
fe1 = [q0*le/2; q0*le^2/12; q0*le/2; -q0*le^2/12];

fe2 = (q0*le/(2*L))*[7/10, 3/10; le/10, le/15; 3/10, 7/10; -le/15, -le/10]*x';

fe = fe1 - fe2;

end
