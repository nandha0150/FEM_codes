function f = beamload_vec(q,x)

le = x(2) - x(1);

f = [   q*le/2;
        q*le^2/12;
        q*le/2;
        -q*le^2/12];

end
