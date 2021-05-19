function [xnume, unume] = postprocessing_beam_gen(nele, coord, connect, un, xi)

xnume=[]; unume=[];
for el=1:nele
    nd1 = connect(el,2);
    nd2 = connect(el,3);
    x_n = [coord(nd1,2), coord(nd2, 2)];
    vec = [2*nd1-1, 2*nd1, 2*nd2-1, 2*nd2];
    u_n = un(vec);
    
    le = x_n(2) - x_n(1);
    Nx= [ (1- xi)/2, (1+xi)/2 ];
    N1=  (2 - 3*xi + xi.^3)./4;
    N2=  (1 - xi - xi.^2 +xi.^3)/4;
    N3=  (2 + 3*xi - xi.^3)./4;
    N4=  (-1 - xi + xi.^2 +xi.^3)/4;
    Nu=  [N1, le*N2/2, N3, le*N4/2];
    xnume = [ xnume; Nx*x_n'];
    unume = [ unume; Nu*u_n];
  
end
    