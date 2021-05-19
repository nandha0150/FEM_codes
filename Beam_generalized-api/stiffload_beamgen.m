% Calculation of element stiffness matrix and load vector

function [K,F] = stiffload_beamgen(nele, ngauss,coord,connect,xivec,wvec,E,Ie, q0,L)
nnode = nele+1;
K = zeros(2*nnode, 2*nnode);
F = zeros(2*nnode,1);

%Calculation of element stiffness matrix

for el = 1:nele
    kele = zeros(4,4);
    fele = zeros(4,1);
%     x= xmat(el,:);
    nd1 = connect(el,2);
    nd2 = connect(el,3);
    x = [coord(nd1,2), coord(nd2,2)];
    vec = [2*nd1-1, 2*nd1, 2*nd2-1, 2*nd2];
    
    for gp = 1:ngauss
        
        xi = xivec(gp); w= wvec(gp);
        kele(1:4,1:4) = kele(1:4,1:4) + beam_elestiff_gen(xi,E,Ie,x)*w;
        fele(1:4) = fele(1:4) + beam_eleload_gen(xi,q0,L,x)*w;
    end
    
     for  ii= 1:4
        for jj =1:4
            K(vec(ii), vec(jj)) = K(vec(ii), vec(jj)) + kele(ii, jj);
        end
            F(vec(ii)) = F(vec(ii)) + fele(ii);
    end
    
end
