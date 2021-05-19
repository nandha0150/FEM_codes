function F = point_ld_moment(F, P_load, P_moment)

% This function applies the point load and moment at the corresponding node
% of the Load vector
% Written by Nandha Kumar S 
% krrish0150@gmail.com


for ii=1:size(P_load,1)
    nd = P_load(ii,1);
    F0 = P_load(ii,2);
    F(2*nd-1) = F(2*nd -1) + F0;
end

for ii = 1:size(P_moment, 1)
    nd = P_moment(ii,1);
    M0 = P_moment(ii,2);
    F(2*nd) = F(2*nd) + M0;
end