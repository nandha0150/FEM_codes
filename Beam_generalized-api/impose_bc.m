function [K,F]= impose_bc(K,F, BC_data)

% Written by Nandha Kumar S
% on 19th May 2021
% krrish0150@gmail.com

%Input 
    % K - Stiffness matrix before applying BC
    % F - Force vector before applying BC with point loads, moments applied
%Output
    % K - Stiffness matrix after applying BC
    % F - Force vector after applying BC with point loads, moments applied
    
% Modification of Load vector
% Based on BC application
supp_dof = [];
for ii = 1:size(BC_data,1)
    nd = BC_data(ii,1);
    dof = BC_data(ii,2);
    val = BC_data(ii,3);
    GDOF = 2*(nd-1) +dof;
    supp_dof = [supp_dof, GDOF];
    if val ~= 0
        for jj=1:2*nnode
            F(jj) = F(jj) - K(jj, GDOF)*val;
        end
    end
end

% Imposing BC on stiffness Matrix
% Elimination of row and column
% of Stiffness and Load vector

supp_dof = sort(supp_dof);
for ii = 1:size(supp_dof,2)
    dof = supp_dof(ii);
    if dof == 1
        % Removing first row
        K = K(dof+1:end, dof+1:end);
        F = F(dof+1:end);
    elseif dof == 2*nnode
        % Removing last row
        K = K(1:dof-1, 1:dof-1)
        F = F(1:dof-1);
    else
        % Removing elsewhere row
        K = K([1:dof-1, dof+1:end], [1:dof-1, dof+1:end]);
        F = F([1:dof-1, dof+1:end]);
    end
    if(ii ~= size(supp_dof,2))
        supp_dof(ii+1:end) = supp_dof(ii+1:end) -1;
        % We do this because after removing the row and column of 
        % the next removable elements are present in rows and cols 
        % one less than their original. 
    end
end

end
