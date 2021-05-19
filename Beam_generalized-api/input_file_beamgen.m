% Generalised Beam program 

% Stiffness matrix calculation- Gauss Quadrature 
% Element Connectivity and Nodal coordinate array given

% Input Parameters and Data

E = 200e9; Ie = 30*10^-6; q0= 2000; L =3;
% Gauss points and weights
xivec = [-0.774597, 0, 0.774597];
wvec  = [ 5/9, 8/9, 5/9];
% xi1 = -0.774597; w1 = 5/9;
% xi2 =  0.774597; w2 = 5/9;
% xi3 = 0;         w3 = 8/9;
ngauss = length(xivec);
% Coordinates for elements
nele = 3;   % Number of elements
nnode = nele+1; % Number of nodes

% xmat = [0 1;
%         1 2;
%         2 3];

coord = [ 1, 0, 0;  % First columnn is node numbers
          2, 1, 0;  % Second column is coordinate
          3, 2, 0;
          4, 3, 0];

connect = [1, 1,2;  % First column is element number
           2, 2,3;  % Second and third column are nodes in sequence
           3, 3,4]; % for that element

% Boundary condition Data
% First column = node no. global, Second col =  local node no. 
% Third col = value

BC_data = [1, 1, 0;
           1, 2, 0];
       
P_load = [4, 50000]; % First column - Node number, 
                     % second column = Load value
P_moment = [4, 15000];