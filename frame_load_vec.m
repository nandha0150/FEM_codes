function f = frame_load_vec(qa, qt,x)

% This program calculcates the elemental load vector for frame element for
% uniform q.
% Written by Nandha Kumar S 15th May 2021

% Always calculate the distributed load in local frame and then premultiply
% with Q to get global f vec. 

lx = x(3)  - x(1);
ly = x(4)  - x(2);
le = sqrt (lx^2 + ly^2);

l = lx/le; m = ly/le;
lm = l*m;

fe_beam = [qt*le/2; qt*le^2/12; qt*le/2; -qt*le^2/12];
fe_bar = (qa*le)*[0.5;0.5];

fep = zeros(6,1);
fep([1;4]) = fe_bar;
fep([2;3;5;6]) = fe_beam;


Q = zeros(6,6);
Q([1:2, 1:2]) = [l, m; m -1];
Q([4:5, 4:5]) = [l, m; m -1];
Q(3,3) = 1;
Q(6,6) = 1;

% Q - Transformation matrix
% Global f vector 
f = Q' * fep;

end