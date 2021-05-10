% This program is to solve a sample truss problem 
% For problem data and description refer notes.
% Nandha Kumar S 

clear all 
clc

E = 29.5e6;
A= 1;

% Element 1
x1 = [0 0 40 0];
k1 = elestiff_truss(E,A,x1);

% Element 2
x2 = [40 0 40 30];
k2 = elestiff_truss(E,A,x2);


% Element 3
x3 = [40 30 0 30];
k3 = elestiff_truss(E,A,x3);


% Element 4
x4 = [40 30 0 0];
k4 = elestiff_truss(E,A,x4);

% Assembly

K = zeros(8,8);
F = zeros(8,1);

K(1:4,1:4) =  k1;
K(3:6, 3:6) = K(3:6, 3:6) + k2(1:4,1:4);
K(5:8, 5:8) = K(5:8, 5:8) + k3(1:4,1:4);
K([ 5 6 1 2], [5 6 1 2]) = K([ 5 6 1 2], [5 6 1 2]) + k4(1:4,1:4);

F(3) = 20000;
F(6) = -25000;

% Imposition of BC

Kreduce = K([3,5:6], [3,5:6]);
Freduce = F([3,5:6]);

%Finding solution
ureduce = Kreduce\Freduce;

%Finding reaction forces
un = zeros(8,1);
un([3,5:6]) = ureduce;

Fr= K*un;


