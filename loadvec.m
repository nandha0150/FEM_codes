function f = loadvec(xi,xvec)
%-------------------------------------------------------------%
% This function calculates the integrand of element load vector 
% for one element (xvec) at one gauss point (xi)
%-------------------------------------------------------------%
N = [ -xi*(1- xi)/2, (1-xi^2), xi*(1+xi)/2 ];
B = [(xi -0.5), -2*xi, (xi+ 0.5)];

J = B*xvec;
xe = N* xvec;

f = (N'*xe^2 /4)* J;

end
