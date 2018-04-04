U = 10;
w = 0.50;
c = 1;
b = c/2;

t = 0:0.001:2.5;

tau = U*t./b;

Cl = 2*pi*(w/U)*(1 - 0.5*exp(-0.13*tau) - 0.5*exp(-tau));
