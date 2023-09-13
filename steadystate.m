function X = steadystate(guess)

global beta delta rho alpha eta theta1 theta2 gamma phic d2 d1 d0 erowcoef;

k = guess(1);
e = guess(2);
mu = guess(3);

i = delta*k;
x = erowcoef*e/(1-eta);
y = (1-d2*(x)^2-d1*(x)-d0)*k^alpha;
z = theta1*mu^theta2*y;
c = y - i - z;

X(1) = e - (1-mu)*y^(1-gamma);
X(2) = -1 + beta*((1-d2*(x)^2-d1*(x)-d0)*alpha*k^(alpha-1)*...
    (1-theta1*mu^theta2-y*theta1*theta2*mu^(theta2-1)*e/y^(2-2*gamma)*...
    (1-gamma)*y^(-gamma))+1-delta);
X(3) = (-(1-theta1*mu^theta2)*(2*d2*(x)+d1)*k^alpha + ...
    y*theta1*theta2*mu^(theta2-1)*(y^(1-gamma)+e*(1-gamma)*y^(-gamma)*...
    (2*d2*(x)+d1)*k^alpha)/y^(2-2*gamma)) - ...
    beta*theta1*theta2*mu^(theta2-1)*eta/y^(1-gamma);