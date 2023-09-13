function X = steadystate_tax(guess)   

global beta delta rho alpha eta theta1 theta2 gamma phic d2 d1 d0;


k = guess(1);
e = guess(2);
%lambda = guess(3);

i = delta*k;
x = 4*e/(1-eta);
y = (1-d2*(x)^2-d1*(x)-d0)*k^alpha;
mu = 1-e/y^(1-gamma);

z = theta1*mu^theta2*y;
c = y - i - z;
tau = theta1*theta2*mu^(theta2-1)*y^gamma;
r = y*alpha*k^(-1)*(1-tau*(1-mu)*(1-gamma)*y^(-gamma)-theta1*mu^theta2);

%derivatives of firm demand equations
zy = (theta2*(1-gamma)-1)/(theta2-1)*(z/y);
ztau = theta2/(theta2-1)*(z/tau);
etau = -1/(theta2-1)*y^(1-gamma)*mu/tau;
ey = (1-gamma)*y^(-gamma) + (gamma/(theta2-1)+gamma-1)*mu*y^(-gamma);
rtau = -alpha*(1-gamma)*y^(1-gamma)/k...
    + alpha*(1-gamma)*(1+1/(theta2-1))*mu*y^(1-gamma)/k...
    - alpha*theta1*theta2/(theta2-1)*y/k*mu^theta2/tau;
ry = alpha/k - alpha*(1-gamma)^2*tau*y^(-gamma)/k...
    + alpha*(1-gamma)*(1-gamma-gamma/(theta2-1))*tau*mu*y^(-gamma)/k...
    - alpha*theta1*(1-theta2*gamma/(theta2-1))*mu^theta2/k;
rk = -alpha*y/k^2 + alpha*(1-gamma)*tau*y^(1-gamma)/k^2 ...
    - alpha*(1-gamma)*tau*mu*y^(1-gamma)/k^2 ...
    + alpha*theta1*y/k^2*mu^theta2;

lambda = (-ztau/etau*(ey+(1-eta*beta)/(k^alpha*(2*d2*x+d1)))-(1-zy))...
    *(phic/c*(1-zy)*(1-1/beta)+ry-(phic/c*ztau*(1/beta-1)+rtau)/etau...
    *(ey+(1-eta*beta)/(k^alpha*(2*d2*x+d1))))^(-1);
zeta = c^(-phic)*(-ztau + lambda*(phic*c^(-1)*ztau*(1/beta-1)+rtau))/etau; 
omega = -zeta*(1-eta*beta)/k^alpha/(2*d2*x+d1);

X(1) = r - 1/beta - delta + 1;

%X(2) = c^(-phic)*(1-zy) + lambda*(phic*c^(-phic-1)*(1-zy)*(1-1/beta)+c^(-phic)*ry)...
%    -zeta*ey + omega;

X(2) = c^(-phic)*(beta*(1-delta)-1)...
    + lambda*(phic*c^(-phic-1)*((1-delta)*(beta-1)+1/beta-1) + beta*c^(-phic)*rk)...
    - alpha*beta*omega*y/k;
