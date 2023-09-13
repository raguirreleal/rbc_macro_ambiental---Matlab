%Program to solve pollution/RBC model with AIM
%New specification of damages to output, not utility
%Decentralized model 
%Asymmetric info - govt does not observe the state, just yt-1
%2/23/2010

clear;

%First define the parameter values
global beta delta rho alpha eta theta1 theta2 gamma phic d2 d1 d0;

% beta = 0.98267; %discount rate
% delta = 0.025; %capital depreciation
% rho = .95; %persistence of TFP shock
% sd = 0.007; %standard deviation of shock
% alpha = 0.36; %curvature of production function
% 
% eta = 0.9979; %pollution depreciation 
% theta1 = .05607; %abatement cost equation parameters, from Nordhaus
% theta2 = 2.8;
% gamma = 1-.696; %1 - elasticity of emissions with respect to output
% phic = 2; %CRRA for consumption
% d2 = 1.4647*10^(-8); %damage function parameters, from Nordhaus
% d1 = -6.6722*10^(-6);
% d0 = 1.395*10^(-3);
% damage_scale = 5.3024; %To scale the pollution levels and get the damage function correct
% d2 = d2/damage_scale^2;
% d1 = d1/damage_scale;
% d = 0.6580; %The government's policy variable

% calibração brasil 
 
beta   = 0.98;      % - 
delta  = 0.025;     % - 
rho    = 0.95;      % - 
sd     = 0.0095;    % - 
alpha  = 0.40;      % - 

eta    = 0.9979;       % de Reilly (1992) e Heutel - pollution depreciation 
theta1 = 0.04183;      % - Nordhaus - abatement cost equation parameters, from Nordhaus 
                       % ...arquivo 'RICE_042510'; planilha 'LatAm'; célula 'C31'
theta2 = 2.8;          % este... da planilha 'Parameters'; célula 'C63' (manteve)
gamma  = 1 - 1.07024;  % via primeira diferença das séries 
phic   = 2;            % - 
d2     = 9.26191*10^(-9);   % - %damage function parameters, from Nordhaus
d1     = -2.16474*10^(-6);   % - 
d0     = -0.0029736;        % - 
erowcoef  = 80;          % - pelo CDIAC (2010) - coefficient relating how rest-of-world emissions compare to domestic emissions (era 4)
                         % CDIAC: em 2010 total mundo = 9.167.000; Brasil = 114.468; US = 1.481.608 (thousand metric tons of carbon)
                         % resultado para Brazil é erowcoef  = 80, mas explode antes
dmg_scl  = 5.3024;       % damage scale - To scale the pollution levels and get the damage function correct (manteve)
d2       = d2/dmg_scl^2; 
d1       = d1/dmg_scl;   

d = 1.4; %The government's policy variable

imprespsize = .01;

%Solve for steady state solution based on first order condition
a_ss = 1;

%guess at the steady state values
%Guess values are from the old model's solution
k_g = 36;
e_g = 4;

guess = [k_g,e_g]; %guess for vector of steady state values
options=optimset('Display','iter','MaxFunEvals',5000);
ss_sol = fsolve(@steadystate_tax,guess,options);
clear options guess k_g e_g;

k_ss = ss_sol(1);
e_ss = ss_sol(2);

i_ss = delta*k_ss;
x_ss = 4*e_ss/(1-eta);
y_ss = (1-d2*(x_ss)^2-d1*(x_ss)-d0)*k_ss^alpha;
mu_ss = 1-e_ss/y_ss^(1-gamma);
z_ss = theta1*mu_ss^theta2*y_ss;
c_ss = y_ss - i_ss - z_ss;
tau_ss = theta1*theta2*mu_ss^(theta2-1)*y_ss^gamma;
r_ss = y_ss*alpha*k_ss^(-1)*(1-tau_ss*(1-mu_ss)*(1-gamma)*y_ss^(-gamma)-theta1*mu_ss^theta2);

%derivatives of firm demand equations
zy = (theta2*(1-gamma)-1)/(theta2-1)*(z_ss/y_ss);
ztau = theta2/(theta2-1)*(z_ss/tau_ss);
etau = -1/(theta2-1)*y_ss^(1-gamma)*mu_ss/tau_ss;
ey = (1-gamma)*y_ss^(-gamma) + (gamma/(theta2-1)+gamma-1)*mu_ss*y_ss^(-gamma);
rtau = -alpha*(1-gamma)*y_ss^(1-gamma)/k_ss...
    + alpha*(1-gamma)*(1+1/(theta2-1))*mu_ss*y_ss^(1-gamma)/k_ss...
    - alpha*theta1*theta2/(theta2-1)*y_ss/k_ss*mu_ss^theta2/tau_ss;
ry = alpha/k_ss - alpha*(1-gamma)^2*tau_ss*y_ss^(-gamma)/k_ss...
    + alpha*(1-gamma)*(1-gamma-gamma/(theta2-1))*tau_ss*mu_ss*y_ss^(-gamma)/k_ss...
    - alpha*theta1*(1-theta2*gamma/(theta2-1))*mu_ss^theta2/k_ss;
rk = -alpha*y_ss/k_ss^2 + alpha*(1-gamma)*tau_ss*y_ss^(1-gamma)/k_ss^2 ...
    - alpha*(1-gamma)*tau_ss*mu_ss*y_ss^(1-gamma)/k_ss^2 ...
    + alpha*theta1*y_ss/k_ss^2*mu_ss^theta2;

lambda_ss = (-ztau/etau*(ey+(1-eta*beta)/(k_ss^alpha*(2*d2*x_ss+d1)))-(1-zy))...
    *(phic/c_ss*(1-zy)*(1-1/beta)+ry-(phic/c_ss*ztau*(1/beta-1)+rtau)/etau...
    *(ey+(1-eta*beta)/(k_ss^alpha*(2*d2*x_ss+d1))))^(-1);
zeta_ss = c_ss^(-phic)*(-ztau + lambda_ss*(phic*c_ss^(-1)*ztau*(1/beta-1)+rtau))/etau; 
omega_ss = -zeta_ss*(1-eta*beta)/k_ss^alpha/(2*d2*x_ss+d1);


%Write out the structural coefficient matrix h
%This is for the firm's problem w/o information uncertainty
%It is the same for tax or quota so just keep the tax equations
%The order of the parameters is [c,e,k,x,mu,y,z,i,tau,r,lambda,zeta,omega,a]
neq = 14;
nlag = 1;
nlead = 1;
h = zeros(neq,neq*(nlag+1+nlead));
%Columns 1-14 are variables at t-1
%Columns 15-28 are variables at t
%Columns 28-42 are variables at t+1

%Eqn 1: evolution of shock a
h(1,14) = -rho;
h(1,28) = 1;

%Eqn 2: Budget constraint
h(2,15) = c_ss;
h(2,22) = i_ss;
h(2,21) = z_ss;
h(2,20) = -y_ss;

%Eqn 3: Capital accumulation
h(3,17) = k_ss;
h(3,3) = -k_ss*(1-delta);
h(3,22) = -i_ss;

%Eqn 4: pollution stock
h(4,4) = -eta;
h(4,16) = -(1-eta);
h(4,18) = 1;

%Eqn 5: Emissions
h(5,16) = e_ss;
h(5,20) = y_ss^(1-gamma)*(1-gamma)*(mu_ss-1);
h(5,19) = mu_ss*y_ss^(1-gamma);

%Eqn 6: Abatement
h(6,21) = 1;
h(6,19) = -theta2;
h(6,20) = -1;

%Eqn 7: Production function
h(7,20) = y_ss;
h(7,28) = -(1-d2*x_ss^2-d1*x_ss-d0)*k_ss^alpha;
h(7,3) = -(1-d2*x_ss^2-d1*x_ss-d0)*k_ss^alpha*alpha;
h(7,18) = k_ss^alpha*(2*d2*x_ss^2+d1*x_ss);

%Eqn 8: Consumer's FOC
h(8,15) = phic;
h(8,29) = beta*(r_ss+1-delta)*(-phic);
h(8,38) = beta*r_ss;

%Eqn 9: Firm's FOC wrt k
h(9,20) = alpha*y_ss/k_ss - alpha*(1-gamma)*tau_ss*y_ss^(1-gamma)/k_ss*(1-gamma)...
    + alpha*(1-gamma)*tau_ss*mu_ss*y_ss^(1-gamma)/k_ss*(1-gamma)...
    - alpha*theta1*y_ss/k_ss*mu_ss^theta2;
h(9,3) = alpha*y_ss/k_ss*(-1) - alpha*(1-gamma)*tau_ss*y_ss^(1-gamma)/k_ss*(-1)...
    + alpha*(1-gamma)*tau_ss*mu_ss*y_ss^(1-gamma)/k_ss*(-1)...
    - alpha*theta1*y_ss/k_ss*mu_ss^theta2*(-1);
h(9,23) = -alpha*(1-gamma)*tau_ss*y_ss^(1-gamma)/k_ss...
    + alpha*(1-gamma)*tau_ss*mu_ss*y_ss^(1-gamma)/k_ss;
h(9,19) = alpha*(1-gamma)*tau_ss*mu_ss*y_ss^(1-gamma)/k_ss...
    - alpha*theta1*y_ss/k_ss*mu_ss^theta2*theta2;
h(9,24) = -r_ss;

%Eqn 10: Firm's FOC wrt mu
h(10,23) = 1;
h(10,20) = -gamma;
h(10,19) = -(theta2-1);


%Eqns 11-14 are the planner's FOCs
%These are dropped in the imperfect info case
%Eqn 11: Planner's FOC wrt tau
h(11,15) = theta2/(theta2-1)*(-c_ss^(-phic))*z_ss/tau_ss*(-phic)...
    - phic*theta2/(theta2-1)*lambda_ss*c_ss^(-phic-1)*z_ss/tau_ss*(-phic-1)...
    + lambda_ss*(-phic)*c_ss^(-phic-1)*(-theta2)/(theta2-1)*z_ss/tau_ss*(r_ss+1-delta)*(-phic-1)...
    + lambda_ss*c_ss^(-phic)*(-alpha)*(1-gamma)*y_ss^(1-gamma)/k_ss*(-phic)...
    + lambda_ss*c_ss^(-phic)*alpha*(1-gamma)*(1+1/(theta2-1))*mu_ss*y_ss^(1-gamma)/k_ss*(-phic)...
    + lambda_ss*c_ss^(-phic)*(-alpha)*theta1*theta2/(theta2-1)*y_ss/k_ss*mu_ss^theta2/tau_ss*(-phic);
h(11,21) = theta2/(theta2-1)*(-c_ss^(-phic))*z_ss/tau_ss...
    - phic*theta2/(theta2-1)*lambda_ss*c_ss^(-phic-1)*z_ss/tau_ss...
    + lambda_ss*(-phic)*c_ss^(-phic-1)*(-theta2)/(theta2-1)*z_ss/tau_ss*(r_ss+1-delta);
h(11,23) = theta2/(theta2-1)*(-c_ss^(-phic))*z_ss/tau_ss*(-1)...
    - phic*theta2/(theta2-1)*lambda_ss*c_ss^(-phic-1)*z_ss/tau_ss*(-1)...
    + lambda_ss*(-phic)*c_ss^(-phic-1)*(-theta2)/(theta2-1)*z_ss/tau_ss*(r_ss+1-delta)*(-1)...
    + lambda_ss*c_ss^(-phic)*(-alpha)*theta1*theta2/(theta2-1)*y_ss/k_ss*mu_ss^theta2/tau_ss*(-1)...
    + zeta_ss/(theta2-1)*mu_ss*y_ss^(1-gamma)/tau_ss*(-1);
h(11,25) = -phic*theta2/(theta2-1)*lambda_ss*c_ss^(-phic-1)*z_ss/tau_ss;
h(11,11) = lambda_ss*(-phic)*c_ss^(-phic-1)*(-theta2)/(theta2-1)*z_ss/tau_ss*(r_ss+1-delta)...
    + lambda_ss*c_ss^(-phic)*(-alpha)*(1-gamma)*y_ss^(1-gamma)/k_ss...
    + lambda_ss*c_ss^(-phic)*alpha*(1-gamma)*(1+1/(theta2-1))*mu_ss*y_ss^(1-gamma)/k_ss...
    + lambda_ss*c_ss^(-phic)*(-alpha)*theta1*theta2/(theta2-1)*y_ss/k_ss*mu_ss^theta2/tau_ss;
h(11,24) = lambda_ss*(-phic)*c_ss^(-phic-1)*(-theta2)/(theta2-1)*z_ss/tau_ss*r_ss;
h(11,20) = lambda_ss*c_ss^(-phic)*(-alpha)*(1-gamma)*y_ss^(1-gamma)/k_ss*(1-gamma)...
    + lambda_ss*c_ss^(-phic)*alpha*(1-gamma)*(1+1/(theta2-1))*mu_ss*y_ss^(1-gamma)/k_ss*(1-gamma)...
    + lambda_ss*c_ss^(-phic)*(-alpha)*theta1*theta2/(theta2-1)*y_ss/k_ss*mu_ss^theta2/tau_ss...
    + zeta_ss/(theta2-1)*mu_ss*y_ss^(1-gamma)/tau_ss*(1-gamma);
h(11,3) = lambda_ss*c_ss^(-phic)*(-alpha)*(1-gamma)*y_ss^(1-gamma)/k_ss*(-1)...
    + lambda_ss*c_ss^(-phic)*alpha*(1-gamma)*(1+1/(theta2-1))*mu_ss*y_ss^(1-gamma)/k_ss*(-1)...
    + lambda_ss*c_ss^(-phic)*(-alpha)*theta1*theta2/(theta2-1)*y_ss/k_ss*mu_ss^theta2/tau_ss*(-1);
h(11,19) = lambda_ss*c_ss^(-phic)*alpha*(1-gamma)*(1+1/(theta2-1))*mu_ss*y_ss^(1-gamma)/k_ss...
    + lambda_ss*c_ss^(-phic)*(-alpha)*theta1*theta2/(theta2-1)*y_ss/k_ss*mu_ss^theta2/tau_ss*theta2...
    + zeta_ss/(theta2-1)*mu_ss*y_ss^(1-gamma)/tau_ss;
h(11,26) = zeta_ss/(theta2-1)*mu_ss*y_ss^(1-gamma)/tau_ss;

%Eqn 12: Planner's FOC wrt x
h(12,26) = zeta_ss;
h(12,40) = -beta*eta*zeta_ss;
h(12,27) = omega_ss*k_ss^alpha*(2*d2*x_ss+d1);
h(12,28) = omega_ss*k_ss^alpha*(2*d2*x_ss+d1);
h(12,3) = omega_ss*k_ss^alpha*(2*d2*x_ss+d1)*alpha;
h(12,18) = omega_ss*k_ss^alpha*2*d2*x_ss;

%Eqn 13: Planner's FOC wrt y
h(13,15) = c_ss^(-phic)*(-phic)...
    - ((1-gamma)*theta2-1)/(theta2-1)*c_ss^(-phic)*z_ss/y_ss*(-phic)...
    + lambda_ss*phic*c_ss^(-phic-1)*(-phic-1)...
    - lambda_ss*phic*c_ss^(-phic-1)*((1-gamma)*theta2-1)/(theta2-1)*z_ss/y_ss*(-phic-1)...
    + lambda_ss*(-phic)*c_ss^(-phic-1)*r_ss*(-phic-1)...
    + lambda_ss*(-phic)*c_ss^(-phic-1)*(-1)*((1-gamma)*theta2-1)/(theta2-1)*z_ss/y_ss*r_ss*(-phic-1)...
    + lambda_ss*(-phic)*c_ss^(-phic-1)*(1-delta)*(-phic-1)...
    + lambda_ss*(-phic)*c_ss^(-phic-1)*(-1)*((1-gamma)*theta2-1)/(theta2-1)*z_ss/y_ss*(1-delta)*(-phic-1)...
    + lambda_ss*c_ss^(-phic)*alpha/k_ss*(-phic)...
    + lambda_ss*c_ss^(-phic)*(-alpha)*(1-gamma)^2*tau_ss*y_ss^(-gamma)/k_ss*(-phic)...
    + lambda_ss*c_ss^(-phic)*alpha*(1-gamma)*(1-gamma-gamma/(theta2-1))*tau_ss*mu_ss*y_ss^(-gamma)/k_ss*(-phic)...
    + lambda_ss*c_ss^(-phic)*(-alpha)*theta1*(1-theta2*gamma/(theta2-1))/k_ss*mu_ss^theta2*(-phic);
h(13,21) = -((1-gamma)*theta2-1)/(theta2-1)*c_ss^(-phic)*z_ss/y_ss...
    - lambda_ss*phic*c_ss^(-phic-1)*((1-gamma)*theta2-1)/(theta2-1)*z_ss/y_ss...
    + lambda_ss*(-phic)*c_ss^(-phic-1)*(-1)*((1-gamma)*theta2-1)/(theta2-1)*z_ss/y_ss*r_ss...
    + lambda_ss*(-phic)*c_ss^(-phic-1)*(-1)*((1-gamma)*theta2-1)/(theta2-1)*z_ss/y_ss*(1-delta);
h(13,20) = -((1-gamma)*theta2-1)/(theta2-1)*c_ss^(-phic)*z_ss/y_ss*(-1)...
    - lambda_ss*phic*c_ss^(-phic-1)*((1-gamma)*theta2-1)/(theta2-1)*z_ss/y_ss*(-1)...
    + lambda_ss*(-phic)*c_ss^(-phic-1)*(-1)*((1-gamma)*theta2-1)/(theta2-1)*z_ss/y_ss*r_ss*(-1)...
    + lambda_ss*(-phic)*c_ss^(-phic-1)*(-1)*((1-gamma)*theta2-1)/(theta2-1)*z_ss/y_ss*(1-delta)*(-1)...
    + lambda_ss*c_ss^(-phic)*(-alpha)*(1-gamma)^2*tau_ss*y_ss^(-gamma)/k_ss*(-gamma)...
    + lambda_ss*c_ss^(-phic)*alpha*(1-gamma)*(1-gamma-gamma/(theta2-1))*tau_ss*mu_ss*y_ss^(-gamma)/k_ss*(-gamma)...
    - zeta_ss*(1-gamma)*y_ss^(-gamma)*(-gamma)...
    - zeta_ss*(gamma/(theta2-1)+gamma-1)*mu_ss*y_ss^(-gamma)*(-gamma);
h(13,25) = lambda_ss*phic*c_ss^(-phic-1)...
    - lambda_ss*phic*c_ss^(-phic-1)*((1-gamma)*theta2-1)/(theta2-1)*z_ss/y_ss;
h(13,11) = lambda_ss*(-phic)*c_ss^(-phic-1)*r_ss...
    + lambda_ss*(-phic)*c_ss^(-phic-1)*(-1)*((1-gamma)*theta2-1)/(theta2-1)*z_ss/y_ss*r_ss...
    + lambda_ss*(-phic)*c_ss^(-phic-1)*(1-delta)...
    + lambda_ss*(-phic)*c_ss^(-phic-1)*(-1)*((1-gamma)*theta2-1)/(theta2-1)*z_ss/y_ss*(1-delta)...
    + lambda_ss*c_ss^(-phic)*alpha/k_ss...
    + lambda_ss*c_ss^(-phic)*(-alpha)*(1-gamma)^2*tau_ss*y_ss^(-gamma)/k_ss...
    + lambda_ss*c_ss^(-phic)*alpha*(1-gamma)*(1-gamma-gamma/(theta2-1))*tau_ss*mu_ss*y_ss^(-gamma)/k_ss...
    + lambda_ss*c_ss^(-phic)*(-alpha)*theta1*(1-theta2*gamma/(theta2-1))/k_ss*mu_ss^theta2;
h(13,24) = lambda_ss*(-phic)*c_ss^(-phic-1)*r_ss...
    + lambda_ss*(-phic)*c_ss^(-phic-1)*(-1)*((1-gamma)*theta2-1)/(theta2-1)*z_ss/y_ss*r_ss;
h(13,3) = lambda_ss*c_ss^(-phic)*alpha/k_ss*(-1)...
    + lambda_ss*c_ss^(-phic)*(-alpha)*(1-gamma)^2*tau_ss*y_ss^(-gamma)/k_ss*(-1)...
    + lambda_ss*c_ss^(-phic)*alpha*(1-gamma)*(1-gamma-gamma/(theta2-1))*tau_ss*mu_ss*y_ss^(-gamma)/k_ss*(-1)...
    + lambda_ss*c_ss^(-phic)*(-alpha)*theta1*(1-theta2*gamma/(theta2-1))/k_ss*mu_ss^theta2*(-1);
h(13,23) = lambda_ss*c_ss^(-phic)*(-alpha)*(1-gamma)^2*tau_ss*y_ss^(-gamma)/k_ss...
    + lambda_ss*c_ss^(-phic)*alpha*(1-gamma)*(1-gamma-gamma/(theta2-1))*tau_ss*mu_ss*y_ss^(-gamma)/k_ss;
h(13,19) = lambda_ss*c_ss^(-phic)*alpha*(1-gamma)*(1-gamma-gamma/(theta2-1))*tau_ss*mu_ss*y_ss^(-gamma)/k_ss...
    + lambda_ss*c_ss^(-phic)*(-alpha)*theta1*(1-theta2*gamma/(theta2-1))/k_ss*mu_ss^theta2*theta2...
    - zeta_ss*(gamma/(theta2-1)+gamma-1)*mu_ss*y_ss^(-gamma);
h(13,26) = -zeta_ss*(1-gamma)*y_ss^(-gamma)...
    - zeta_ss*(gamma/(theta2-1)+gamma-1)*mu_ss*y_ss^(-gamma);
h(13,27) = omega_ss;

%Eqn 14: Planner's FOC wrt k
h(14,15) = -c_ss^(-phic)*(-phic)...
    + lambda_ss*(-phic)*c_ss^(-phic-1)*(-phic-1)...
    + lambda_ss*phic*c_ss^(-phic-1)*r_ss*(-phic-1)...
    + lambda_ss*phic*c_ss^(-phic-1)*(1-delta)*(-phic-1);
h(14,29) = beta*(1-delta)*c_ss^(-phic)*(-phic)...
    + lambda_ss*phic*c_ss^(-phic-1)*(1-delta)*beta*(-phic-1)...
    + lambda_ss*beta*(-phic)*c_ss^(-phic-1)*(1-delta)*r_ss*(-phic-1)...
    + lambda_ss*beta*(-phic)*c_ss^(-phic-1)*(1-delta)^2*(-phic-1)...
    + lambda_ss*beta*c_ss^(-phic)*(-alpha)*y_ss*k_ss^(-2)*(-phic)...
    + lambda_ss*beta*c_ss^(-phic)*alpha*(1-gamma)*tau_ss*y_ss^(1-gamma)*k_ss^(-2)*(-phic)...
    + lambda_ss*beta*c_ss^(-phic)*(-alpha)*(1-gamma)*tau_ss*mu_ss*y_ss^(1-gamma)*k_ss^(-2)*(-phic)...
    + lambda_ss*beta*c_ss^(-phic)*alpha*theta1*y_ss*k_ss^(-2)*mu_ss^theta2*(-phic);
h(14,39) = lambda_ss*phic*c_ss^(-phic-1)*(1-delta)*beta;
h(14,25) = lambda_ss*(-phic)*c_ss^(-phic-1)...
    + lambda_ss*beta*(-phic)*c_ss^(-phic-1)*(1-delta)*r_ss...
    + lambda_ss*beta*(-phic)*c_ss^(-phic-1)*(1-delta)^2 ...
    + lambda_ss*beta*c_ss^(-phic)*(-alpha)*y_ss*k_ss^(-2)...
    + lambda_ss*beta*c_ss^(-phic)*alpha*(1-gamma)*tau_ss*y_ss^(1-gamma)*k_ss^(-2)...
    + lambda_ss*beta*c_ss^(-phic)*(-alpha)*(1-gamma)*tau_ss*mu_ss*y_ss^(1-gamma)*k_ss^(-2)...
    + lambda_ss*beta*c_ss^(-phic)*alpha*theta1*y_ss*k_ss^(-2)*mu_ss^theta2;
h(14,38) = lambda_ss*beta*(-phic)*c_ss^(-phic-1)*(1-delta)*r_ss;
h(14,20) = lambda_ss*beta*c_ss^(-phic)*(-alpha)*y_ss*k_ss^(-2)...
    + lambda_ss*beta*c_ss^(-phic)*alpha*(1-gamma)*tau_ss*y_ss^(1-gamma)*k_ss^(-2)*(1-gamma)...
    + lambda_ss*beta*c_ss^(-phic)*(-alpha)*(1-gamma)*tau_ss*mu_ss*y_ss^(1-gamma)*k_ss^(-2)*(1-gamma)...
    + lambda_ss*beta*c_ss^(-phic)*alpha*theta1*y_ss*k_ss^(-2)*mu_ss^theta2;
h(14,3) = lambda_ss*beta*c_ss^(-phic)*(-alpha)*y_ss*k_ss^(-2)*(-2)...
    + lambda_ss*beta*c_ss^(-phic)*alpha*(1-gamma)*tau_ss*y_ss^(1-gamma)*k_ss^(-2)*(-2)...
    + lambda_ss*beta*c_ss^(-phic)*(-alpha)*(1-gamma)*tau_ss*mu_ss*y_ss^(1-gamma)*k_ss^(-2)*(-2)...
    + lambda_ss*beta*c_ss^(-phic)*alpha*theta1*y_ss*k_ss^(-2)*mu_ss^theta2*(-2);
h(14,23) = lambda_ss*beta*c_ss^(-phic)*alpha*(1-gamma)*tau_ss*y_ss^(1-gamma)*k_ss^(-2)...
    + lambda_ss*beta*c_ss^(-phic)*(-alpha)*(1-gamma)*tau_ss*mu_ss*y_ss^(1-gamma)*k_ss^(-2);
h(14,19) = lambda_ss*beta*c_ss^(-phic)*(-alpha)*(1-gamma)*tau_ss*mu_ss*y_ss^(1-gamma)*k_ss^(-2)...
    + lambda_ss*beta*c_ss^(-phic)*alpha*theta1*y_ss*k_ss^(-2)*mu_ss^theta2*theta2;
h(14,11) = lambda_ss*phic*c_ss^(-phic-1)*r_ss...
    + lambda_ss*phic*c_ss^(-phic-1)*(1-delta);
h(14,24) = lambda_ss*phic*c_ss^(-phic-1)*r_ss;
h(14,41) = omega_ss*(-beta)*(1-d0)*alpha*k_ss^(alpha-1)...
    + omega_ss*(-beta)*(-d2*x_ss^2)*alpha*k_ss^(alpha-1)...
    + omega_ss*(-beta)*(-d1*x_ss)*alpha*k_ss^(alpha-1);
h(14,42) = omega_ss*(-beta)*(1-d0)*alpha*k_ss^(alpha-1)...
    + omega_ss*(-beta)*(-d2*x_ss^2)*alpha*k_ss^(alpha-1)...
    + omega_ss*(-beta)*(-d1*x_ss)*alpha*k_ss^(alpha-1);
h(14,17) = omega_ss*(-beta)*(1-d0)*alpha*k_ss^(alpha-1)*(alpha-1)...
    + omega_ss*(-beta)*(-d2*x_ss^2)*alpha*k_ss^(alpha-1)*(alpha-1)...
    + omega_ss*(-beta)*(-d1*x_ss)*alpha*k_ss^(alpha-1)*(alpha-1);
h(14,32) = omega_ss*(-beta)*(-d2*x_ss^2)*alpha*k_ss^(alpha-1)*2 ...
    + omega_ss*(-beta)*(-d1*x_ss)*alpha*k_ss^(alpha-1);



%These are parameters for the AIM solution method
condn = .0000001;
uprbnd = 1;

[b,rts,ia,nexact,nnumeric,lgroots,mcode,q] = aim_eig(h,neq,nlag,nlead,condn,uprbnd);

%Get the observable structure

scof = obstruct(h,b,neq,nlag,nlead);

%Use this to get the reduced form of the structural model

S_minus1 = scof(:,1:neq); %This is for nlag = 1
S_zero = scof(:,neq+1:2*neq);

B_minus1 = -inv(S_zero)*S_minus1;
B_zero = inv(S_zero);

%%%%%%%%%%%%%%%%%%%%%%%%%
%Generate impulse response functions

T = 100; %periods for response
state = zeros(neq,T);
eps = zeros(neq,T);

eps(:,1) = [imprespsize;0;0;0;0;0;0;0;0;0;0;0;0;0]; %shock to a 

state(:,1) = B_minus1*[0;0;0;0;0;0;0;0;0;0;0;0;0;0] + B_zero*eps(:,1);
for t = 2:T
    state(:,t) = B_minus1*state(:,t-1) + B_zero*eps(:,t);
end;

%Examine how this affects emissions in each period, capital, abatement, and
%production

c = state(1,:);
e = state(2,:);
k = state(3,:);
x = state(4,:);
mu = state(5,:);
y = state(6,:);
z = state(7,:);
i = state(8,:);
tau = state(9,:);
r = state(10,:);
lambda = state(11,:);
zeta = state(12,:);
omega = state(13,:);
a = state(14,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Now re-solve with the planner's information asymmetry
%Now variables are [c,e,k,x,mu,y,z,i,q,r,lambda,zeta,omega,a]
%Eqn 11: Planner's policy option (q as function of yt-1)
h(11,:) = zeros(1,neq*(nlag+1+nlead));
h(11,23) = 1;
h(11,6) = -d;

%Eqns 12-14: Set the Lagrangian multiplier to be constant
h(12:14,:) = zeros(3,neq*(nlag+1+nlead));
h(12,25) = 1;
h(13,26) = 1;
h(14,27) = 1;

%Eqn 9: Firm's foc wrt k
h(9,:) = zeros(1,neq*(nlag+1+nlead));
h(9,20) = alpha*y_ss/k_ss...
    + theta1*(-1+theta2*(1-gamma))*alpha*y_ss/k_ss*mu_ss^theta2...
    - theta1*theta2*(1-gamma)*alpha*y_ss/k_ss*mu_ss^(theta2-1);
h(9,3) = alpha*y_ss/k_ss*(-1)...
    + theta1*(-1+theta2*(1-gamma))*alpha*y_ss/k_ss*mu_ss^theta2*(-1)...
    - theta1*theta2*(1-gamma)*alpha*y_ss/k_ss*mu_ss^(theta2-1)*(-1);
h(9,19) = theta1*(-1+theta2*(1-gamma))*alpha*y_ss/k_ss*mu_ss^theta2*theta2...
    - theta1*theta2*(1-gamma)*alpha*y_ss/k_ss*mu_ss^(theta2-1)*(theta2-1);
h(9,24) = -r_ss;

%Eqn 10: Firm's constraint
h(10,:) = zeros(1,neq*(nlag+1+nlead));
h(10,16) = 1;
h(10,23) = -1;

[b,rts,ia,nexact,nnumeric,lgroots,mcode,q] = aim_eig(h,neq,nlag,nlead,condn,uprbnd);

%Get the observable structure

scof = obstruct(h,b,neq,nlag,nlead);

%Use this to get the reduced form of the structural model

S_minus1 = scof(:,1:neq); %This is for nlag = 1
S_zero = scof(:,neq+1:2*neq);

B_minus1 = -inv(S_zero)*S_minus1;
B_zero = inv(S_zero);

%%%%%%%%%%%%%%%%%%%%%%%%%
%Generate impulse response functions

T = 100; %periods for response
state_2 = zeros(neq,T);
eps = zeros(neq,T);

eps(:,1) = [imprespsize;0;0;0;0;0;0;0;0;0;0;0;0;0]; %shock to a 

state_2(:,1) = B_minus1*[0;0;0;0;0;0;0;0;0;0;0;0;0;0] + B_zero*eps(:,1);
for t = 2:T
    state_2(:,t) = B_minus1*state(:,t-1) + B_zero*eps(:,t);
end;

%Examine how this affects emissions in each period, capital, abatement, and
%production

c_2 = state_2(1,:);
e_2 = state_2(2,:);
k_2 = state_2(3,:);
x_2 = state_2(4,:);
mu_2 = state_2(5,:);
y_2 = state_2(6,:);
z_2 = state_2(7,:);
i_2 = state_2(8,:);
tau_2 = state_2(9,:);
r_2 = state_2(10,:);
lambda_2 = state_2(11,:);
zeta_2 = state_2(12,:);
omega_2 = state_2(13,:);
a_2 = state_2(14,:);

diff = (c-c_2)*(c-c_2)'; %If we just want to match consumption

%%%%
%What I want is to solve the function IRFdiff for the optimal d

%[d_opt,val] = fminsearch(@IRFdiff_quota,d);

%%%%%%%%%%%%%
%Simulate business cycle and do welfare analysis
state_2 = zeros(neq,T);
eps = zeros(neq,T);

load commonshocks; %This contains the matrix commonshocks, which has a series
                    %of shocks, that are common for comparison of
                    %different policies

for t = 1:T
    eps(1,t) = sd*commonshocks(t); %shock to a 
end

state_2(:,1) = B_minus1*[0;0;0;0;0;0;0;0;0;0;0;0;0;0] + B_zero*eps(:,1);
for t = 2:T
    state_2(:,t) = B_minus1*state_2(:,t-1) + B_zero*eps(:,t);
end;

%Examine how this affects emissions in each period, capital, abatement, and
%production

c_2 = state_2(1,:);
e_2 = state_2(2,:);
k_2 = state_2(3,:);
x_2 = state_2(4,:);
mu_2 = state_2(5,:);
y_2 = state_2(6,:);
z_2 = state_2(7,:);
i_2 = state_2(8,:);
tau_2 = state_2(9,:);
r_2 = state_2(10,:);
lambda_2 = state_2(11,:);
zeta_2 = state_2(12,:);
omega_2 = state_2(13,:);
a_2 = state_2(14,:);

Cyclematrix = [y_2',e_2',k_2',x_2'];
plot(Cyclematrix);

%Evaluate welfare
U_ss = c_ss^(1-phic)/(1-phic);
for t = 1:T
    c_level(t) = c_ss*(1+c_2(t));
    U(t) = beta^t*c_level(t)^(1-phic)/(1-phic);
end;

welfare = sum(U);