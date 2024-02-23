%Program to solve pollution/RBC model with AIM
%New specification of damages to output, not utility
%12/1/2010

clear;

%First define the parameter values
global beta delta rho alpha eta theta1 theta2 gamma phic d2 d1 d0 erowcoef;


%{
beta = 0.98267; %discount rate
delta = 0.025; %capital depreciation
rho = .95; %persistence of TFP shock
sd = 0.007; %standard deviation of shock
alpha = 0.36; %curvature of production function

eta = 0.9979; %pollution depreciation 
theta1 = .05607; %abatement cost equation parameters, from Nordhaus
theta2 = 2.8;
gamma = 1-.696; %1 - elasticity of emissions with respect to output
phic = 2; %CRRA for consumption
d2 = 1.4647*10^(-8); %damage function parameters, from Nordhaus
d1 = -6.6722*10^(-6);
d0 = 1.395*10^(-3);
damage_scale = 5.3024; %To scale the pollution levels and get the damage function correct
d2 = d2/damage_scale^2;
d1 = d1/damage_scale;
erowcoef = 4; %coefficient relating how rest-of-world emissions compare to domestic emissions
damage_scale = 5.3024*(1-.9979)/(1-eta); %new values of eta mean this needs to be rescaled
d2 = 1.4647*10^(-8)/damage_scale^2;
d1 = -6.6722*10^(-6)/damage_scale;
%}


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


imprespsize = .01;

%Solve for steady state solution based on first order condition
a_ss = 1;

%guess at the steady state values
%Guess values are from the old model's solution

k_g = 27;
e_g = 2;
mu_g =0.001;

guess = [k_g,e_g,mu_g]; %guess for vector of steady state values
options=optimset('Display','iter','TolFun',1e-15,'MaxFunEvals',5000);
ss_sol = fsolve(@steadystate,guess,options);
clear options guess;

k_ss = ss_sol(1);
e_ss = ss_sol(2);
mu_ss = ss_sol(3);
i_ss = delta*k_ss;
x_ss = erowcoef*e_ss/(1-eta);
y_ss = (1-d2*(x_ss)^2-d1*(x_ss)-d0)*k_ss^alpha;
z_ss = theta1*mu_ss^theta2*y_ss;
c_ss = y_ss - i_ss - z_ss;

%Write out the structural coefficient matrix h
%The order of the parameters is [c,e,k,x,mu,y,z,i,a]
neq = 9;
nlag = 1;
nlead = 1;
h = zeros(neq,neq*(nlag+1+nlead));
%Columns 1-9 are variables at t-1
%Columns 10-18 are variables at t
%Columns 19-27 are variables at t+1

%Eqn 1: evolution of shock a
h(1,9) = -rho;
h(1,18) = 1;

%Eqn 2: Budget constraint
h(2,10) = c_ss;
h(2,17) = i_ss;
h(2,16) = z_ss;
h(2,15) = -y_ss;

%Eqn 3: Capital accumulation
h(3,12) = k_ss;
h(3,3) = -k_ss*(1-delta);
h(3,17) = -i_ss;

%Eqn 4: pollution stock
h(4,4) = -eta;
h(4,11) = -(1-eta);
h(4,13) = 1;

%Eqn 5: Emissions
h(5,11) = e_ss;
h(5,15) = y_ss^(1-gamma)*(1-gamma)*(mu_ss-1);
h(5,14) = mu_ss*y_ss^(1-gamma);

%Eqn 6: Abatement
h(6,16) = 1;
h(6,14) = -theta2;
h(6,15) = -1;

%Eqn 7: Production function
h(7,15) = y_ss;
h(7,18) = -(1-d2*x_ss^2-d1*x_ss-d0)*k_ss^alpha;
h(7,3) = -(1-d2*x_ss^2-d1*x_ss-d0)*k_ss^alpha*alpha;
h(7,13) = k_ss^alpha*(2*d2*x_ss^2+d1*x_ss);

%Eqn 8: FOC for choice of kt
C81 = alpha*y_ss*k_ss^(-1)*(1-theta1*mu_ss^theta2...
    -theta1*theta2*(1-gamma)*y_ss^(gamma-1)*mu_ss^(theta2-1)*e_ss);
C82 = alpha*y_ss*k_ss^(-1)*...
    (-theta1*theta2*(1-gamma)*y_ss^(gamma-1)*mu_ss^(theta2-1)*e_ss);
h(8,10) = phic;
h(8,19) = -beta*phic*(C81+1-delta);
h(8,24) = beta*C81 + beta*C82*(gamma-1);
h(8,12) = -beta*C81;
h(8,23) = beta*(alpha*y_ss*k_ss^(-1)*(-theta1*mu_ss^theta2*theta2))...
    +beta*C82*(theta2-1);
h(8,20) = beta*C82;

%Eqn 9: FOC for choice of xt
h(9,10) = -phic*(-(1-theta1*mu_ss^theta2)*(2*d2*x_ss+d1)*k_ss^alpha...
    +theta1*theta2*mu_ss^(theta2-1)*y_ss^gamma...
    +theta1*theta2*(1-gamma)*mu_ss^(theta2-1)*y_ss^(gamma-1)*e_ss...
    *(2*d2*x_ss+d1)*k_ss^alpha);
h(9,14) = (theta1*mu_ss^theta2)*(2*d2*x_ss+d1)*k_ss^alpha*theta2...
    +theta1*theta2*mu_ss^(theta2-1)*y_ss^gamma*(theta2-1)...
    +theta1*theta2*(1-gamma)*mu_ss^(theta2-1)*y_ss^(gamma-1)*e_ss*(2*d2*x_ss+d1)*k_ss^alpha*(theta2-1);
h(9,13) = (-(1-theta1*mu_ss^theta2)*2*d2*x_ss*k_ss^alpha)...
    +theta1*theta2*(1-gamma)*mu_ss^(theta2-1)*y_ss^(gamma-1)*e_ss*2*d2*x_ss*k_ss^alpha;
h(9,18) = (-(1-theta1*mu_ss^theta2)*(2*d2*x_ss+d1)*k_ss^alpha)...
    +theta1*theta2*(1-gamma)*mu_ss^(theta2-1)*y_ss^(gamma-1)*e_ss*(2*d2*x_ss+d1)*k_ss^alpha;
h(9,3) = (-(1-theta1*mu_ss^theta2)*(2*d2*x_ss+d1)*k_ss^alpha)*alpha...
    +theta1*theta2*(1-gamma)*mu_ss^(theta2-1)*y_ss^(gamma-1)*e_ss*(2*d2*x_ss+d1)*k_ss^alpha*alpha;
h(9,15) = theta1*theta2*mu_ss^(theta2-1)*y_ss^gamma*gamma...
    +theta1*theta2*(1-gamma)*mu_ss^(theta2-1)*y_ss^(gamma-1)*e_ss*(2*d2*x_ss+d1)*k_ss^alpha*(gamma-1);
h(9,11) = theta1*theta2*(1-gamma)*mu_ss^(theta2-1)*y_ss^(gamma-1)*e_ss*(2*d2*x_ss+d1)*k_ss^alpha;
h(9,19) = beta*theta1*theta2*eta*mu_ss^(theta2-1)*y_ss^(gamma-1)*phic;
h(9,23) = -beta*theta1*theta2*eta*mu_ss^(theta2-1)*y_ss^(gamma-1)*(theta2-1);
h(9,24) = -beta*theta1*theta2*eta*mu_ss^(theta2-1)*y_ss^(gamma-1)*(gamma-1);

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

eps(:,1) = [imprespsize;0;0;0;0;0;0;0;0]; %shock to a 

state(:,1) = B_minus1*[0;0;0;0;0;0;0;0;0] + B_zero*eps(:,1);
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
a = state(9,:);

%Create matrices to input in the functions to create figures
IRFAmatrix = [a',y',k',c'];
plot(IRFAmatrix);
%createfigure_IRF_A(IRFAmatrix);

IRFBmatrix = [a',z',e',x',mu'];
plot(IRFBmatrix);
%createfigure_IRF_B(IRFBmatrix);


%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generate business cycle simulations

T = 100; %periods for response
state = zeros(neq,T);
eps = zeros(neq,T);

load commonshocks; %This contains the matrix commonshocks, which has a series
                    %of shocks, that are common for comparison of
                    %different policies

for t = 1:T
    eps(1,t) = sd*commonshocks(t); 
end;

state(:,1) = B_minus1*[0;0;0;0;0;0;0;0;0] + B_zero*eps(:,1);
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
a = state(9,:);

%Create matrices to input in the functions to create figures
Cyclematrix = [y',e',k',x'];
plot(Cyclematrix);
%createfigure_Cycle(Cyclematrix);

%Find moments of simulation
std(y); %standard deviation of detrended output
std(e);
corrcoef(y,e); %correlations between output and emissions
