tic
s0 = 100; r=0.05; theta=0.348; sigma=0.39; kappa=1.15; rho=-0.64; v0=0.09; t=0;
M = 50000;       % number of paths 
T = 1;
k = 90;
h = 0.001;      % steps
N = T/h;
b = 1;

S_Ari1 = s0*ones(M,1);
S_Ari2 = s0*ones(M,1);
V1 = v0*ones(M,1);
V2 = v0*ones(M,1);

cumsum_Ari1 = s0*ones(M,1);
cumsum_Geo1 = log(s0)*ones(M,1);
for i = 1:N
    %Euler scheme
    dW1 = randn(M,1)*sqrt(h);
    dW2 = rho*dW1 + sqrt(1-rho^2)*randn(M,1)*sqrt(h);
    S_Ari1 = S_Ari1.*(1 + r*h + sqrt(abs(V1)).*dW1);
    V1= V1 + kappa*(theta-V1)*h + sigma*sqrt(abs(V1)).*dW2;
    S_Geo1 = log(S_Ari1);
    cumsum_Ari1 = cumsum_Ari1 + S_Ari1;
    cumsum_Geo1 = cumsum_Geo1 + S_Geo1;
end
payoff_Ari1 = exp(-r*T)*max(0,cumsum_Ari1/(N+1)-k);
stdpayoff_Ari1 = std(payoff_Ari1);
call1 = mean(payoff_Ari1);

% control variates
payoff_Geo1 = exp(-r*T)*max(0,exp(cumsum_Geo1/(N+1))-k);
payoff_Ari1_cv=payoff_Ari1-b*(payoff_Geo1- GeoAsianOption(s0,k,v0,theta,sigma,r,kappa,rho,T,t));
stdpayoff_Ari1_cv = std(payoff_Ari1_cv);
call1_cv = mean(payoff_Ari1_cv)
call1_cv_left = call1_cv - 1.96*stdpayoff_Ari1_cv/sqrt(M)
call1_cv_right = call1_cv + 1.96*stdpayoff_Ari1_cv/sqrt(M)
toc 
tic
cumsum_Ari2 = s0*ones(M,1);
cumsum_Geo2 = log(s0)*ones(M,1);
for i = 1:N
    %Milstein scheme
    dW1 = randn(M,1)*sqrt(h);
    dW2 = rho*dW1 + sqrt(1-rho^2)*randn(M,1)*sqrt(h);
    S_Ari2 = S_Ari2.*(1 + r*h + sqrt(abs(V2)).*dW1)+1/2*V2.*S_Ari2.*(dW1.^2-h);
    V2 = V2 + kappa*(theta-V2)*h + sigma*sqrt(abs(V2)).*dW2+1/4*sigma^2.*(dW2.^2-h);
    S_Geo2 = log(S_Ari2);
    cumsum_Ari2 = cumsum_Ari2 + S_Ari2;
    cumsum_Geo2 = cumsum_Geo2 + S_Geo2;
end
payoff_Ari2 = exp(-r*T)*max(0,cumsum_Ari2/(N+1)-k);
stdpayoff_Ari2 = std(payoff_Ari2);
call2 = mean(payoff_Ari2);

% control variates
payoff_Geo2 = exp(-r*T)*max(0,exp(cumsum_Geo2/(N+1))-k);
payoff_Ari2_cv=payoff_Ari2-b*(payoff_Geo2- GeoAsianOption(s0,k,v0,theta,sigma,r,kappa,rho,T,t));
stdpayoff_Ari2_cv = std(payoff_Ari2_cv);
call2_cv = mean(payoff_Ari2_cv)
call2_cv_left = call2_cv - 1.96*stdpayoff_Ari2_cv/sqrt(M)
call2_cv_right = call2_cv + 1.96*stdpayoff_Ari2_cv/sqrt(M)
toc