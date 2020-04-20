%k=90;t=0; s0 =100; v0=0.09; theta=0.348; sigma=0.39; r=0.05; kappa=1.15; rho=-0.64; T=1.5;
function call = GeoAsianOption(s0,k,v0,theta,sigma,r,kappa,rho,T,t)
a1 = (2*v0)/sigma^2;
a2 = (2*kappa*theta)/sigma^2;
a3 =(T-t)/T*log(s0) + ((r*sigma - kappa*theta*rho)*(T-t)^2)/(2*sigma*T) - rho*(T-t)/(sigma*T)*v0;
a4 = log(s0) - rho/sigma*v0 + (r - rho*kappa*theta/sigma)*(T-t);
a5 = (kappa*v0 + kappa^2*theta*(T-t))/sigma^2;
call = exp(-r*(T-t))*((Characteristic(0,a1,a2,a3,a4,a5,sigma,kappa,rho,T,t,1) -...
    k)/2 + GeoAsianIntegral(a1,a2,a3,a4,a5,k,sigma,kappa,rho,T,t));
end

function result = GeoAsianIntegral(a1,a2,a3,a4,a5,k,sigma,kappa,rho,T,t)
result = 1/pi*quadl(@RealSumCharacteristic,0,100000,[],[],a1,a2,a3,a4,a5,k,sigma,kappa,rho,T,t);
end

function result = RealSumCharacteristic(xi,a1,a2,a3,a4,a5,k,sigma,kappa,rho,T,t)
result = real((Characteristic(xi,a1,a2,a3,a4,a5,sigma,kappa,rho,T,t,2) -...
    k*Characteristic(xi,a1,a2,a3,a4,a5,sigma,kappa,rho,T,t,3)).*(exp(-1i*xi*log(k)))./(1i*xi));
end

function phi = Characteristic(xi,a1,a2,a3,a4,a5,sigma,kappa,rho,T,t,type)
[H,H_Tilde,s,w] = H_calculation(xi,sigma,kappa,rho,T,t,type);
phi = exp(-a1*H_Tilde./H - a2*log(H) + a3*s + a4*w + a5);
end

function [H_Sum,H_Tilde_Sum,s,w]=H_calculation(xi,sigma,kappa,rho,T,t,type)
if type == 1
    s = 1;
    w = 0;
elseif type == 2
    s = 1 +1i*xi;
    w = 0;
else
    s =  1i*xi;
    w = 0;
end

E = -s.^2*sigma^2*(1 - rho^2)*(T-t)^2;
F =  (s*sigma*T*(sigma - 2*rho*kappa ) - 2*s*w*sigma^2*T*(1 - rho^2))*(T-t);
G = T*(kappa^2*T - 2*s*rho*sigma - w*(2*rho*kappa - sigma)*sigma*T - w^2*(1 - rho^2)*sigma^2*T);

m = size(E,2);
n = 10;
h = zeros(m,n+3);
h_Tilde = zeros(m,n);
h(:,1) = 0;
h(:,2) = 0;
h(:,3) = 1;
h(:,4) = (T-t)*(kappa - w*rho*sigma)/2;
h_Tilde(:,1) = 1/(T-t)*h(:,4);

for j = 2:n
    k = j + 3;
    h(:,k) = (T-t)^2/(4*j*(j-1)*T^2)*(E'.*h(:,k-4) + F'.*h(:,k-3) + G'.*h(:,k-2));
    h_Tilde(:,j) = j/(T-t)*h(:,k);
end

H_Sum = (sum(h(:,3:end),2))';
H_Tilde_Sum = (sum(h_Tilde(:,1:end),2))';
end