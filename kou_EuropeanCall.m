%s0=100; r=0.05; k= 98; sigma=0.16; T=0.5; eta1=10; eta2=5; p=0.4; lambda=1;
function call_kou = kou_EuropeanCall(s0, k, sigma, r, T, eta1, eta2, p, lambda)
limit = 15; % kou suggest that series upper limit 10 to 15
zeta = (1-p)*eta2/(eta2+1) + p*eta1/(eta1-1) - 1;
p_tilde = p/(1+zeta)*eta1/(eta1-1);
lambda_tilde = lambda*(zeta + 1);
eta1_tilde = eta1 - 1;
eta2_tilde = eta2 + 1;
mu1 = r + 0.5*sigma*sigma - lambda*zeta;
mu2 = r - 0.5*sigma*sigma - lambda*zeta;
a = log(k/s0);

call_kou = s0*phi(mu1, sigma, lambda_tilde, p_tilde, eta1_tilde, eta2_tilde, a, T, limit)...
    - k*exp(-r*T)*phi(mu2, sigma, lambda, p, eta1, eta2, a, T, limit);
end

function result = phi(mu, sigma, lambda, p, eta1, eta2, a, T, limit)
coef1 = exp((sigma*eta1)^2*T/2)/(sigma*sqrt(2*pi*T));
coef2 = exp((sigma*eta2)^2*T/2)/(sigma*sqrt(2*pi*T));
P_outter = 0;
Q_outter = 0;
for n=1 : limit
    pi_n = exp(-lambda*T)*(lambda*T)^n/factorial(n);
    P_inner = 0;
    Q_inner = 0;
    for k=1 : n
        P_value = 0;
        Q_value = 0;
        if k < n
            for m = k : n-1
                x = n - k - 1;
                y = m - k;
                P_prob = p^m*(1-p)^(n - m);
                Q_prob = p^(n-m)*(1-p)^m;
                P_value = P_value + P_prob*factorial(x)/(factorial(y)*factorial(x - y))...
                    *factorial(n)/(factorial(m)*factorial(n-m))*(eta1/(eta1 + eta2))^(m - k)*(eta2/(eta1 + eta2))^(n - m);
                Q_value = Q_value + Q_prob*factorial(x)/(factorial(y)*factorial(x - y))...
                    *factorial(n)/(factorial(m)*factorial(n-m))*(eta1/(eta1 + eta2))^(n - m)*(eta2/(eta1 + eta2))^(m - k);
            end
        else
            P_value = P_value + p^n;
            Q_value = Q_value + (1-p)^n;
        end
        P_inner = P_inner + P_value*(sigma*eta1*sqrt(T))^k...
            *I_function(a-mu*T, -eta1, -1/(sigma*sqrt(T)), -sigma*eta1*sqrt(T), k-1, limit);
        Q_inner = Q_inner + Q_value*(sigma*eta2*sqrt(T))^k...
            *I_function(a-mu*T, eta2, 1/(sigma*sqrt(T)), -sigma*eta2*sqrt(T), k-1, limit);     
    end
    P_outter = P_outter + pi_n*P_inner;
    Q_outter = Q_outter + pi_n*Q_inner;
end
result = coef1*P_outter + coef2*Q_outter ...
    + exp(-lambda*T)*(lambda*T)^0/factorial(0)*normcdf(-(a-mu*T)/(sigma*sqrt(T)));
end

function result = I_function(c, alpha, beta, delta, n, limit)
coef = -exp(alpha*c)/alpha;
if (beta>0) && (alpha ~=0)
    const1 = (beta/alpha)^(n+1)*sqrt(2*pi)/beta*exp(alpha*delta/beta+alpha^2/(2*beta^2))...
        *normcdf(-beta*c+delta+alpha/beta);
    result = coef*sum_Hh(beta*c-delta, n, alpha, beta, limit) + const1;
elseif (beta<0) && (alpha<0)
    const2 = -(beta/alpha)^(n+1)*sqrt(2*pi)/beta*exp(alpha*delta/beta+alpha^2/(2*beta^2))...
        *normcdf(beta*c-delta-alpha/beta);
    result = coef*sum_Hh(beta*c-delta, n, alpha, beta, limit) + const2;
end
end

function result = sum_Hh(x, n, alpha, beta, limit)
count = 0;
Hh = zeros(1, limit+2);
Hh(1,1) = exp(-x^2/2);
Hh(1,2) = sqrt(2*pi)*normcdf(-x);
for k=0 : n
    coef = (beta/alpha)^(n - k);
    for m=1 : k
        h = m+2;
        Hh(1,h) = 1/m*(Hh(1,h-2) - x*Hh(1,h-1));
    end
    count = count + coef*Hh(1,k+2);
end
result = count;
end