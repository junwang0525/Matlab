s0=100; r=0.05; k= 98; sigma=0.16; T=0.5; eta1=10; eta2=5; p=0.4; lambda=1;
nPath=10000;
nStep=1000;

dt = T/nStep;
s = s0*ones(nPath,nStep+1);
dW = randn(nPath,nStep)*sqrt(dt);
Nt = poissrnd(lambda*dt,[nPath,nStep]);

Bt = binornd(1,p,[nPath,nStep]);
Yt = exprnd((1/eta1),[nPath,nStep]).*Bt + exprnd((1/eta2),[nPath,nStep]).*(Bt-1);
J = Nt.*(exp(Yt)-1);
zeta = (1-p)*eta2/(eta2+1) + p*eta1/(eta1-1) - 1; %To make the process a martingale

for i = 1:nStep
    s(:,i+1) = s(:,i).*(1+(r - lambda*zeta)*dt + sigma*dW(:,i) + J(:,i));
end

t_space=linspace(0,T,nStep+1);
payoff_kou = exp(-r*T)*mean(max(s(:,end)-k,0));

%figure(1)
%plot(t_space,s);
%figure(2)
%histogram(Yt);
%hold on
%line([0,0], [0,5000], 'Color', 'red', 'LineWidth', 2)
%hold off
disp(payoff_kou);
