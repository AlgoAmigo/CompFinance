function [price,subordinator]=EX5_VG(simN, stepsN, S0, T, r, q, params)

%%Variance Gamma Levy subordinated model


%simN: number of simulations
% stepsN: number of time steps
%S0: starting price
%T: time horizon of simulation
%r: riskfree rate
%q: dividend yield

%Author: Lorenzo Torricelli

rng(123, 'twister'); %seeding the Mersenne Twister

kappa=params(1);  %%variance of the variance gamma subordinator:  variance subordinator
sigma=params(2);  %% drift of the subordinated Brownian motion: smile 
theta=params(3);  %%volatility of the subordinated brownian motion: skew 



compensator=-1/kappa*log(1+sigma^2*kappa/2+theta*kappa);   %%this is the Levy compensator to be included in the drift in order for the price model
%%to be a martisngale, i.e the characteristic exponent computed at -i

price=S0*ones(simN, stepsN+1); 

dt=T/(stepsN);


dWt=normrnd(0,1,simN, stepsN ); %Brownian motion subordinating

%IG parameters
lambda=dt/kappa;
beta=kappa;

 %%Generates gamma process subordinator
grnd=gamrnd( lambda,  beta, simN , stepsN );  
subordinator=[ zeros(simN,1) ,cumsum(grnd,2)];  

%compensator=0;

increment=(r-q-compensator)*dt+theta.*grnd+sigma*dWt.*sqrt(grnd);

price(:, 2:end)=S0*exp(cumsum(increment, 2));
%mart=exp(-r*T)*mean(price(:,end));  %%check of martingale property (q=0)

return