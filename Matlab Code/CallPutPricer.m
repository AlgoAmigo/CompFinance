function [CallPrice, PutPrice, ForwardPrice]=CallPutPricer(S0, price, K, T, r, q)

%A pricer for call and purt option, given paths as inputs.

%price: array of terminal prices
%K: strike price
%T: expiration
%r: riskfree rate
%q: dividend yield
%sigma: volatility

%Author: Lorenzo Torricelli


df=exp(-r.*T);
CallPrice=df.*mean(max(price(:,end)-K, 0));
ForwardPrice=S0.*exp(-q).*T-df.*K; %no arbitrage model-free value. Needs not to be simulated.
PutPrice=CallPrice-ForwardPrice;
return