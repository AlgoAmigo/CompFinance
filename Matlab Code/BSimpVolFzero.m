function [impliedVol] = BSimpVolFzero(S0,K,T,r, marketValue,q) 
% Black-Scholes implied volatility. 

%S0: spot price
%K: strike price 
%r: riskFree rate
%T: time to maturity
%marketValue: observed call option market value
%q: stock dividend yield

%using fzero. 

%Adapted by Lorenzo Torricelli
%Modified by Aaron Wittmann and Jan Keesen

impliedVol = ones(length(T),length(K));
for i = 1:length(T)
    it = BSimpVolAnalyticBrSu(S0, T(i),marketValue);
    for j = 1:length(K)
        fun = @(x) 	BSprice(S0, K(j), r, T(i), x, q)-marketValue;
        impliedVol(i,j) = fzero(fun,it);
    end
end

return