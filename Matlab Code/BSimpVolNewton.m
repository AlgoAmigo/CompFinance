function [impliedVol,n] = BSimpVolNewton(S0,K,T,r, marketValue,q) 
% Black-Scholes implied volatility. 

%S0: spot price
%K: strike price 
%r: riskFree rate
%T: time to maturity
%marketValue: observed call option market value
%q: stock dividend yield

%using Newton-Raphson method. 

%Adapted by Lorenzo Torricelli

%%Tolerance and max iterations



    maxiter = 50;

	epsilon = 1e-6;

% Use the Newton-Raphson Method to find Implied Volatility
n = 0;
sigma = BSimpVolAnalyticBrSu(S0, T, marketValue);
%we can use Brenner Subrahmanyam formula! to speed up
previousSigma = 0;   


while abs(sigma - previousSigma ) >epsilon,
  
    previousSigma = sigma;
	
	BSEvaluation = BSprice(S0, K, r, T, previousSigma, q);

    vegaEvaluation = BSvega(S0, K, r, T, previousSigma, q);
	
	sigma = previousSigma  - ((BSEvaluation- marketValue) ./ vegaEvaluation ); %newton iterator
	
	n = n + 1;
	if n > maxiter,
		error('No solution is found');
	end
end
if n == 0,
		error('No iteration performed');
end

impliedVol=sigma;


return