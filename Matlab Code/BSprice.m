
function price= BSprice(S0, K, r, T, sigma, q)
%%Black Scholes price implied volatility of a dividend-payinng asset 

%Author: Lorenzo Torricelli

%S0: spot price: scalar
%K: strike price:  Matrix
%r: riskFree rate:  scalar
%T: time to maturity: Matrix
%sigma: black scholes volatility:  Matrix
%q: stock dividend yield:   scalar

realVol=sigma.*sqrt(T); %realised volatility 
 forward=S0.*exp((r-q).*T); %forward asset price
d1=(log(forward./K)+realVol.*realVol/2)./realVol;
  d2=d1-realVol;
 price=exp(-r*T).*(forward.*normcdf(d1)-K.*normcdf(d2)) ;
end