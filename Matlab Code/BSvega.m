function vega= BSvega(S0, K, r, T, sigma, q)
%Black scholes vega:
%Author: Lorenzo Torricelli

%S0: spot price
%r: riskfree rate
%sigma: volatility
%q: dividend yield
%T: time-to-maturity

 realVol=sigma.*sqrt(T);  %realised volatiliy variable
 forward=S0.*exp((r-q)*T);  %forward price
d1=(log(forward./K)+(sigma.*sigma/2).*T)./realVol;
 vega=exp(-q*T).*sqrt(T)*S0.*normpdf(d1);
end