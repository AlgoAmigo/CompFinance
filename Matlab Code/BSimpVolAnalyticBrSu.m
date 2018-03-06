function implVol =BSimpVolAnalyticBrSu(S0 , T,marketPrice) 
%Brenner and Subrahmanyan ATM implied volatility approximation

%S0: spot price
%T: maturity
%marketPrice: Market Option Price

%Author: Lorenzo Torricelli

implVol = marketPrice./S0.*sqrt(2*pi./T)  ; 

return;
end