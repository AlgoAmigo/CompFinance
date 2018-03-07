function localVol=localVolFromPrices(K,T, r , q, marketValues)

%%Dupire Local Volatility functions from matrix of prices. Uses central
%%differences for partial derivatives calculation. Prices and maturities
%%are assumed to be equally spaced

%K: Strike
%T: Maturity
%r: riskfree rate
%q dividend yield
%marketValues: option prices

%Author: Lorenzo Torricelli

numStrikes=length(K);
numMaturities=length(T);
localVol=zeros(numMaturities-1, numStrikes-1);

%%Initialise the Log-forward-moneyness variable

deltaT=T(2)-T(1);
deltaK=K(2)-K(1);


for i=2:size(localVol,1)-1
    for j=2:size(localVol,2)-1
        timeDerivative=(marketValues(i+1, j)- marketValues(i-1, j))/ (2*deltaT) ; %%Central differnces  for time derivative
        
        strikeDerivative=(marketValues(i, j+1)- marketValues(i, j-1))/ (2*deltaK);  %%Central differnces  for strike derivative
       
        strikeSecondDerivative=(marketValues(i, j+1)- 2*marketValues(i, j) + marketValues(i, j-1))/ (deltaK*deltaK) ; %%Central differnces  for strike derivative
        
         

   %Dupire formula
        localVol(i, j)=sqrt((timeDerivative+K(j)*(r-q)*strikeDerivative ...
            +q*marketValues(i,j) ) / (strikeSecondDerivative*K(j)*K(j)) )*sqrt(2);
    end
end
localVol=localVol(2:size(localVol,1)-1 , (2:size(localVol,2)-1) ) ;

return