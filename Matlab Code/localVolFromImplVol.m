function localVol=localVolFromImplVol( S0, K,T, r , q, marketValues)

%%Local Volatility functions from matrix of implied volatilities. BBF formula

%%Uses central differences for the time derivative and general second
%%derivative for a discret function passing thorugh 3 unevenly spaced
%%points 


%K: Strike
%T: Maturity
%r: riskfree rate
%q dividend yield implied volatilities prices

%Author: Lorenzo Torricelli




numStrikes=length(K);
numMaturities=length(T);
localVol=zeros(numMaturities-1, numStrikes-1);

%%Initialise the Log-forward-moneyness variable
deltaT=T(2)-T(1);
deltaK=K(2)-K(1);


for i=2:size(localVol,1)-1
    for j=2:size(localVol,2)-1
        timeDerivative=(marketValues(i+1, j)- marketValues(i-1, j))/ (2*deltaT);  %%Central differnces  for time derivative
        
        strikeDerivative=(marketValues(i, j+1)- marketValues(i, j-1))/ (2*deltaK);  %%Central differnces  for strike derivative
         strikeSecondDerivative=(marketValues(i, j+1)- 2*marketValues(i, j) + marketValues(i, j-1))/ (deltaK*deltaK);  %%Central differnces  for strike derivative
        
         sigma=marketValues(i,j);
         
         forwardLogMoneynessAdjusted=(log(S0/K(j))+(r-q+sigma*sigma/2)*T(i))/sigma; %the d+ !!

         
         
         num=sigma*sigma + 2*T(i)*sigma*(timeDerivative + (r-q)*K(j)*strikeDerivative);
           denom=(1+K(j)*forwardLogMoneynessAdjusted*strikeDerivative)^2 + K(j)*K(j)*T(i)*sigma*sigma*(strikeSecondDerivative ...
          -forwardLogMoneynessAdjusted*strikeDerivative*strikeDerivative);
         
         %Local volatility formula from implied vols
        localVol(i, j)=sqrt(num/denom);
    end
end
localVol=localVol(2:size(localVol,1)-1 , (2:size(localVol,2)-1) ) ;
return