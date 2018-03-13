%% Local Volatility in a stochastic volatility model
function surf = locVolSurf(K, T, r, q, inputData, method, S0)
% This function computes the local volatility surface by using either
% Dupires formula or the BBF formula.

%If the Dupire-formula is used, the following parameters are used:
%K: strike price
%T: time horizon
%r: riskfree rate
%q: dividend yield
%inputData: price-matrix or implied volatility matrix

%If the BBF formula is used, the additional parameter is:
%S0: starting price

%method: the used method. Either BBF or DUPIRE.

%Author: Aaron Wittmann and Jan Keesen

switch method
    case 'DUPIRE'
        surf = localVolFromPrices(K,T, r , q, inputData);
    case 'BBF'
        if(~(nargin == 7))
            error('The amount of passed parameters is not correct.')
        end
        surf = localVolFromImplVol(S0, K, T, r, q, inputData);
    otherwise
        error('The passed method parameter is not valid. It has to be either DUPIRE or BBF.')
end
end