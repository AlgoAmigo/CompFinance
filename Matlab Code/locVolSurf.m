%% Local Volatility in a stochastic volatility model
function surf = locVolSurf(K, T, r, q, inputData, method, S0)
% This function computes the local volatility surface by using either
% Dupires formula or the BFF formula. (Problem 4b)

%If the Dupire-formula is used, the following parameters are used:
%   K: strike price
%   T: time horizon of simulation
%   r: riskfree rate
%   q: dividend yield
%   inputData: price-matrix
%If the BFF formula is used, the additional parameter is:
%   S0: starting price

%method: the used method. Either BFF or DUPIRE.

switch method
    case 'DUPIRE'
        surf = localVolFromPrices(K,T, r , q, inputData);
    case 'BFF'
        if(~(nargin == 7))
            error('The amount of passed parameters is not correct.')
        end
        surf = localVolFromImplVol(S0, K, T, r, q, inputData);
    otherwise
        error('The passed method parameter is not valid. It has to be either DUPIRE or BFF.')
end
end