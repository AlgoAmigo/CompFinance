function [surf, prices] = ImplVolSurf(S0, K, T, r, q, params, model)
%Creating the implied volatility surface for a prespecified model.

%S0: starting price
%K: strike of simulation
%T: time horizon of simulation
%r: riskfree rate
%q: dividend yield
%params: parameter of the model
%model: the used model, at the moment only 'HESTON' is available

%Author: Aaron Wittmann and Jan Keesen

% preallocating for speed
prices = zeros(length(T),length(K));
surf = zeros(length(T),length(K));

% Determine the model
switch(model)
    case 'HESTON'
        % The Transform is regular for im>1, so 1.1 is a reasonable choice
        im = 1.1;
        for i=1:length(T)
            for j=1:length(K)
                % calculating the price for different strikes K and
                % maturities T with Lewis Method
                prices(i,j) = CallPricingLewis(S0, K(j), T(i), r, q, params, im, model);
                surf(i,j) = BSimpVolFzero(S0, K(j), T(i), r, prices(i,j), q);
            end
        end
    otherwise
        error('This function does not support the passed variable.')
end
end