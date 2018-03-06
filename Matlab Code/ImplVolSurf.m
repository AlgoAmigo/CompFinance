%% Implied Volatility in a stochastic volatility model
function [surf, prices] = ImplVolSurf(S0, K, T, r, q, params, model)
%Creating the implied volatility surface for a prespecified model.

%simN: number of simulations
%stepsN: number of time steps
%S0: starting price
%T: time horizon of simulation
%r: riskfree rate
%q: dividend yield
%params: parameter of the model
%model: the used model, at the moment only 'HESTON' is available

% preallocating for speed
prices = zeros(length(K),length(T));
surf = zeros(length(K),length(T));

% Determine the model
switch(model)
    case 'HESTON'
        % The Transform is regular for im>1, so 1.1 is a reasonable choice
        im = 1.1;
        for j=1:length(T)
            for i=1:length(K)
                % calculating the price for different strikes K and
                % maturities T with Lewis Method
                prices(i,j) = CallPricingLewis(S0, K(i),T(j),r,q, params, im,  model);
                surf(i,j) = BSimpVolNewton(S0, K(i), T(j), r, prices(i,j), q);
            end
        end
    otherwise
        error('This function does not support the passed variable.')
end
end