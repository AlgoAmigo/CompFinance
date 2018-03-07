function call = CallPricingLewis(S0, K,T,r,q, params, im, model)



%Fourier inversion method using Lewis approach

%S0: spot price
%K: strike price
%T: time to matruity
%r: riskfree rate
%q: dividend yield

%params: model params (see the cf .m file)
%model: stochastic model used
%trunc: value of numerical integral truncation

%%Bonus: im>1 call price: im<0 putprice

%Author: Lorenzo Torricelli

%%other possible implementation: 1) take out the constant K of thhe intergal:
%% 2) replace the intergal from -trunc to trunc with 2*Re(integrand)



i=1i;
integrand = @(x)cfLibrary(-(x+i*im),S0,T,r,q,params, model).*K.^(1+i.*(x+i*im))./(i.*(x+i*im)-(x+i*im).^2);

% Relative error of 0.1%
call = .5*exp(-r*T)/pi * integral(integrand,-Inf,Inf,'AbsTol',1e-6,'RelTol',1e-4);
call=abs(call);
 

end