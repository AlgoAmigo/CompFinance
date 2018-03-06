function cf = cfLibrary(om,S0,T,r,q,params, model)

%Characterisitc function of a jump diffusion to be used in various
%Fourier/Laplace valuation methods

%S0: spot price
%T: time to matruity
%r: riskfree rate
%q: dividend yield

%om: inversion variable

%Variance process parameters: kappa=mean reversion level; theta=long run
%variance; eta=volatility of variance, rho=correlation between the Brownian
%motions, lambda=intensity of jumps,  
%mu: average of normally distributed jump size in return
%delta: variance of normally distributed jump size

%Gilli and Schumann 2012. Modified by Lorenzo Torricelli

switch(model)

    
    
    case 'BATES'

        
v0=params(1); kappa = params (2); theta = params (3);
eta = params (4);  rho = params (5); lambda=params(6); mu=params(7);
delta=params(8);

        
        
d = sqrt( (rho * eta * 1i*om - kappa).^2 + ...
           eta^2 * (1i*om + om .^ 2) );
g = (kappa - rho*eta*1i*om - d) ./ (kappa - rho*eta*1i*om + d);
cf1 = 1i*om .* (log(S0) + (r - q) * T);
cf2 = theta*kappa / (eta^2) * ((kappa - rho*eta*1i*om - d) * T - ...
    2 * log((1 - g .* exp(-d * T)) ./ (1 - g)));
cf3 = v0/eta^2*(kappa-rho*eta*1i*om-d).*(1-exp(-d*T)) ./ ...
       (1-g.*exp(-d*T));
% jump
cf4 = -lambda*mu*1i*T*om + lambda*T* ...
      ((1+mu).^(1i*om) .* exp( delta*(1i*om/2) .* (1i*om-1) )-1);

cf  = exp(cf1 + cf2 + cf3 + cf4);


    case 'HESTON'

v0=params(1); kappa = params (2); theta = params (3);
eta = params (4);  rho = params (5); 
        
        
d = sqrt( (rho * eta * 1i*om - kappa).^2 + ...
           eta^2 * (1i*om + om .^ 2) );
       
g = (kappa - rho*eta*1i*om - d) ./ (kappa - rho*eta*1i*om + d);

cf1 = 1i*om .* (log(S0) + (r - q) * T);

cf2 = theta*kappa / (eta^2) * ((kappa - rho*eta*1i*om - d) * T - ...
    2 * log((1 - g .* exp(-d * T)) ./ (1 - g)));

cf3 = v0/eta^2*(kappa-rho*eta*1i*om-d).*(1-exp(-d*T)) ./ ...
       (1-g.*exp(-d*T));
% jump
cf  = exp(cf1 + cf2 + cf3);

        
        
    case 'MERTON'

        
v0=params(1);  lambda=params(2); mu=params(3);
delta=params(4);

        
        cf1= 1i*om*log(S0) + 1i*om*T*(r-q-0.5*v0-lambda*mu) ...
    - 0.5*(om.^2)*v0*T;
        
cf2 = lambda*T*(exp(1i*om*log(1+mu) ...
    -0.5*1i*om*delta-0.5*delta*om.^2) - 1); %%somewhat awkward definition by Bates
cf  = exp(cf1 + cf2);
  
        
    case 'BS'

v0=params(1);
        cf1= 1i*om*log(S0) + 1i*om*T*(r-q-0.5*v0) ...
    - 0.5*(om.^2)*v0*T;
        
cf  = exp(cf1);
  


 
    case 'NIG'

  kappa= params(1); sigma = params(2); mu = params(3);
   
        
cf1= (log(S0) + (r-q)*T).*1i.*om;
cf2=T.*(1-sqrt(1+om.^2.*sigma.^2.*kappa-2.*1i.*mu.*om.*kappa))./kappa;        
cf3=-T.*1i*om.*(1-sqrt(1 - sigma.*sigma.*kappa - 2.*mu.*kappa) )./kappa;
cf  = exp((cf1+cf2 +cf3));


    case 'VG'

    nu= params(1); sigma = params(2); theta = params(3);
    
    omega = (1/nu)*( log(1-theta*nu-sigma*sigma*nu/2) );
    tmp = 1 - 1i * theta * nu * om + 0.5 * sigma * sigma * om.* om * nu;
    %tmp = tmp.^(-T / nu);
    %phi = exp( i * u * (lnS + (r + omega - d) * T )) .* tmp;
    %phi = 1i * u * (lnS + (r + omega - d) * T ) - T*log(tmp)/nu;
    cf = exp(1i * om * (log(S0) +(r-q+omega)* T)  - T*log(tmp)/nu);
   


 case '32'
 cf=cf32(om,S0,T,r,q,params);

     otherwise 
        error('Characteristic function not available');
   
 
 end
 
return 