%This script generates synthetic prices for options and reprices them in a Monte
%Carlo Euler simulator.

rng(123, 'twister'); %seeding the Mersenne Twister

r=0.05; %risk free rate
S0=100; %underlying price 1
q=0;

n=40; %no. of strike steps

Tmin=0.5;
Tmax=3; %maximum value of time axis
Kmin=60;
Kmax=140; %maximum value of strike price axis

Handle = 'NONE';  %How to handle S>Kmax and S<Kmin

dt=(Tmax-Tmin)/(n-1);
T=(Tmin:dt:Tmax)';
dk=(Kmax-Kmin)/(n-1);
K=(Kmin:dk:Kmax)';

%%Local volatility surface from ivols
[Ktab, Ttab]=meshgrid(K, T);
%%Implied volatility surface
alpha=0.33;
ImpliedVolSurface=ones(n).*0.15 + 0.15.*(ones(n).*0.5 + alpha*Ttab).*(.5*Ktab./S0 - 1.5*ones(n)).^2./ ... 
    (2*(Ktab./S0)+ones(n)) ;  %%Synthetically generated surface from the last maturity and last stirke to first maturity and first stirke
ImpliedVolSurface=flipud(ImpliedVolSurface); %%Just to read the synthetic data in the right order

%Computing prices:
PriceSurface = BSprice(S0, Ktab, r, Ttab, ImpliedVolSurface, q);

%Now compute the local volatility surface.
locSurf = locVolSurf(K,T,r,q,ImpliedVolSurface,'BBF',S0);

simN = 100000;
SimPrices = BSEulerMod(simN,dt,S0,K,T,r,q,locSurf,Handle);



simPriceSurf = zeros(length(T)-3,length(K)-3);
for i = 3:length(T)-1
    for j = 2:length(K)-2
        simPriceSurf(i-2,j-1) = CallPutPricer(S0,SimPrices(:,i-1),K(j),T(i),r,q);
    end
end

RelError = abs(simPriceSurf-PriceSurface(3:39,2:38))./PriceSurface(3:39,2:38);