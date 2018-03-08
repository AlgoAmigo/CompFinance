%This function generates synthetic prices for options and reprices them in a Monte
%Carlo Euler simulator.

r=0.05; %risk free rate
S0=100; %underlying price 1
q=0;

n=40; %no. of strike steps

Tmin=0.5;
Tmax=3; %maximum value of time axis
Kmin=60;
Kmax=140; %maximum value of strike price axis


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
PriceSurface=BSprice(S0, Ktab, r, Ttab, ImpliedVolSurface, q);
%Now compute the local volatility surface.
locSurf = locVolSurf(K,T,r,q,ImpliedVolSurface,'BBF',S0);

%Create first stock price
simN = 10000;
Kind = find(abs(S0-K(2:38))==min(abs(S0-K(2:38))),1,'first');
sigma = locSurf(1,Kind);
simP = BSEuler(simN, dt, S0, T(2:38, r, q, locSurf);