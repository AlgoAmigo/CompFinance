function price = BSEuler(simN, dt, S0, T, r, q, sigmaMatr)
%Euler scheme for the Black Scholes model

%simN: number of simulations
% stepsN: number of time steps
%S0: starting price
%T: time horizon of simulation
%r: riskfree rate
%q: dividend yield
%sigma: volatility

%Author: Lorenzo Torricelli

%rng(123, 'twister'); %seeding the Mersenne Twister

price=S0*ones(simN, length(T)+1);

drift=(r-q)*dt*ones(simN, length(T)); %at time 0 there is no evolution
dWt=normrnd(0,1,simN, length(T));
diffusion=dWt*sqrt(dt);
for i=2:length(T)+1
    [~,J] = find(abs(price(:, i-1)-K(2:38)')==min(abs(price(:, i-1)-K(2:38)'),[],2),simN,'first');
    price(:, i)=price(:, i-1)+ price(:, i-1).*(drift(:,i-1)+K(J).*diffusion(:,i-1));
end
return