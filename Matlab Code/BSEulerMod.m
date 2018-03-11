function price = BSEulerMod(simN, dt, S0,K, T, r, q, locSurf, Handle)
%Euler scheme for the Black Scholes model

%simN: number of simulations
%stepsN: number of time steps
%S0: starting price
%T: time horizon of simulation
%r: riskfree rate
%q: dividend yield
%sigma: volatility

%Author: Lorenzo Torricelli

%rng(123, 'twister'); %seeding the Mersenne Twister


%Create first stock price
price=S0*ones(simN, length(T)-2);
drift=(r-q)*T(3)*ones(simN,1); %at time 0 there is no evolution
dWt=normrnd(0,1,simN, 1);
Kind = find(abs(S0-K(2:end-2))==min(abs(S0-K(2:end-2))),1,'first');
sigma = locSurf(1,Kind);
diffusion=sigma.*sqrt(T(3)).*dWt;
price(:, 2)=price(:, 1)+ price(:, 1).*(drift+diffusion);

if strcmp(Handle, 'TRUNC')==1
    price(price(:,2)>150,2) = 150;
    price(price(:,2)<50,2) = 50;
end

%Create other stock prices
drift=(r-q)*dt*ones(simN,1);
for i = 3:length(T)-2
    [I,~] = find((abs(price(:, i-1)-K(2:end-2)')-min(abs(price(:, i-1)-K(2:end-2)'),[],2))'==0,simN,'first');
    sigma = locSurf(i-1,I)';
    dWt=normrnd(0,1,simN, 1);
    diffusion=sigma.*sqrt(dt).*dWt;
    price(:, i)=price(:, i-1)+ price(:, i-1).*(drift+diffusion);
    if strcmp(Handle, 'TRUNC')==1
        price(price(:,i)>170,i) = 170;
        price(price(:,i)<40,i) = 40;
    end
end

if strcmp(Handle,'STOPS')==1
    [I,~] = find(price>150|price<50);
    REM = unique(I);
    IND = setdiff(1:simN,unique(REM));
else
    IND = 1:simN;
end
price = price(IND,:);

return