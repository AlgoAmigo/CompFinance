T=[(7:7:70)]./365;
S0 = 2726;
K = 2705:5:2750;
q = 0.0179;
r = 0.0165;
%price(time,strike)
prices = [
    34.20,30.60,27.10,23.80,20.90,18.00,15.30,12.90,10.60,8.60;
    43.50,39.80,36.40,33.10,30.10,27.10,24.30,21.60,19.10,16.80;
    50.00,46.50,43.10,39.80,36.70,33.60,30.70,27.90,25.30,22.70;
    56.60,53.10,49.60,46.30,43.00,39.90,36.80,33.60,31.10,28.50;
    62.10,58.20,55.20,51.90,48.70,45.50,42.40,39.50,36.60,33.90;
    70.20,66.90,63.30,59.70,56.80,53.60,50.60,47.60,44.60,41.80;
    75.70,72.30,69.00,65.50,62.40,59.20,56.10,52.70,50.00,47.10;
    81.00,77.40,74.00,70.70,67.60,64.40,61.20,58.10,55.10,52.10;
    84.70,81.30,78.00,74.60,71.30,68.10,64.90,62.00,58.80,55.90;
    88.40,84.90,81.50,78.60,75.20,72.10,68.90,65.70,62.80,59.80
    ];