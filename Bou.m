function v = Bou(T, P, x_H2, x_CO, x_CO2, Cratio, eps_bed)

a_c = exp(2300/T - 0.92 + (3860/T)*Cratio + log(Cratio/(1-Cratio)) );

P = P/101325; %atm, convert from Pa

P_H2 = x_H2 * P;
P_CO = x_CO * P;
P_CO2 = x_CO2 * P;


R = 8.314;

Keq = 10^(9141/T + 0.000224*T-9.595); %from wikipedia bou

k1 = 1.8*exp(-27200/(R*T));
k2 = 2.2*exp(-8800/(R*T));

v = ((k1*(P_H2^0.5) + k2) * (P_CO^2 - P_CO2*a_c/Keq)*(1-eps_bed));