function v = Bou(X, T, P, Vp, x_H2, x_CO, x_CO2, Cratio)


P = P/101325; %atm, convert from Pa

P_H2 = x_H2 * P;
P_CO = x_CO * P;
P_CO2 = x_CO2 * P;


R = 8.314;

% Modelling carbon formation using thermodynamic and kinetic methods in a steam methane reformer over nickel catalysts 
G = (-188030.19 + 402.82*T - 0.00524*T^2 + 828509.9/T - 32.026*T*log(T));
Keq = exp(-G/(R*T))/P;


if Cratio > 0.5
    v = 1e-10;

elseif Cratio < 1e-6
    v = 1e-10;

elseif X < 0.4
    v = 1e-20;

else
    a_c = exp(2300/T - 0.92 + (3860/T)*Cratio + log(Cratio/(1-Cratio)) );

    k1 = 1.8*exp(-27200/(R*T));
    k2 = 2.2*exp(-8800/(R*T));
    
        v = Vp * ( (k1*(P_H2^0.5) + k2) * (P_CO^2 - P_CO2*a_c/Keq) );
    
end