function v = WGS(X, T, P, Vp, x_H2, x_H2O, x_CO, x_CO2)

P = P/101325; %atm, convert from Pa
%Vp = Vp*100^3; %convert to cm^3

P_H2 = x_H2 * P;
P_H2O = x_H2O * P;
P_CO = x_CO * P;
P_CO2 = x_CO2 * P;


R = 8.314;

Keq = exp(3863.7/T - 3.5414);

Z = 1000/T - 1;
Keq = (exp(-0.2935*Z^3 + 0.6351*Z^2 + 4.1788*Z + 0.3169));

if X < 0.5
    kf = 1.83e-5*exp(7.84e-3/(R*1e-3*T))*100^3;
    %kf = 8.83e-2*exp(-4.56/(R*1e-3*T))*100^3;
else

    kf = 9.33e1*exp(-7.32/(R*1e-3*T))*100^3;
    
    %kf = (0.0171*exp(-103191/(R*T)) ) * 101325^2;


end


v = (Vp*kf*(P_CO*P_H2O - P_CO2*P_H2/Keq));

