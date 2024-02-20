function v = SMR(X, T, P, Vp, x_H2, x_H2O, x_CO, x_CH4)

P = P/101325; %atm, convert from Pa
%Vp = Vp*100^3; %convert to cm^3

% x_H2 = x_H2 * x_normalizer;
% x_H2O = x_H2O * x_normalizer;
% x_CO = x_CO * x_normalizer;
% x_CH4 = x_CH4 * x_normalizer;

P_H2 = x_H2 * P;
P_H2O = x_H2O * P;
P_CO = x_CO * P;
P_CH4 = x_CH4 * P;

R = 8.314;

%Keq = exp(23700/T - 8.424*log(T) - 1.571)

% Propertiesofthermodynamicequilibrium-based methaneautothermalreformingtogenerate hydrogen
Keq  = 1/exp(-22790*T^(-1) + 8.156*log(T) -  4.421e3*T^(-2) - 26.030);


if X < 0.25
    kf = 1e-10; %mol/s-m^3-atm^4
    %kf = 1.15e-9*exp(4.18/(R*1e-3*T))*100^3; %mol/s-m^3-atm^4
    %kf = 392*exp(6770/(R*T))*100^3;

else

    kf = 1.15e-9*exp(4.18/(R*1e-3*T))*100^3; %mol/s-m^3-atm^4
    %kf = 1.15e-8*exp(4.18/(R*1e-3*T))*100^3; %mol/s-m^3-atm^4

    %kf = 1.15e-9*exp(3/(R*1e-3*T))*100^3; %mol/s-m^3-atm^4
    %kf = 392*exp(6770/(R*T))*100^3;
end

% v = (Vp*kf*(P_CO*P_H2^3 - (P_CH4*P_H2O)/Keq))/1e13;
%kf = 1.15e-9*exp(4.18/(R*1e-3*T))*100^3 %mol/s-m^3-atm^4
v = (Vp*kf*(P_CO*P_H2^3 - (P_CH4*P_H2O)/Keq));

%kf = 392*exp(6770/(R*T))*100^3
%v = -kf *(P_CH4*P_H2O - (P_CO*P_H2^3)/(1/Keq)) * (1-0.5);



