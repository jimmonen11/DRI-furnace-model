function [v, Keq] = SMR(X, T, P, Vp, x_H2, x_H2O, x_CO, x_CH4)


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

%Keq = 1/exp(23700/T - 8.424*log(T) - 1.571); %atm^-2
%Keq = (exp(30.53 - 4.8486e4/T +2.421748e6/T^2 + 2.49e9/T^3 ))

z = 1000/T - 1;
%Keq = ((1.0267e10* exp(-0.2513*z^4 + 0.3665*z^3 + 0.5810*z^2 - 27.134*z + 3.2770)/101325^2)) %atm^2

%Keq = ((1.0267e10* exp(-0.2513*z^4 + 0.3665*z^3 + 0.5810*z^2 - 27.134*z + 3.2770))/101325^2) %atm^2

Keq = (10266.76*10^6*exp (-26830/T + 30.11))/101325^2; %atm^2
%Keq = exp(-26830/T +30.114)  %/101325^2 %atm^2


%Keq  =exp(-26830/(T) + 30.114) 

% Propertiesofthermodynamicequilibrium-based methaneautothermalreformingtogenerate hydrogen
%Keq  = (exp(-22790/T + 8.156*log(T)  - 4.421e-3*T -  4.330e3*T^(-2) - 26.030)*P^2);

%Keq = exp((-125323-166.53*T)/(R*T));
% 
if X < 0.25
    kf = 1e-20; %mol/s-m^3-atm^4
    %kf = 1.15e-9*exp(4.18/(R*1e-3*T))*100^3; %mol/s-m^3-atm^4
    %kf = 392*exp(6770/(R*T))*100^3;
 
else

    %kf = 1.15e-9*exp(4.18/(R*1e-3*T))*100^3; %mol/s-m^3-atm^4
    %kf = 5e-5*exp(4.18/(R*1e-3*T))*100^3; %mol/s-m^3-atm^4

    %kf = 1e-4*exp(6.77/(R*1e-3*T))*100^3; %mol/s-m^3-atm^2
    %kf = 1.15e-7*exp(4.18/(R*1e-3*T))*100^3; %mol/s-m^3-atm^4

    %kf = 1.15e-3*exp(4.18/(R*1e-3*T))*100^3; %mol/s-m^3-atm^2
    %kf = 392*exp(-6.77/(R*1e-3*T))*100^3; %mol/s-m^3-atm^2
    %kb = 1.15e-9*exp(4.18/(R*1e-3*T))*100^3; %mol/s-m^3-atm^4

    %kf = 1.15e-9*exp(3/(R*1e-3*T))*100^3; %mol/s-m^3-atm^4
    %kf = 392*exp(6770/(R*T))*100^3;

    %kf = exp(21.55-258e3/(R*T));

    %kf = 4.225e16
    
kf = (2395*exp(-231266/(R*T)) ) * 101325^2;
    %kf = (2395*exp(-231266/(R*T)) ) * 101325^2;

    % Three-dimensional simulation of chemically reacting gas f lows in the porous support structure of an integrated-planar solid oxide fuel cell

end

% v = (Vp*kf*(P_CO*P_H2^3 - (P_CH4*P_H2O)/Keq))/1e13;
%kf = 1.15e-9*exp(4.18/(R*1e-3*T))*100^3 %mol/s-m^3-atm^4
%v = (Vp*kf*(P_CO*P_H2^3 - (P_CH4*P_H2O)/Keq));

%kf = 392*exp(6770/(R*T))*100^3
%v = Vp* kf * P^2*(x_CH4*x_H2O - (x_CO*x_H2^3*Keq));

%v = Vp* kf * P^2 * (x_CH4*x_H2O - ((x_CO*x_H2^3)*Keq));
%v = Vp* kf * P^2 * ((x_CH4*x_H2O/Keq) - x_CO*x_H2^3);

%v = (Vp*kf*(P_CH4*P_H2O - (P_CO*(P_H2^3))*Keq));

v = Vp* kf * (P_CH4*P_H2O - P_CO*P_H2^3/Keq);


%v = -Vp* kf * P^4 * ((x_CO*x_H2^3) - x_CH4*x_H2O/Keq);
%v = Vp* kf * P^2 * (x_CH4*x_H2O - x_CO*x_H2^3*Keq);

%m_Fe = rho_Fe*Vp*1000;

%v = m_Fe* kf * 1/(1+Keq*(P_H2O/P_H2))*(P_CH4/P_H2^0.5);
