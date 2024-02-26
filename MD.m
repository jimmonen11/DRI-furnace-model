function v = MD(X, T, P, x_H2, x_CH4, Cratio, eps_bed)

% if X < 0.35
%     v = 0;
% else

if Cratio > 0.5
    v = 0;

elseif X < 0.1
    v = 0;

elseif Cratio < 1e-6
    v = 0;
else
    a_c = exp(2300/T - 0.92 + (3860/T)*Cratio + log(Cratio/(1-Cratio)) );
    
    P = P/101325; %atm, convert from Pa
    
    P_H2 = x_H2 * P;
    P_CH4 = x_CH4 * P;
    
    R = 8.314;
    
    G = 26694 - 24.77*T; % Theeffectofmethanedecompositiononthe formationandmagneticpropertiesof ironcarbide preparedfromoolitichematite
    Keq = exp(-G/(R*T));
    
    k = 16250*exp(-55000/(R*T));
    
    v = ((k/(P_H2^0.5)) * (P_CH4 - P_H2^2*a_c/Keq)*(1-eps_bed));
end
% end