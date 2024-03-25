function v = MD(X, T, P, Vp, x_H2, x_CH4, Cratio)

% if X < 0.35
%     v = 0;
% else

if Cratio > 0.5
    v = 1e-10;

elseif X < 0.1
    v = 1e-10;

elseif Cratio < 1e-6
    v = 1e-10;
else
    a_c = exp(2300/T - 0.92 + (3860/T)*Cratio + log(Cratio/(1-Cratio)) );
    
    P = P/101325; %atm, convert from Pa
    
    P_H2 = x_H2 * P;
    P_CH4 = x_CH4 * P;
    
    R = 8.314;
   
    G = 89658.88 - 102.27*T - 0.00428*T^2 - 2499358.99*T^(-1); %thermocatalytic hydrogen production through decompositon of methane - a review
    % G = G + 26485.7;
    % Kgraph = 5e-7*exp(0.0177*T);
    % Ggraph = 8.314*T*log(Kgraph);
    Ggraph = 21950.5; % estimated as graphite Del G from Equilibria of decomposition reactions of carbon monoxide and methane over nickel catalysts
    %Ggraph = 25000;
    %Ggraph = 0;

    G = G + Ggraph;
    %5708.51; % J/mol

    %G = 26694 - 24.77*T; % Theeffectofmethanedecompositiononthe formationandmagneticpropertiesof ironcarbide preparedfromoolitichematite
    Keq = exp(-G/(R*T))*P;
    
    k = 16250*exp(-55000/(R*T));

    %v = ((k/(P_H2^0.5)) * (P_CH4 - P_H2^2*a_c/Keq)*(1-eps_bed));
    v = Vp*((k/(P_H2^0.5)) * (P_CH4 - P_H2^2*a_c/Keq));

end
% end