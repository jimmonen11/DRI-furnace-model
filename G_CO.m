S_CO = 197.66;
S_CO2 = 213.79;


S_Fe2O3 = 87.2;
S_Fe3O4 = 146.4;
S_FeO = 60.75;
S_Fe = 27.3;

S1 = (2*S_Fe3O4 + S_CO2) - (3*S_Fe2O3 + S_CO);
S2 = (3*S_FeO + S_CO2) - (S_Fe3O4 + S_CO);
S3 = (S_Fe + S_CO2) - (S_FeO + S_CO);

R = 8.314;

% 
% H1(i) = -52.87e3;
% H2(i) = 36.25e3;
% H3(i) = -17.31e3;

% H1 = -12636*4.184;
% H2 = 8664*4.184;
% H3 = -4136*4.184;

T = 570:1:1400;
% 
% for i = 1:length(T)
%     [H1(i), H2(i), H3(i)] = heat_rxn_CO(T(i));
% end

% 
% H_CO2 = py.CoolProp.CoolProp.PropsSI('Hmolar','T',T,'P',101325,"CO2")...
%     - py.CoolProp.CoolProp.PropsSI('Hmolar','T',298,'P',101325,"CO2")+ (-393509);
% 
% H_CO = py.CoolProp.CoolProp.PropsSI('Hmolar','T',T,'P',101325,"CO")...
%     - py.CoolProp.CoolProp.PropsSI('Hmolar','T',298,'P',101325,"CO") + (-110525);
% 
% 
% [~, H_irons] = iron_props(T);
% 
% H_Fe2O3  = H_irons(1);
% H_Fe3O4  = H_irons(2);
% H_FeO  = H_irons(3);
% H_Fe  = H_irons(4);
% 
% % Overall
% % 3Fe2O3 + 9H2 ====> 6Fe + 9H2O 
% %  Fe2O3 + 3H2 ====> 2Fe + 3H2O 
% delH = ((2*H_Fe + 3*H_CO2) - (H_Fe2O3 + 3*H_CO))/3;
% 
% % First
% % 3Fe2O3 + H2 ====> 2Fe3O4 + H2O
% delH1 = ((2*H_Fe3O4 + H_CO2) - (3*H_Fe2O3 + H_CO));
% 
% 
% % Second
% % 2Fe3O4 + 2H2 ====> 6FeO + 2H2O
% %  Fe3O4 +  H2 ====> 3FeO + H2O
% 
% 
% delH2 = ((3*H_FeO + H_CO2) - (H_Fe3O4 + H_CO));
% 
% 
% % Third
% % 6FeO + 6H2 ====> 6Fe + 6H2O
% delH3 = ((H_Fe + H_CO2) - (H_FeO + H_CO));

% G1 = -96.744*(T) - 33259;
% G2 = 1e-5.*T.^3 - 0.0173.*T.^2 + 5.104.*T + 4825.9;
% G3 = 29.737*(T) - 15088;

G1 = -52870 - 47.33*(T);
G2 = 36250 - 51.98*(T);
G3 = -17310 + 17.32*(T);

% 
% G1 =  -47184 - 47.33*(T);
% G2 =  19416 - 51.98*(T);
% G3 = -10984  + 17.32*(T);


% K1 = exp(-G1./(R*(T)));
% K2 = exp(-G2./(R*(T)));
% K3 = exp(-G3./(R*(T)));

%K2 = exp(-1373./T  - 0.3411*log(T) + 0.41e-3.*T+2.303);
%K3 = exp(381./T  - 2.110*log(T) + 0.395e-3.*T+5.373);

% K1 = exp((362/T)+10.32);
% K2 = exp((-8580/T)+8.98);
% K3 = exp((-2070/T)+1.30);

K1 = exp((5128.6./T)+5.7);
K2 = exp((-3132.5./T)+3.661);
K3 = exp((2240.6./T)-2.667);


x_COeq1 = (1./(K1+1));
x_COeq2 = (1./(K2+1));
x_COeq3 = (1./(K3+1));

close all
plot(T-273, x_COeq1)
hold on
plot(T-273, x_COeq2)
plot(T-273, x_COeq3)
legend()
