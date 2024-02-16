function [Cps, Hs] = iron_props(T)

% Molar mass of species in kg/mol
MM_Fe2O3 = 159.69/1000;
MM_Fe3O4 = 251.53/1000;
MM_FeO = 71.844/1000;
MM_Fe = 55.845/1000;


% From Perry's Handbook in cal/mol-K converted to J/kg-K
Cp_Fe2O3 = (24.72 + 0.01604*T - 423400/T^2)*4.184/MM_Fe2O3; %
Cp_Fe3O4 = (41.17 + 0.01882*T - 979500/T^2)*4.184/MM_Fe3O4;
Cp_FeO = (12.62 + 0.001492*T - 76200/T^2)*4.184/MM_FeO;
Cp_Fe = (4.13 + 0.00638*T)*4.184/MM_Fe;

Cp_Gan = Cp_Fe2O3; %estimate
Cp_C = 890; %estimate from NIST, doens't change with T


% enthalpy's - enthalpy of fomfation from Smith, Van Ness et al.
% Introduction to Chemical Engineering Thermodynamics
% J/mol

T0 = 298; 

H_Fe2O3 = ( (24.72*T + (0.01604/2)*T^2 + 423400/T) - (24.72*T0 + (0.01604/2)*T0^2 + 423400/T0))*4.184 + (-824200); %825500
H_Fe3O4 = ( (41.17*T + (0.01882/2)*T^2 + 979500/T) - (41.17*T0 + (0.01882/2)*T0^2 + 979500/T0) )*4.184 + (-1118400); %1120890
H_FeO = ( (12.62*T + (0.001492/2)*T^2 + 76200/T) - (12.62*T0 + (0.001492/2)*T0^2 + 76200/T0) )*4.184 + (-272000);
H_Fe = ( (4.13*T + (0.00638/2)*T^2) - (4.13*T0 + (0.00638/2)*T0^2) )*4.184;

% S_Fe2O3 = 87.2;
% S_Fe3O4 = 146.4;
% S_FeO = 60.75;
% S_Fe = 27.3;

S_Fe2O3 = ( (24.72*log(T) + 0.01604*T - 423400/(2*T^2)) - (24.72*log(T0) + 0.01604*T0 - 423400/(2*T0^2)) )*4.184 + (87.2);
S_Fe3O4 = ( (41.17*log(T) + 0.01882*T - 979500/(2*T^2) ) - (41.17*log(T0) + 0.01882*T0 - 979500/(2*T0^2)) )*4.184 + (146.4);
S_FeO = ( (12.62*log(T) + 0.001492*T - 76200/(2*T^2)) - (12.62*log(T0) + 0.001492*T0 - 76200/(2*T0^2)) )*4.184 + (60.75);
S_Fe = ( (4.13*log(T) + 0.00638*T ) - (4.13*log(T0) + 0.00638*T0) )*4.184 + (27.3);



Cps = [Cp_Fe2O3, Cp_Fe3O4, Cp_FeO, Cp_Fe, Cp_C, Cp_Gan]; %J/kg-K
Hs = [H_Fe2O3, H_Fe3O4, H_FeO, H_Fe]; %J/mol
Ss = [S_Fe2O3, S_Fe3O4, S_FeO, S_Fe]; %J/mol-K ]
