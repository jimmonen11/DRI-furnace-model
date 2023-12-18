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


% enthalpy's - enthalpy of fomfation from Smith, Van Ness et al.
% Introduction to Chemical Engineering Thermodynamics
% J/mol

T0 = 298; 

H_Fe2O3 = ( (24.72*T + (0.01604/2)*T^2 + 423400/T) - (24.72*T0 + (0.01604/2)*T0^2 + 423400/T0))*4.184 + (-824200);
H_Fe3O4 = ( (41.17*T + (0.01882/2)*T^2 + 979500/T) - (41.17*T0 + (0.01882/2)*T0^2 + 979500/T0) )*4.184 + (-1118400);
H_FeO = ( (12.62*T + (0.001492/2)*T^2 + 76200/T) - (12.62*T0 + (0.001492/2)*T0^2 + 76200/T0) )*4.184 + (-272000);
H_Fe = ( (4.13*T + (0.00638/2)*T^2) - (4.13*T0 + (0.00638/2)*T0^2) )*4.184;


Cps = [Cp_Fe2O3, Cp_Fe3O4, Cp_FeO, Cp_Fe]; %J/mol-K
Hs = [H_Fe2O3, H_Fe3O4, H_FeO, H_Fe]; %J/mol