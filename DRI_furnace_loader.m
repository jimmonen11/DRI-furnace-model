%% Number of discretization modes and type of operation

n_furnace = 75;

changenodes = false; %set to true if plan on changing n_furnace to something other than 75
H2only = false; % set to true for a furnace that only uses H2

steady_state = true; % set to true if want constant inlet conditions to be constant
stepcase1 = false;
stepcase2 = false;
stepcase3 = false;

%% Furnace Geometry and Solids Parameters
% Geometry of furnace, pellets, etc. come from data from a plant in Quebec taken from 
% Hamadeh, H. (2017). Modélisation mathématique détaillée du procédé de réduction directe du minerai de fer. 1–143. https://theses.hal.science/tel-01740462v1/document
rho_p = 3528; % kg/m^3 density of a pellet
eps_bed = 0.5; % m^3 gas/m^3 total

r_furnace = 2.75; %m, radius of reducing section of furnace
h_furnace = 10; %m, height of reducing section of furnace
r_p = 15e-3/2; %m, radius of iron ore pellets

A_furnace = pi*r_furnace^2; %m^2, c.s. area of flow
A_furnace_pel = A_furnace*(1-eps_bed); %c.s. area for pellet flow - excludes gas lanes
A_furnace_gas = A_furnace*(eps_bed); %c.s. area for gas flow - excludes solids lanes

V_p = 4/3*pi*r_p^3; %m^3, volume of a pellet
V_furnace = A_furnace*h_furnace; %m^3, volume of reducing section of furnace
V_pellet_bed = ((4/3)*pi*r_p^3)/(1-eps_bed); %m^3, volume of pellets in furnace

dz = (h_furnace/(n_furnace-1)); %m, spacing of nodes

n_pellets = V_furnace/V_pellet_bed; % no. of pellets in reducing section
n_pellets_dz = n_pellets*dz/h_furnace; % no. of pellets per node

% Volume of gas and solid in each spatial node - constant
V_g = (4/3)*pi*r_p^3*n_pellets_dz*(eps_bed/(1-eps_bed)); % m^3, volume of gas
V_s = (4/3)*pi*r_p^3*n_pellets_dz; % m^3, volume of solid

a_b = 6*(1-eps_bed)/(r_p*2); %m^2/m^3, suface area for gas solid heat exchange, Wagner thesis
A_wall = (2*pi*r_furnace*dz); %m^2 surface area of furnace

tau = 60*5; % seconds, time constant for 1st order step changes

%% Constants

% Molar mass of species in kg/mol
MM_Fe2O3 = 159.69/1000;
MM_Fe3O4 = 231.53/1000;
MM_FeO = 71.844/1000;
MM_Fe = 55.845/1000;
MM_C = 12.011/1000; 
MM_Gan = 60.08/1000;

% Molar mass of species in g/mol
MM_H2O = 18.015;
MM_H2 = 2.016;
MM_N2 = 28.013;
MM_CO = 28.01;
MM_CO2 = 44.01;
MM_CH4 = 16.04; 

% gas diffusion properties - MM and diffusion volumes
H2O_diff_prop = [MM_H2O, 13.1];
H2_diff_prop = [MM_H2, 2.31*2];

CO_diff_prop = [MM_CO, 18.9];
CO2_diff_prop = [MM_CO2, 26.9];

% ideal gas constant
R = 8.314; %(m^3*Pa)/(K*mol)

%% Solid Inlet Conditions
Solids_In_Flow = 45.1; %kg/s, flow of solids in

T_sin = -10 + 273; %K, solids temperature in
P_gin = 101325*1.785; %Pa, pressure of gas in (found via trial and error)

% weight fraction of pellets in
w_Fe2O3in = 0.9665;
w_Fe3O4in = 0;
w_FeOin = 0;
w_Fein = 0;
w_Cin = 0;
w_Ganin = 0.0335;

%% Gas Inlet Conditions

if H2only == true
    load('initcond_H2.mat')
    %load('initcond_H2_10m.mat') % can load 10 m results if want but change h_furnace
    h_furnace = 5;
    
    x_CH4in = 0.0;
    x_H2in = 0.962059;
    x_COin = 0.0;
    x_H2Oin = 0.019314;
    x_CO2in = 0.0;
    x_N2in = 0.018627;
    Gas_In_Flow = 2078.874*(x_H2in*MM_H2 + x_H2Oin*MM_H2O + x_N2in*MM_N2)/1000; %kg/s
    T_gin = 942 + 273;

    
else
    load('initcond_NG.mat')
    
    % Steady state NG-DRI conditions
    x_CH4in = 0.105824828;
    x_H2in = 0.489526511;
    x_COin = 0.320711862;
    x_H2Oin = 0.042557156;
    x_CO2in = 0.024;
    x_N2in = 0.017379642;
    Gas_In_Flow = 2228.1*(x_H2in*MM_H2 + x_COin*MM_CO + x_H2Oin*MM_H2O +x_CO2in*MM_CO2 + x_N2in*MM_N2 + x_CH4in*MM_CH4)/1000; %kg/s
    T_gin = 942 + 273;


end

if changenodes == true %this 'helps' change nodes, but will still required sim time to sort itself out

    T_ginit = interp1([1:1:length(T_ginit)],T_ginit, linspace(1,length(T_ginit),n_furnace));
    T_sinit = interp1([1:1:length(T_sinit)],T_sinit, linspace(1,length(T_sinit),n_furnace));
    
    c_H2Oinit = interp1([1:1:length(c_H2Oinit)],c_H2Oinit, linspace(1,length(c_H2Oinit),n_furnace));
    c_H2init = interp1([1:1:length(c_H2init)],c_H2init, linspace(1,length(c_H2init),n_furnace));
    c_N2init = interp1([1:1:length(c_N2init)],c_N2init, linspace(1,length(c_N2init),n_furnace));
    c_COinit = interp1([1:1:length(c_COinit)],c_COinit, linspace(1,length(c_COinit),n_furnace));
    c_CO2init =interp1([1:1:length(c_CO2init)],c_CO2init, linspace(1,length(c_CO2init),n_furnace));
    c_CH4init =interp1([1:1:length(c_CH4init)],c_CH4init, linspace(1,length(c_CH4init),n_furnace));
    
    c_Feinit = interp1([1:1:length(c_Feinit)],c_Feinit, linspace(1,length(c_Feinit),n_furnace));
    c_FeOinit = interp1([1:1:length(c_FeOinit)],c_FeOinit, linspace(1,length(c_FeOinit),n_furnace));
    c_Fe3O4init = interp1([1:1:length(c_Fe3O4init)],c_Fe3O4init, linspace(1,length(c_Fe3O4init),n_furnace));
    c_Fe2O3init = interp1([1:1:length(c_Fe2O3init)],c_Fe2O3init, linspace(1,length(c_Fe2O3init),n_furnace));
    c_Cinit = interp1([1:1:length(c_Cinit)],c_Cinit, linspace(1,length(c_Cinit),n_furnace));
    
    
    nr1init = interp1([1:1:length(nr1init)],nr1init, linspace(1,length(nr1init),n_furnace));
    nr2init = interp1([1:1:length(nr2init)],nr2init, linspace(1,length(nr2init),n_furnace));
    nr3init = interp1([1:1:length(nr3init)],nr3init, linspace(1,length(nr3init),n_furnace));
    ndotinit = interp1([1:1:length(ndotinit)],ndotinit, linspace(1,length(ndotinit),n_furnace));

end


x_sumin = x_CH4in  + x_H2in + x_COin +x_H2Oin + x_CO2in + x_N2in; %check to make sure equal to 1


if steady_state == true
    
    %inlet values are the same throughout simulation
    x_CH4step = x_CH4in;
    x_H2step = x_H2in;
    x_COstep = x_COin;
    x_H2Ostep = x_H2Oin;
    x_CO2step = x_CO2in;
    x_N2step = x_N2in;
    Gas_In_Flow_Step = Gas_In_Flow*1;
    Solids_In_Flow_Step = Solids_In_Flow*1;
    T_ginstep = T_gin;

end

if stepcase1 == true

    % 25% turndown of syngas with CH4/N2
    x_CH4step = 0.11138;
    x_H2step = 0.76842;
    x_COstep = 0.08439;
    x_H2Ostep = 0.01120;
    x_CO2step = 0.00632;
    x_N2step = 0.01829;
    Gas_In_Flow_Step = 14.16536754;
    T_ginstep = 942+273;

end


if stepcase2 == true 
    
    %Move to 80% of flow
    x_CH4step = x_CH4in;
    x_H2step = x_H2in;
    x_COstep = x_COin;
    x_H2Ostep = x_H2Oin;
    x_CO2step = x_CO2in;
    x_N2step = x_N2in;
    Gas_In_Flow_Step = Gas_In_Flow*0.8;
    Solids_In_Flow_Step = Solids_In_Flow*1;
    T_ginstep = T_gin;

end

if stepcase3 == true 
    
    %Move to 80% of flow
    x_CH4step = x_CH4in;
    x_H2step = x_H2in;
    x_COstep = x_COin;
    x_H2Ostep = x_H2Oin;
    x_CO2step = x_CO2in;
    x_N2step = x_N2in;
    Gas_In_Flow_Step = Gas_In_Flow*1;
    Solids_In_Flow_Step = Solids_In_Flow*1;
    T_ginstep = T_gin-100;

end

x_sumstep = x_CH4step  + x_H2step + x_COstep +x_H2Ostep + x_CO2step + x_N2step; %check to make sure equal to 1

