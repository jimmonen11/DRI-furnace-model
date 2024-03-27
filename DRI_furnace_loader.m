
changenodes = false;

H2only = false;
NGstart = true;

n_furnace = 100;

% Geometry of furnace, pellets, etc. come from data from a plant in Quebec taken from 
% Hamadeh, H. (2017). Modélisation mathématique détaillée du procédé de réduction directe du minerai de fer. 1–143. https://theses.hal.science/tel-01740462v1/document

rho_p = 3528; % kg/m^3 density of a pellet
eps_bed = 0.5; % m^3 gas/m^3 total

r_furnace = 2.75; %m, radius of reducing section of furnace
h_furnace = 10; %m, height of reducing section of furnace
r_p = 15e-3/2; %m, radius of iron ore pellets

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

R = 8.314; %(m^3*Pa)/(K*mol)

Solids_In_Flow = 45.1;
Solids_In_Flow_Step = 45.1;

T_sin = -10 + 273; %K, solids temperature in
P_gin = 101325*1.8; %Pa, pressure of gas in


w_Fe2O3in = 0.9665;
w_Fe3O4in = 0;
w_FeOin = 0;
w_Fein = 0;
w_Cin = 0;
w_Ganin = 0.0335;


if H2only == true
    load('initcond_H2.mat')
    %load('initcond_H2-SOE.mat')

    h_furnace = 5;
    
    x_CH4in = 0.0;
    x_H2in = 0.963133748;
    x_COin = 0.0;
    x_H2Oin = 0.019591391;
    x_CO2in = 0.0;
    x_N2in = 0.017274862;
    Gas_In_Flow = 2241.614*(x_H2in*MM_H2 + x_H2Oin*MM_H2O + x_N2in*MM_N2)/1000; %kg/s
    T_gin = 947 + 273; %K, temperature of gas in

     
    % x_CH4in = 0.0;
    % x_H2in = 0.9598;
    % x_COin = 0.0;
    % x_H2Oin = 0.0196;
    % x_CO2in = 0.0;
    % x_N2in = 0.0206;
    % Gas_In_Flow = 1880.858*(x_H2in*MM_H2 + x_H2Oin*MM_H2O + x_N2in*MM_N2)/1000; %kg/s
    % T_gin = 947 + 273; %K, temperature of gas in

   
    
elseif NGstart == true
    load('initcond_NG.mat')
    
    % 0.10 slip 
    x_CH4in = 0.105824828;
    x_H2in = 0.489526511;
    x_COin = 0.320711862;
    x_H2Oin = 0.042557156;
    x_CO2in = 0.024;
    x_N2in = 0.017379642;
    Gas_In_Flow = 2228.1*(x_H2in*MM_H2 + x_COin*MM_CO + x_H2Oin*MM_H2O +x_CO2in*MM_CO2 + x_N2in*MM_N2 + x_CH4in*MM_CH4)/1000; %kg/s
    T_gin = 947 + 273; %K, temperature of gas in

    % 0.05 slip
    % x_CH4in = 0.098397835;
    % x_H2in = 0.493566787;
    % x_COin = 0.323869613;
    % x_H2Oin = 0.042677198;
    % x_CO2in = 0.024;
    % x_N2in = 0.017488568;
    % Gas_In_Flow = 2203.05*(x_H2in*MM_H2 + x_COin*MM_CO + x_H2Oin*MM_H2O +x_CO2in*MM_CO2 + x_N2in*MM_N2 + x_CH4in*MM_CH4)/1000; %kg/s
    % T_gin = 952 + 273; %K, temperature of gas in



else
    load('initcond_NG-H2.mat')

    x_CH4in = 0.11138;
    x_H2in = 0.76842;
    x_COin = 0.08439;
    x_H2Oin = 0.01120;
    x_CO2in = 0.00632;
    x_N2in = 0.01829;
    Gas_In_Flow = 14.16536754;
    T_gin = 1047+273;

end

if changenodes == true

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


%Reducerflow = 2178*(x_H2in*MM_H2 + x_COin*MM_CO + x_H2Oin*MM_H2O +x_CO2in*MM_CO2 + x_N2in*MM_N2 + x_CH4in*MM_CH4)/1000;

% 25% with CH4/N2 still -  GOOD!
x_CH4step = 0.11138;
x_H2step = 0.76842;
x_COstep = 0.08439;
x_H2Ostep = 0.01120;
x_CO2step = 0.00632;
x_N2step = 0.01829;
Gas_In_Flow_Step = 14.16536754;
T_ginstep = 947+273;

% 50% with CH4/N2 still -  GOOD!
% x_CH4step = 0.10947;
% x_H2step = 0.67225;
% x_COstep = 0.16588;
% x_H2Ostep = 0.02201;
% x_CO2step = 0.01241;
% x_N2step = 0.01798;
% %Gas_In_Flow_Step = 14.16536754;
% Gas_In_Flow_Step = 2154.0*(x_H2in*MM_H2 + x_COin*MM_CO + x_H2Oin*MM_H2O +x_CO2in*MM_CO2 + x_N2in*MM_N2 + x_CH4in*MM_CH4)/1000; %kg/s
% T_ginstep = 947+273;

% Step back to normal NG operation!
% x_CH4step = 0.105824828;
% x_H2step = 0.489526511;
% x_COstep = 0.320711862;
% x_H2Ostep = 0.042557156;
% x_CO2step = 0.024;
% x_N2step = 0.017379642;
% Gas_In_Flow_Step = 2228.1*(x_H2step*MM_H2 + x_COstep*MM_CO + x_H2Ostep*MM_H2O +x_CO2step*MM_CO2 + x_N2step*MM_N2 + x_CH4step*MM_CH4)/1000; %kg/s
% T_ginstep = 947 + 273; %K, temperature of gas in


% 50 with CH4/N2 still
% x_CH4step = 0.10947;
% x_H2step = 0.67225;
% x_COstep = 0.16588;
% x_H2Ostep = 0.02201;
% x_CO2step = 0.01241;
% x_N2step = 0.01798;
% Gas_In_Flow_Step = 19.82445819;
% Solids_In_Flow_Step = Solids_In_Flow;


% SOE step
% x_CH4step = 0.00;
% x_H2step = 0.882955;
% x_COstep = 0.00;
% x_H2Ostep = 0.098106;
% x_CO2step = 0.0;
% x_N2step = 0.018939;
% Gas_In_Flow_Step = 3293.942*(x_H2step*MM_H2 + x_COstep*MM_CO + x_H2Ostep*MM_H2O +x_CO2step*MM_CO2 + x_N2step*MM_N2 + x_CH4step*MM_CH4)/1000; %kg/s
% T_ginstep = T_gin; %K, temperature of gas in


%
% x_CH4step = x_CH4in;
% x_H2step = x_H2in;
% x_COstep = x_COin;
% x_H2Ostep = x_H2Oin;
% x_CO2step = x_CO2in;
% x_N2step = x_N2in;
% Gas_In_Flow_Step = 12.5;
% T_ginstep = 800+273;

% H2 step
x_CH4step = x_CH4in;
x_H2step = x_H2in;
x_COstep = x_COin;
x_H2Ostep = x_H2Oin;
x_CO2step = x_CO2in;
x_N2step = x_N2in;
Gas_In_Flow_Step = Gas_In_Flow*1;
Solids_In_Flow_Step = Solids_In_Flow*1;
T_ginstep = T_gin;


% zpts = 1:1:n_furnace+2;
% b = log(h_furnace+1)/(n_furnace+1);
% a = 1/exp(b);
% dzfun = a*exp(b.*zpts)-1;
% dz = flip(diff(dzfun));


tau = 60*5; % seconds, time constant for 1st order step changes

A_furnace = pi*r_furnace^2; %m^2, c.s. area of flow
A_furnace_pel = A_furnace*(1-eps_bed); %c.s. area for pellet flow - excludes gas lanes
A_furnace_gas = A_furnace*(eps_bed); %c.s. area for gas flow - excludes solids lanes

V_p = 4/3*pi*r_p^3; %m^3, volume of a pellet
V_furnace = A_furnace*h_furnace; %m^3, volume of reducing section of furnace
V_pellet_bed = ((4/3)*pi*r_p^3)/(1-eps_bed); %m^3, volume of pellets in furnace

dz = ones(1, n_furnace+1)* (h_furnace/(n_furnace+1)); %m, spacing of nodes
n_pellets = V_furnace/V_pellet_bed; % no. of pellets in reducing section
n_pellets_dz = n_pellets.*dz/h_furnace; % no. of pellets per node

% Volume of gas and solid in each spatial node - constant
V_g = (4/3)*pi*r_p^3*n_pellets_dz*(eps_bed/(1-eps_bed)); % m^3, volume of gas
V_s = (4/3)*pi*r_p^3*n_pellets_dz; % m^3, volume of solid

a_b = 6*(1-eps_bed)/(r_p*2); %m^2/m^3, suface area for gas solid heat exchange, Wagner thesis
A_wall = (2*pi*r_furnace*dz); %m^2 surface area of furnace


V_ref = 193*1000/3600;

ndot_ref = 101325*V_ref/(R*273.15)

x_H2val = 0.4028;
x_COval = 0.1958;
x_H2Oval = 0.1903;
x_CO2val = 0.1709;
x_CH4val = 0.0295;
x_N2val = 0.0102;

ndot_ref*(x_H2val*MM_H2 + x_COval*MM_CO + x_H2Oval*MM_H2O +x_CO2val*MM_CO2 + x_N2val*MM_N2 + x_CH4val*MM_CH4)/1000 %kg/s