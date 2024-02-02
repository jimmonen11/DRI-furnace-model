n_furnace = 10;

eps_bed = 0.4; %S. Yu, L. Shao, Z. Zou, and H. Saxén, "A numerical study on the performance of the h2 shaft furnace with dual-row top gas recycling," Processes, vol. 9, no. 12, 2021, doi: 10.3390/pr9122134.
rho_p = 4750; %kg/m^3 density of a pellet 10% porosity, from Da Costa thesis
%A. Ranzani Da Costa, D. Wagner, and F. Patisson, “Modelling a new, low CO2 emissions, hydrogen steelmaking process,” J. Clean. Prod., vol. 46, pp. 27–35, 2013, doi: 10.1016/j.jclepro.2012.07.045.
rho_bed = rho_p*(1-eps_bed);

% Molar mass of species in kg/mol
MM_Fe2O3 = 159.69/1000;
MM_Fe3O4 = 251.53/1000;
MM_FeO = 71.844/1000;
MM_Fe = 55.845/1000;

% Molar mass of species in g/mol
MM_H2O = 18.015;
MM_H2 = 2.016;
MM_N2 = 28.013;
MM_CO = 28.01;
MM_CO2 = 44.01;

R = 8.314; %(m^3*Pa)/(K*mol)

load('initcond.mat')

T_ginit = interp1([1:1:length(T_ginit)],T_ginit, linspace(1,length(T_ginit),n_furnace));
%T_ginit = ones(size(T_ginit))*1073;
T_sinit = interp1([1:1:length(T_sinit)],T_sinit, linspace(1,length(T_sinit),n_furnace));

c_H2Oinit = interp1([1:1:length(c_H2Oinit)],c_H2Oinit, linspace(1,length(c_H2Oinit),n_furnace));
c_H2init = interp1([1:1:length(c_H2init)],c_H2init, linspace(1,length(c_H2init),n_furnace));
c_N2init = 0*ones(1, n_furnace);
c_COinit = 0*ones(1, n_furnace);
c_CO2init = 0*ones(1, n_furnace);


c_Feinit = interp1([1:1:length(c_Feinit)],c_Feinit, linspace(1,length(c_Feinit),n_furnace));
c_FeOinit = interp1([1:1:length(c_FeOinit)],c_FeOinit, linspace(1,length(c_FeOinit),n_furnace));
c_Fe3O4init = interp1([1:1:length(c_Fe3O4init)],c_Fe3O4init, linspace(1,length(c_Fe3O4init),n_furnace));
c_Fe2O3init = interp1([1:1:length(c_Fe2O3init)],c_Fe2O3init, linspace(1,length(c_Fe2O3init),n_furnace));

%c_Fe2O3init = rho_p*ones(1, n_furnace);

nr1init = interp1([1:1:length(nr1init)],nr1init, linspace(1,length(nr1init),n_furnace));
nr2init = interp1([1:1:length(nr2init)],nr2init, linspace(1,length(nr2init),n_furnace));
nr3init = interp1([1:1:length(nr3init)],nr3init, linspace(1,length(nr3init),n_furnace));

% nr1init = 0*ones(1, n_furnace);
% nr2init = 0*ones(1, n_furnace);
% nr3init = 0*ones(1, n_furnace);

DRIflow = 50.46; %kg/s
Reducerflow = 8.4309; %kg/s
%Reducerflow = 12.675; % kg/s

% Midrex newer paper
r_furnace = 2.75; % m
h_furnace = 9.68; % m

% da Costa reference
DRIflow = 52; %kg/s
Reducerflow = 8.49; %kg/s 

Reducerflow = 15;

% Midrex reference
r_furnace = 6.6/2; %m, radius of furnace
h_furnace = 6; % m, height of furnace

A_furnace = pi*r_furnace^2; %m^2, c.s. area of flow

A_furnace_pel = A_furnace*(1-eps_bed); %c.s. area for pellet flow - excludes gas lanes

%r_p = 6e-3; %m

r_p = 14e-3/2; %m, da Costa thesis

V_furnace = A_furnace*h_furnace;
V_pellet_bed = ((4/3)*pi*r_p^3)/(1-eps_bed);

zpts = 1:1:n_furnace+2;

b = log(h_furnace+1)/(n_furnace+1);
a = 1/exp(b);

dzfun = a*exp(b.*zpts)-1;
dz = flip(diff(dzfun));


dz = ones(1, n_furnace+1)* (h_furnace/(n_furnace+1));

%V_g = (4/3)*pi*r_p^3*n_pellets_dz*(eps_bed/(1-eps_bed));

n_pellets = V_furnace/V_pellet_bed;
% n_pellets_dz = n_pellets/(n_furnace+1);

n_pellets_dz = n_pellets.*dz/h_furnace;

% Volume of gas and solid in each spatial node - constant
V_g = (4/3)*pi*r_p^3*n_pellets_dz*(eps_bed/(1-eps_bed)); % m^3, volume of gas
V_s = (4/3)*pi*r_p^3*n_pellets_dz; % m^3, volume of solid
% A_g = V_g/dz;

a_b = 6*(1-eps_bed)/(r_p*2); %m^2/m^3, suface area for gas solid heat exchange, Wagner thesis
%a_sh = (2*pi*r_furnace*dz)/(A_furnace*dz); %m^2/m^3, surface area to lose heat to environment

A_wall = (2*pi*r_furnace*dz); %m^2 surface area of furnace
