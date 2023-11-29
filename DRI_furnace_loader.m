n_furnace = 10;

eps_bed = 0.4; %S. Yu, L. Shao, Z. Zou, and H. Saxén, “A numerical study on the performance of the h2 shaft furnace with dual-row top gas recycling,” Processes, vol. 9, no. 12, 2021, doi: 10.3390/pr9122134.
rho_p = 4750; %kg/m^3 %A. Ranzani Da Costa, D. Wagner, and F. Patisson, “Modelling a new, low CO2 emissions, hydrogen steelmaking process,” J. Clean. Prod., vol. 46, pp. 27–35, 2013, doi: 10.1016/j.jclepro.2012.07.045.
rho_bed = rho_p*(1-eps_bed);

c_H2Oinit = 1*ones(1, n_furnace);
c_H2init = 12*ones(1, n_furnace);
%c_H2Oinit = 2*ones(1, n_furnace);
%c_H2Oinit = 2*ones(1, n_furnace);

c_Fe2O3init = rho_bed*1.25*ones(1, n_furnace);
c_Fe3O4init = rho_bed/4*ones(1, n_furnace);
c_FeOinit = rho_bed/4*ones(1, n_furnace);
c_Feinit = rho_bed/4*ones(1, n_furnace);

T_sinit = (500+273)*ones(1, n_furnace);
T_ginit = (900+273)*ones(1, n_furnace);

nrinit = zeros(1, n_furnace);

c_H2in = 12;
c_H2Oin = 1;

%DRIflow = 50.33; %kg/s
DRIflow = 52; %kg/s
Reducerflow = 8.4309; %kg/s
%Reducerflow = 10.685; % kg/s


reducerflowmole = Reducerflow/4 * 1000;

r_furnace = 6.6/2; %m, radius of furnace
h_furnace = 6; %m, height of furnace
A_furnace = pi*r_furnace^2; %m^2, c.s. area of flow

V_furnace = A_furnace*h_furnace;

%rho_bed = 4750; %kg/m^3, density of bed
%rho_gas = 0.0374; %kg/m^3, density of gas


%u_s = DRIflow/rho_bed/A_furnace;
%u_g = Reducerflow/rho_gas/A_furnace;

dz = h_furnace/n_furnace;


r_p = 6e-3; %m


a_b = 6*(1-eps_bed)/(r_p*2); %m^2/m^3, suface area for gas solid heat exchange


pellet_flow = (3*DRIflow)/(4*pi*r_p^3*rho_bed); %s^-1
tau_furnace = h_furnace*A_furnace/(DRIflow/rho_bed)
n_pellets = pellet_flow*tau_furnace/(n_furnace);


rho_p = 5275; %kg/m^3, Fe2O3
V_pellet = (4/3)*pi*r_p^3;
m_pellet = rho_p*V_pellet;

M_Fe2O3 = 159.69; %g/mol

molesO2pellet = 0.63*m_pellet/M_Fe2O3*1000;
kappa = molesO2pellet%m_pellet;

T = 800+273;
P = 101325*1.16; %thesis de Costa
R = 8.314; %(m^3*Pa)/(K*mol)
ct = P/(R*T)

MM_H2 = 2.016;
MM_H2O = 18.015;
MM_N2 = 28.02;

x_H2in = 0.98;
x_H2Oin = 1 - x_H2in;

%c_H2in = ct*x_H2in;
%c_H2Oin = ct*x_H2Oin;

n_gas_in = 3634; %mol/s
flow = (n_gas_in*MM_H2*x_H2in + n_gas_in*MM_H2O*x_H2Oin)/1000;

% c_H2Oinit = n_gas_in*x_H2Oin*ones(1, n_furnace);
% c_H2init = n_gas_in*x_H2Oin*ones(1, n_furnace);

%Actually flows
% c_H2in = n_gas_in*x_H2in;
% c_H2Oin = n_gas_in*x_H2Oin;
% 
% c_Fe2O3init = 52*ones(1, n_furnace);
% c_Fe3O4init = 1e-20*ones(1, n_furnace);
% c_FeOinit = 1e-20*ones(1, n_furnace);
% c_Feinit = 1e-20*ones(1, n_furnace);
