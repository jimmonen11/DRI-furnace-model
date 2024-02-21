function [V1, V2, V3, x_COeq1, x_COeq2, x_COeq3] = USCM_CO(X1, X2, X3, r0, T, c_CO, ct, x_og, rxns, kf)

%% Clean Up

% %Makes sure X values are positive

tol = 1e-10;
if X1 <= tol && X2 <= tol
    rxns = 1;
elseif X2 <= tol
    rxns = 2;
else
    rxns = 3;
end

%rxns = 3;

%% Effective Diffusion Relations
%From:
%S. Yu, L. Shao, Z. Zou, and H. Saxén, “A numerical study on the performance of the h2 shaft furnace with dual-row top gas recycling,” Processes, vol. 9, no. 12, 2021, doi: 10.3390/pr9122134.
%which references:
%R. Takahashi, Y. Takahashi, J. Yagi, and Y. Omori, “Operation and Simulation of Pressurized Shaft Furnace for Direct Reduction,” Trans. Iron Steel Inst. Japan, vol. 26, no. 9, pp. 765–774, 1986.
D_CO_1 = exp(8.76 - 14.1e3/T)/(100^2); %m^2/s
D_CO_2 = exp(2.77 - 7.2e3/T)/(100^2); %m^2/s
D_CO_3 = exp(5.09 - 8.8e3/T)/(100^2);%m^2/s

%kf = 0.01; %m/s
%kf = 0.47; %m/s

% D_H2_1 = 4.767751e-5;
% D_H2_2 = 4.767751e-5;
% D_H2_3 = 4.767751e-5;


%% Equilibrium Constants and Parameters
%From "The Effect of Gas and Solids Maldistribution on the Performance of
%Moving-bed Reactors: The Reduction of Iron Oxide Pelletswith Hydrogen"
%T in Kelvin
R = 8.314; %J/mol*K

% G values from "The production and application of hydrogen in Steel industry"
% Liu et al. 
% G1 = -52870 - 47.33*(T);
% G2 = 36250 - 51.98*(T);
% G3 = -17310 + 17.32*(T);
% 
% K1 = exp(-G1/(R*T));
% K2 = exp(-G2/(R*T));
% K3 = exp(-G3/(R*T));


% K1 = exp((362/T)+10.32);
% K2 = exp((-8580/T)+8.98);
% K3 = exp((-2070/T)+1.30);

%Thermodynamic Analyses of Iron Oxides Redox Reactions 
K1 = exp((5128.6./T)+5.7);
K2 = exp((-3132.5./T)+3.661);
K3 = exp((2240.6./T)-2.667);

K1 = exp(6.0259 - 0.000525*T + 5635.3/T);
K2 = exp(12.2909 - 0.006183*T - 5388.3/T);
K3 = exp(-1.5304 + 0.000729*T +1610.03/T);

c_COeq1 = ((x_og)/(K1+1))*ct;
c_COeq2 = ((x_og)/(K2+1))*ct;
c_COeq3 = ((x_og)/(K3+1))*ct;

x_COeq1 = c_COeq1/ct;
x_COeq2 = c_COeq2/ct;
x_COeq3 = c_COeq3/ct;

%% Reaction rate Coefficient

%wagner thesis m/s
% k1 = 6.3587e-2*exp(-43276/(R*T));
% k2 = 1.1290e-3*exp(-21570/(R*T));
% k3 = 3.8337e-2*exp(-50129/(R*T));


%da costa thesis m/s
% k1 = 7.79e-4*exp(-27000/(R*T));
% k2 = 1.11e-2*exp(-55000/(R*T));
% k3 = 16*exp(-136000/(R*T));

%%From "The Effect of Gas and Solids Maldistribution on the Performance of
%Moving-bed Reactors: The Reduction of Iron Oxide Pelletswith Hydrogen"
% k1 = 1.44e5*exp(-6650/T)/3600;
% k2 = 2.88e5*exp(-8000/T)/3600;
% k3 = 2.45e7*exp(-1400/T)/3600;

% Takahashi, Yagi, Operation for Direct and Simulation Reduction*

% k1 = exp(3.16 - 50.2/(R*1e-3*T))/100;
% k2 = exp(2.09 - 40/(R*1e-3*T))/100;
% k3 = exp(5.42 - 61.4/(R*1e-3*T))/100;

%Hara
% k1 = 41.7*exp(-66512/(R*T))/5; %Hara
% k2 = 22.2*exp(-74826/(R*T))/5;
% k3 = 4083.3*exp(-116396/(R*T))/5;

k1 = 29.2*exp(-66974/(R*T))/5; %Takenaka
k2 = 15.6*exp(-75345/(R*T))/5;
k3 = 2858.3*exp(-66974/(R*T))/5;


%French latest thesis
% k1 = 13*exp(-113859/(R*T));
% k2 = 0.053*exp(-73674/(R*T));
% k3 = 0.0027*exp(-69488/(R*T));
% 
% k1 = -26177.6 -92.125*T + 0.096808*T^2;
% k2 = 7744.5 +83.5404*T - 0.09669*T^2;
% k3 = -4527.4 -34.5157*T + 0.023612*T^2;
% 
% k1 = 2700*exp(-113859/(R*T));
% k2 = 25*exp(-73674/(R*T));
% k3 = 17*exp(-69488/(R*T));
% % 
% 
% %French latest thesis 
% k1 = 13*exp(-113859/(R*T));
% k2 = 0.053*exp(-73674/(R*T));
% k3 = 0.0027*exp(-69488/(R*T));

%%
A1 = 1/(X1^2*(k1*(1+1/K1)));
A2 = 1/(X2^2*(k2*(1+1/K2)));
A3 = 1/(X3^2*(k3*(1+1/K3)));

B1 = (X2-X1)*r0/(D_CO_1*X1*X2);
B2 = (X3-X2)*r0/(D_CO_2*X2*X3);
B3 = (1-X3)*r0/(D_CO_3*X3);

F = 1/kf;

W1 = (A1+B1)*(A3*(A2+B2+B3+F) + (A2+B2)*(B3+F)) + A2*(A3*(B2+B3+F) + B2*(B3+F));
W2 = (A2+B2)*(A3+B3+F) + A3*(B3+F);
W3 = A3 + B3 + F;

%%

if rxns == 3

    V1 = 4*pi*r0^2*( ((A3*(A2+B2+B3+F)) + (A2+B2)*(B3+F))*(c_CO-c_COeq1)...
        - (B3+F)*A2*(c_CO-c_COeq3) - (A3*(B2+B3+F)+B2*(B3+F))*(c_CO-c_COeq2)) / W1;

    V2 = 4*pi*r0^2*( ((A1+B1+B2)*(A3+B3+F) + A3*(B3+F))*(c_CO-c_COeq2)...
        - (B2*(A3+B3+F)+(A3*(B3+F)))*(c_CO-c_COeq1) - ((A1+B1)*(B3+F))*(c_CO-c_COeq3)) / W1;

    V3 = 4*pi*r0^2*( ((A1+B1)*(A2+B2+B3+F) + A2*(B2+B3+F))*(c_CO-c_COeq3)...
        - (A2*(B3+F))*(c_CO-c_COeq1) - ((A1+B1)*(B3+F))*(c_CO-c_COeq2)) / W1;

elseif rxns == 2

    V1 = 0;

    V2 = 4*pi*r0^2*( (A3+B3+F)*(c_CO-c_COeq2)...
        - (B3+F)*(c_CO-c_COeq1) )/ W2 ;

    V3 = 4*pi*r0^2*( ((A1+B1)*(A2+B2+B3+F) + A2*(B2+B3+F))*(c_CO-c_COeq3)...
        - (A2*(B3+F))*(c_CO-c_COeq1) - ((A1+B1)*(B3+F))*(c_CO-c_COeq2)) / W2;
else
    
    V1 = 0;
    
    V2 = 0;
    
    V3 = 4*pi*r0^2*(c_CO-c_COeq3)/ W3; 
    
    
end


if X2 > 0.99
    V3 = 0;
end
