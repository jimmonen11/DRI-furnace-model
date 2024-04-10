function [V1, V2, V3] = USCM_CO(X1, X2, X3, r0, T, c_CO, ct, x_og, kf)
% Unreacted shrinking core model using CO

%% Clean Up

% %Makes sure X values are positive

tol = 1e-15;
if X1 <= tol && X2 <= tol
    rxns = 1;
elseif X2 <= tol
    rxns = 2;
else
    rxns = 3;
end

%% Effective Diffusion Relations
%From:
%S. Yu, L. Shao, Z. Zou, and H. Saxén, “A numerical study on the performance of the h2 shaft furnace with dual-row top gas recycling,” Processes, vol. 9, no. 12, 2021, doi: 10.3390/pr9122134.
%which references:
%R. Takahashi, Y. Takahashi, J. Yagi, and Y. Omori, “Operation and Simulation of Pressurized Shaft Furnace for Direct Reduction,” Trans. Iron Steel Inst. Japan, vol. 26, no. 9, pp. 765–774, 1986.
D_CO_1 = exp(8.76 - 14.1e3/T)/(100^2); %m^2/s
D_CO_2 = exp(2.77 - 7.2e3/T)/(100^2); %m^2/s
D_CO_3 = exp(5.09 - 8.8e3/T)/(100^2);%m^2/s



%% Equilibrium Constants and Parameters

%T in Kelvin
R = 8.314; %J/mol*K

%Thermodynamic Analyses of Iron Oxides Redox Reactions 
K1 = exp((5128.6/T)+5.7);
K2 = exp((-3132.5/T)+3.661);
K3 = exp((2240.6/T)-2.667);

c_COeq1 = ((x_og)/(K1+1))*ct;
c_COeq2 = ((x_og)/(K2+1))*ct;
c_COeq3 = ((x_og)/(K3+1))*ct;

%% Reaction rate Coefficient

% Takahashi, Yagi, Operation for Direct and Simulation Reduction*
k1 = exp(3.16 - 50.2/(R*1e-3*T))/100;
k2 = exp(2.09 - 40/(R*1e-3*T))/100;
k3 = exp(5.42 - 61.4/(R*1e-3*T))/100;


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


if X2 >= 0.99
    V3 = 0;
end
