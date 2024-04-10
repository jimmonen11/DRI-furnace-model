function [V1, V2, V3] = USCM_H2(X1, X2, X3, r0, T, c_H2, ct, x_og, kf)
% Unreacted shrinking core model using H2

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
D_H2_1 = exp(3.43 - 4.2e3/T)/(100^2); %m^2/s
D_H2_2 = exp(5.64 - 6.8e3/T)/(100^2); %m^2/s
D_H2_3 = exp(4.77 - 5.9e3/T)/(100^2);%m^2/s


%% Equilibrium Constants and Parameters

%T in Kelvin
R = 8.314; %J/mol*K

%Thermodynamic Analyses of Iron Oxides Redox Reactions 
K1 = exp((1433.37/T)+9.083);
K2 = exp((-7393.9/T)+7.563);
K3 = exp((-2023.8/T)+1.239);

c_H2eq1 = ((x_og)/(K1+1))*ct;
c_H2eq2 = ((x_og)/(K2+1))*ct;
c_H2eq3 = ((x_og)/(K3+1))*ct;


%% Reaction rate Coefficient

% UyS. et al. A numerical study on the performance of the h2 shaft furnace with dual-row top gas recycling
% Takahashi, Yagi, Operation for Direct and Simulation Reduction*
 
k1 = exp(4.49 - 33.4/(R*1e-3*T))/100;
k2 = exp(6.7 - 58.2/(R*1e-3*T))/100;
k3 = exp(6.97 - 57.1/(R*1e-3*T))/100;

%%
A1 = 1/(X1^2*(k1*(1+1/K1)));
A2 = 1/(X2^2*(k2*(1+1/K2)));
A3 = 1/(X3^2*(k3*(1+1/K3)));

B1 = (X2-X1)*r0/(D_H2_1*X1*X2);
B2 = (X3-X2)*r0/(D_H2_2*X2*X3);
B3 = (1-X3)*r0/(D_H2_3*X3);

F = 1/kf;

W1 = (A1+B1)*(A3*(A2+B2+B3+F) + (A2+B2)*(B3+F)) + A2*(A3*(B2+B3+F) + B2*(B3+F));
W2 = (A2+B2)*(A3+B3+F) + A3*(B3+F);
W3 = A3 + B3 + F;

%%

if rxns == 3

    V1 = 4*pi*r0^2*( ((A3*(A2+B2+B3+F)) + (A2+B2)*(B3+F))*(c_H2-c_H2eq1)...
        - (B3+F)*A2*(c_H2-c_H2eq3) - (A3*(B2+B3+F)+B2*(B3+F))*(c_H2-c_H2eq2)) / W1;

    V2 = 4*pi*r0^2*( ((A1+B1+B2)*(A3+B3+F) + A3*(B3+F))*(c_H2-c_H2eq2)...
        - (B2*(A3+B3+F)+(A3*(B3+F)))*(c_H2-c_H2eq1) - ((A1+B1)*(B3+F))*(c_H2-c_H2eq3)) / W1;

    V3 = 4*pi*r0^2*( ((A1+B1)*(A2+B2+B3+F) + A2*(B2+B3+F))*(c_H2-c_H2eq3)...
        - (A2*(B3+F))*(c_H2-c_H2eq1) - ((A1+B1)*(B3+F))*(c_H2-c_H2eq2)) / W1;

elseif rxns == 2

    V1 = 0;

    V2 = 4*pi*r0^2*( (A3+B3+F)*(c_H2-c_H2eq2)...
        - (B3+F)*(c_H2-c_H2eq1) )/ W2 ;

    V3 = 4*pi*r0^2*( ((A1+B1)*(A2+B2+B3+F) + A2*(B2+B3+F))*(c_H2-c_H2eq3)...
        - (A2*(B3+F))*(c_H2-c_H2eq1) - ((A1+B1)*(B3+F))*(c_H2-c_H2eq2)) / W2;
else
    
    V1 = 0;
    
    V2 = 0;
    
    V3 = 4*pi*r0^2*(c_H2-c_H2eq3)/ W3; 
    
    
end


if X2 >= 0.99
    V3 = 0;
end

