function [V1, V2, V3] = USCM(X1, X2, X3, r0, T, c_H2, ct, rxns, kf)

%% Clean Up

% %Makes sure X values are positive
% if X1 < 1e-5
%     X1 = 1e-5;
% end
% 
% if X2 < 1e-5
%     X2 = 1e-5;
% end
% 
% if X3 < 1e-5
%     X3 = 1e-5;
% end

%% Effective Diffusion Relations
%From:
%S. Yu, L. Shao, Z. Zou, and H. Saxén, “A numerical study on the performance of the h2 shaft furnace with dual-row top gas recycling,” Processes, vol. 9, no. 12, 2021, doi: 10.3390/pr9122134.
%which references:
%R. Takahashi, Y. Takahashi, J. Yagi, and Y. Omori, “Operation and Simulation of Pressurized Shaft Furnace for Direct Reduction,” Trans. Iron Steel Inst. Japan, vol. 26, no. 9, pp. 765–774, 1986.
D_H2_1 = exp(3.43 - 4.2e3/T)/(100^2); %m^2/s
D_H2_2 = exp(5.65 - 6.8e3/T)/(100^2); %m^2/s
D_H2_3 = exp(4.77 - 5.9e3/T)/(100^2);%m^2/s

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

G1 = -155547 - 74.7*(T);
G2 = 71490 - 73.12*(T);
G3 = 23430 - 16.16*(T);

K1 = exp(-G1/(R*T));
K2 = exp(-G2/(R*T));
K3 = exp(-G3/(R*T));


% K1 = exp((362/T)+10.32);
% K2 = exp((-8580/T)+8.98);
% K3 = exp((-2070/T)+1.30);


c_H2eq1 = (1/(K1+1))*ct;
c_H2eq2 = (1/(K2+1))*ct;
c_H2eq3 = (1/(K3+1))*ct;


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
k1 = 1.44e5*exp(-6650/T)/3600;
k2 = 2.88e5*exp(-8000/T)/3600;
k3 = 2.45e7*exp(-1400/T)/3600;

% %UyS. et al. A numerical study on the performance of the h2 shaft furnace with dual-row top gas recycling 
% k1 = exp(4.49 - 33.4/(R*T))/100;
% k2 = exp(6.7 - 58.2/(R*T))/100;
% k3 = exp(6.97 - 57.1/(R*T))/100;


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

