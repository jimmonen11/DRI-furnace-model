S_CO = 197.66;
S_CO2 = 213.79;

S_Fe2O3 = 87.2;
S_Fe3O4 = 146.4;
S_FeO = 60.75;
S_Fe = 27.3;

S1 = (2*S_Fe3O4 + S_CO2) - (3*S_Fe2O3 + S_CO);
S2 = (3*S_FeO + S_CO2) - (S_Fe3O4 + S_CO);
S3 = (S_Fe + S_CO2) - (S_FeO + S_CO);


% H1 = -12636*4.184;
% H2 = 8664*4.184;
% H3 = -4136*4.184;

H1 = -52.87e3;
H2 = 36.25e3;
H3 = -17.31e3;

T = 570:1:1400;

G1 = H1 - S1*(T);
G2 = H2 - S2*(T);
G3 = H3 - S3*(T);

K1 = exp(-G1./(R*T));
K2 = exp(-G2./(R*T));
K3 = exp(-G3./(R*T));


% K1 = exp((362/T)+10.32);
% K2 = exp((-8580/T)+8.98);
% K3 = exp((-2070/T)+1.30);


x_COeq1 = (1./(K1+1));
x_COeq2 = (1./(K2+1));
x_COeq3 = (1./(K3+1));

close all
plot(T-273, x_COeq1)
hold on
plot(T-273, x_COeq2)
plot(T-273, x_COeq3)

