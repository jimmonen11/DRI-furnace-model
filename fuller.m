function D = fuller(Speciesi, Speciesj, T, P)
%species i values
%M, molecular weight of species i (g/mol)
%V, diffusion volume of species i (unitless)
Mi = Speciesi(1);
Vi = Speciesi(2);

%species j values
%M, molecular weight of species j (g/mol)
%V, diffusion volume of species j (unitless)
Mj = Speciesj(1);
Vj = Speciesj(2);

Cf = 1.013e-2; %constant

D = Cf*T^1.75*(sqrt((Mi+Mj)/(Mi*Mj)))/(P*(Vi^(1/3)+Vj^(1/3))^2);
