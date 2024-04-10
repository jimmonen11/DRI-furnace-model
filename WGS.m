function v = WGS(X, T, P, Vp, x_H2, x_H2O, x_CO, x_CO2)

P = P/101325; %atm, convert from Pa

P_H2 = x_H2 * P;
P_H2O = x_H2O * P;
P_CO = x_CO * P;
P_CO2 = x_CO2 * P;

R = 8.314;

Keq  = exp(4400/T - 4.063); %from Hamadeh who cites J. Xu and G. F. Froment, “Methane steam reforming, methanation and water‐gas shift: I. Intrinsic kinetics,” AIChE Journal, vol. 35, no. 1, pp. 88–96, 1989, doi: 10.1002/aic.690350109.

% Takahashi, R., Takahashi, Y., Yagi, J., & Omori, Y. (1986). Operation and Simulation of Pressurized Shaft Furnace for Direct Reduction. Transactions of the Iron and Steel Institute of Japan, 26(9), 765–774.
% Changed to 0.25 to better validate
if X < 0.25
    kf = 1.83e-5*exp(7.84e-3/(R*1e-3*T))*100^3;
else

    kf = 9.33e1*exp(-7.32/(R*1e-3*T))*100^3;
    
end


v = (Vp*kf*(P_CO*P_H2O - P_CO2*P_H2/Keq));


