function v = SMR(X, T, P, Vp, x_H2, x_H2O, x_CO, x_CH4)


P = P/101325; %atm, convert from Pa

P_H2 = x_H2 * P;
P_H2O = x_H2O * P;
P_CO = x_CO * P;
P_CH4 = x_CH4 * P;

R = 8.314;


 %from Hamadeh who cites J. Xu and G. F. Froment, “Methane steam reforming, methanation and water‐gas shift: I. Intrinsic kinetics,” AIChE Journal, vol. 35, no. 1, pp. 88–96, 1989, doi: 10.1002/aic.690350109.
Keq = (10266.76*10^6*exp (-26830/T + 30.11))/101325^2; %atm^2
 %G = -236.58*T + 187438;
%Keq = exp(-G/(R*T));

% changed to 0.25 to better validate
if X < 0.25
    kf = 1e-20; %mol/s-m^3-atm^4

else

    % Three-dimensional simulation of chemically reacting gas f lows in the porous support structure of an integrated-planar solid oxide fuel cell
    % Modified to give better validation
    kf = (2395*exp(-180000/(R*T)) ) * 101325^2;

end

v = Vp* kf * (P_CH4*P_H2O - P_CO*P_H2^3/Keq);


