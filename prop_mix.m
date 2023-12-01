function [rho, mu, k, Cp] = prop_mix(rhos, mus, ks, Cps, x_fracs)

%Function that gives weighted average of properties

%using mole weighted average because gas properties
rho = dot(rhos, x_fracs);
mu = dot(mus,x_fracs);
k = dot(ks,x_fracs);
Cp = dot(Cps, x_fracs);