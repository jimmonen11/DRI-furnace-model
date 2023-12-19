function [rho, mu, k, Cp] = prop_mix(rhos, mus, ks, Cps, w_fracs)

%Function that gives weighted average of properties

%using mass weighted average
rho = dot(rhos, w_fracs);
mu = dot(mus,w_fracs);
k = dot(ks,w_fracs);
Cp = dot(Cps, w_fracs);