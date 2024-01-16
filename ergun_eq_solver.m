syms n R T u_i u_im mu_i mu_im eps d rho_i rho_im Pim Per_i Per_im A_f dz

Per_im = -dz*( (150*(1-eps)^2*mu_im*u_im)/(d^2*eps^2)...
    + ((1.75*rho_im*(1-eps)*u_im^2)/(d*eps) ) ) ;


Per_i = -dz*( (150*(1-eps)^2*mu_i*u_i)/(d^2*eps^2)...
    + ((1.75*rho_i*(1-eps)*u_i^2)/(d*eps) ) ) ;

eqn = u_i - ((n*R*T)/ ...
    ((Pim - (Per_i - Per_im))*A_f*eps)...
    ) == 0 ;

Sa = solve(eqn,u_i, "MaxDegree", 3);

simplify(Sa)