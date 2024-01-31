syms n R T u_i u_im mu_i mu_im eps d rho_i rho_im Pi Pim A_f dz

%Per_im = -dz*( (150*(1-eps)^2*mu_im*u_im)/(d^2*eps^2)...
%    + ((1.75*rho_im*(1-eps)*u_im^2)/(d*eps) ) ) ;
% 
% 

% Per_i = dz*( (150*(1-eps)^2*mu_i*u_i)/(d^2*eps^3)...
%     + ((1.75*rho_i*(1-eps)*u_i^2)/(d*eps^3) ) ) ;

%delP = dz*( (150*(1-eps)^2*mu_i*u_i)/(d^2*eps^2)...
%     + ( (1.75*rho_i*(1-eps)*u_i^2)/(d*eps) ) ) ;

% eqn = u_i == ((n*R*T)/ ...
%     ( (Pim -  (Per_im-Per_i) )*A_f*eps)...
%     ) ;

u_i = (n*R*T)/(Pi*A_f*eps);

t1 =  (150*(1-eps)^2*mu_i*u_i)/(d^2*eps^2);

t2 = (1.75*rho_i*(1-eps)*u_i^2)/(d*eps)  ;


% eqn = u_i == ((n*R*T)/ ...
%     ( (-Per_i + 2*Pim) *A_f*eps)...
%     ) ;

eqn = (Pi - Pim)/dz == -t1 - t2;


Sa = solve(eqn, Pi, "MaxDegree", 3);
%Sa = solve(eqn, Pi);


simplify(Sa)
