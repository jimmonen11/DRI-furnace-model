% Need to change initcond_H2.mat or init_cond.mat depending on what initial
% conditions you want to save!

c_Fe2O3init = out.Fe2O3_conc.data(end,:);
c_Fe3O4init = out.Fe3O4_conc.data(end,:);
c_FeOinit = out.FeO_conc.data(end,:);
c_Feinit = out.Fe_conc.data(end,:);
c_Cinit = out.C_conc.data(end,:);

c_H2init = out.H2_conc.data(end,:);
c_H2Oinit = out.H2O_conc.data(end,:);
c_N2init = out.N2_conc.data(end,:);
c_COinit = out.CO_conc.data(end,:);
c_CO2init = out.CO2_conc.data(end,:);
c_CH4init = out.CH4_conc.data(end,:);

T_ginit = out.T_g.data(end,2:end);
T_sinit = out.T_s.data(end,1:end-1);

nr1init = out.nr1.data(end,:);
nr2init = out.nr2.data(end,:);
nr3init = out.nr3.data(end,:);

ndotinit = out.ndot.data(end,2:end);


save("initcond_NG.mat","c_Fe2O3init","c_Fe3O4init","c_FeOinit","c_Feinit", "c_Cinit", "c_H2init", "c_H2Oinit", ...
    "c_N2init", "c_COinit", "c_CO2init", "c_CH4init", "T_ginit","T_sinit", "nr1init", "nr2init","nr3init", "ndotinit")