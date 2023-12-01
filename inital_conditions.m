

c_Fe2O3init = out.Fe2O3_conc.data(end,:);
c_Fe3O4init = out.Fe3O4_conc.data(end,:);
c_FeOinit = out.FeO_conc.data(end,:);
c_Feinit = out.Fe_conc.data(end,:);

c_H2init = out.H2_conc.data(end,:);
c_H2Oinit = out.H2O_conc.data(end,:);

T_ginit = out.T_g.data(end,:);
T_sinit = out.T_s.data(end,:);


nr1init = out.nr1.data(end,:);
nr2init = out.nr2.data(end,:);
nr3init = out.nr3.data(end,:);


save("initcond.mat","c_Fe2O3init","c_Fe3O4init","c_FeOinit","c_Feinit","c_H2init", "c_H2Oinit", ...
    "T_ginit","T_sinit", "nr1init", "nr2init","nr3init")