function [delH1, delH2, delH3, delH, delG1, delG2, delG3] = heat_rxn_CO(T)

H_CO2 = py.CoolProp.CoolProp.PropsSI('Hmolar','T',T,'P',101325,"CO2")...
    - py.CoolProp.CoolProp.PropsSI('Hmolar','T',298,'P',101325,"CO2")+ (-393509);

H_CO = py.CoolProp.CoolProp.PropsSI('Hmolar','T',T,'P',101325,"CO")...
    - py.CoolProp.CoolProp.PropsSI('Hmolar','T',298,'P',101325,"CO") + (-110525);

S_CO2 = py.CoolProp.CoolProp.PropsSI('Smolar','T',T,'P',101325,"CO2")...
    - py.CoolProp.CoolProp.PropsSI('Smolar','T',298,'P',101325,"CO2")+ (213.79);

S_CO = py.CoolProp.CoolProp.PropsSI('Smolar','T',T,'P',101325,"CO")...
    - py.CoolProp.CoolProp.PropsSI('Smolar','T',298,'P',101325,"CO") + (197.66);


[~, H_irons, S_irons] = iron_props(T);

H_Fe2O3  = H_irons(1);
H_Fe3O4  = H_irons(2);
H_FeO  = H_irons(3);
H_Fe  = H_irons(4);

S_Fe2O3  = S_irons(1);
S_Fe3O4  = S_irons(2);
S_FeO  = S_irons(3);
S_Fe  = S_irons(4);

% Overall
% 3Fe2O3 + 9H2 ====> 6Fe + 9H2O 
%  Fe2O3 + 3H2 ====> 2Fe + 3H2O 
delH = ((2*H_Fe + 3*H_CO2) - (H_Fe2O3 + 3*H_CO))/3;
delS = ((2*S_Fe + 3*S_CO2) - (S_Fe2O3 + 3*S_CO))/3;


% First
% 3Fe2O3 + H2 ====> 2Fe3O4 + H2O
delH1 = ((2*H_Fe3O4 + H_CO2) - (3*H_Fe2O3 + H_CO));
delS1 = ((2*S_Fe3O4 + S_CO2) - (3*S_Fe2O3 + S_CO));

delG1 = delH1 - T*delS1;

% Second
% 2Fe3O4 + 2H2 ====> 6FeO + 2H2O
%  Fe3O4 +  H2 ====> 3FeO + H2O

delH2 = ((3*H_FeO + H_CO2) - (H_Fe3O4 + H_CO));
delS2 = ((3*S_FeO + S_CO2) - (S_Fe3O4 + S_CO));

delG2 = delH2 - T*delS2;


% Third
% 6FeO + 6H2 ====> 6Fe + 6H2O
delH3 = ((H_Fe + H_CO2) - (H_FeO + H_CO));
delS3 = ((S_Fe + S_CO2) - (S_FeO + S_CO));

delG3 = delH3 - T*delS3;

