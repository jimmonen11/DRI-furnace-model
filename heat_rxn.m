function [delH1, delH2, delH3, delH] = heat_rxn(T)

H_H2O = py.CoolProp.CoolProp.PropsSI('Hmolar','T',T,'P',101325,"water")...
    - py.CoolProp.CoolProp.PropsSI('Hmolar','T',298,'P',101325,"water")+ (-285830);

H_H2 = py.CoolProp.CoolProp.PropsSI('Hmolar','T',T,'P',101325,"hydrogen")...
    - py.CoolProp.CoolProp.PropsSI('Hmolar','T',298,'P',101325,"hydrogen");

[~, H_irons] = iron_props(T);

H_Fe2O3  = H_irons(1);
H_Fe3O4  = H_irons(2);
H_FeO  = H_irons(3);
H_Fe  = H_irons(4);

% Overall
% 3Fe2O3 + 9H2 ====> 6Fe + 9H2O 
%  Fe2O3 + 3H2 ====> 2Fe + 3H2O 
delH = ((2*H_Fe + 3*H_H2O) - (H_Fe2O3 + 3*H_H2))/3;

% First
% 3Fe2O3 + H2 ====> 2Fe3O4 + H2O
delH1 = ((2*H_Fe3O4 + H_H2O) - (3*H_Fe2O3 + H_H2));


% Second
% 2Fe3O4 + 2H2 ====> 6FeO + 2H2O
%  Fe3O4 +  H2 ====> 3FeO + H2O

delH2 = ((3*H_FeO + H_H2O) - (H_Fe3O4 + H_H2));


% Third
% 6FeO + 6H2 ====> 6Fe + 6H2O
delH3 = ((H_Fe + H_H2O) - (H_FeO + H_H2));

