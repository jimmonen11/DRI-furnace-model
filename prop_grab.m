function [rho, mu, k, C] = prop_grab(T, P, species)

%Function that returns properties from CoolProp

coder.extrinsic('py.CoolProp.CoolProp.PropsSI')

rho = py.CoolProp.CoolProp.PropsSI('D','T',T,'P',P,species); %density, kg/m^3
mu  = py.CoolProp.CoolProp.PropsSI('V','T',T,'P',P,species); %viscocity, Pa*s
k   = py.CoolProp.CoolProp.PropsSI('L','T',T,'P',P,species); %conductive heat trans coef (W/mK)
C   = py.CoolProp.CoolProp.PropsSI('C','T',T,'P',P,species); %specific heat capacity (J/kgK)
    
end
