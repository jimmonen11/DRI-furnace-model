# DRI-furnace-model
DRI furnace model

This is a model of a direct reduced iron (DRI) furnace.

To use you need to have a Python installed and Matlab needs to know how to access it. See documentation below:
https://www.mathworks.com/help/matlab/matlab_external/install-supported-python-implementation.html#buialof-39

Once you have Matlab found by Python you need to install CoolProp on your given version of Python. See documentation below:
http://www.coolprop.org/coolprop/wrappers/MATLAB/index.html#matlab-wrapper

Check that CoolProp is working with this Matlab command: py.CoolProp.CoolProp.PropsSI('T','P',101325,'Q',0,'Water')

