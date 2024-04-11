# DRI-furnace-model
DRI furnace model

This is a dynamic model of a direct reduced iron (DRI) furnace.

With all the files in this repository in the same folder use the file 'DRI_project_loader.m' to load supporting parameters. There are
options to run steady state simulations and three default dynamic simulations to run. Full customization of dynamic simulations can be
performed by changing values in the loader file. Once the values are loaded in the workspace then can run the simulations through the
'DRI_furnace.slx' Simulink file. Varibles of interest are saved in the variable 'out' once the simulation is finished. Supporting 'plotter'
files are provided to visualize the results.


If you want to use the 'CoolPropGrabber.mlx' file to see where gas properties come from you need to have a Python installed and Matlab needs to know how to access it.
This is NOT required to run the simulation. See documentation below:
https://www.mathworks.com/help/matlab/matlab_external/install-supported-python-implementation.html#buialof-39

Once you have Matlab found by Python you need to install CoolProp on your given version of Python. See documentation below:
http://www.coolprop.org/coolprop/wrappers/MATLAB/index.html#matlab-wrapper

Check that CoolProp is working with this Matlab command: py.CoolProp.CoolProp.PropsSI('T','P',101325,'Q',0,'Water')

