import numpy as np
import pybindingcurve as pbc
my_system = pbc.BindingCurve("1:1")

system_parameters = {"p": np.linspace(0, 900), "l": 100, "kdpl": 10}

my_system.add_curve(system_parameters, name= 'Curve 1')

my_system.show_plot()

print(my_system.query({'p':5, 'l':100, 'kdpl':10}))

