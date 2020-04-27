# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 16:54:11 2017

@author: sm905
"""

# DEMO SCRIPT TO SHOW BASIC USE OF THE MULTILITHMELTING CODE
# More detailed documentation is procided within the module, and should be
# accessible by the spyder object inspector.


import multilith as m

# Define lithologies. Calling the functions will use the default parameters.
# You can change the constants used in the functions, or the thermodynamic constants
# by setting them when you create the lithology object, or later by calling the
# constant directly. Once the lithology object is defined you can call functions
# e.g. TSolidus to directly utilise its properties.
lz = m.LithologyKatz()
g2 = m.LithologySimple()
kg1 = m.LithologyShorttle()
hz = m.LithologyNonMelting()

# Define a mantle object. A mantle object is made up of one or more lithologies,
# in a specified proportion, with specified names.
mantle = m.mantle([kg1,lz,hz],[1,10,2],['px','lz','hz'])

# You can plot many of the important curves that control the behaviour of the model
# by calling this function. It will also return a figure object, should you want to
# customise the plot further.
mantle.PlotBoundaries()

# Create a melting column object. Currently the code can only deal with 1D adiabatic
# upwelling. The number passed is the Tp in degC. You can further specify how the
# calculation should be run, if you wish. See documentation. To directly extract
# properties from this object use r['lithology name'] (for F), r.P and r.Temperature etc.
# r.F_total will give the total melt fraction.
r = mantle.AdiabaticMelt_1D(1400)

# Create a default plot of the thermal gradient in the melting region, and the
# melt fraction of each lithology.
r.plot()

# Integrate the melt fraction assuming passive upwelling over a triangular corner-flow
# melting region. Adds various properties to the melting column object.
r.integrate_tri()

# Calls one of the new properties of the melting column object, the crustal thickness,
# and prints it.
print('The crustal thickness is: ' + str(r.tc) + ' km')