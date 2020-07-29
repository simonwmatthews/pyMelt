# pyMelt
A python package for calculating the melting behaviour of multi-lithologic mantle. The module implements the equations developed by Phipps Morgan (2001) to calculate the melting behaviour of mantle composed of any chosen lithology. The lithologies included in the module are:
- KLB-1 lherzolite (Matthews et al., in review)
- KG1 silica-undersaturated pyroxenite (Matthews et al., in review) 
- Silica-saturated pyroxenite (Matthews et al., in review) 
- KG2 silica-undersaturated pyroxenite (Shorttle et al. 2014)
- Lherzolite (Katz et al., 2003)
- G2 pyroxenite (Pertermann & Hirschmann, 2002)
- Harzburgite (non-melting)

Currently supported calculations:
- Adiabatic decompression melting
- Isobaric melting

Parameters that can be calculated:
- The geotherm for decompressing mantle
- Melt fractions for each lithology
- Crustal thickness for passive-upwelling at a mid-ocean ridge
- Crystallisation temperatures (following the method in Matthews et al., 2016)

## Getting started
First the pyMelt module must be imported. If it is in the same directory as you are working in:
```
import pyMelt as m
```
Each lithology is treated as a python object, and an instance can be assigned to a variable:
```
lz = m.LithologyKLB1()
px = m.LithologyKG1()
```  
Any constant used in the model can be changed in this call. See the documentation for each Lithology class for more information. To find more lithologies, look through the pyMelt module with a text editor, or if you're in an ipython environment press tab after typing `m.`.

Next a mantle object must be constructed, which is built up of lithologies in specified proportions:
``` 
mantle = m.mantle([lz,px],[8,2],['Lz','Px'])
```
The first list provides the lithology objects you just defined. The second list supplies their relative proportions (in this case 80% lherzolite, 20% pyroxenite). The final list provides the names you wish to call the lithologies. This will be used in data outputs.

To calculate the consequences of adiabatic decompression melting for this mantle object, the following method can be called, which will return a new `column` object:
```
column = mantle.AdiabticMelt_1D(1400.0)
```
The float provided is the mantle potential temperature (Tp) in °C. See the documentation for the method for how to change the starting pressure, final pressure, and number of steps. By default, the calculation starts at 8 GPa and proceeds all the way to the surface. If the calculation starts above the solidus, an interval of isobaric melting (conserving entropy) will be performed, as described in Matthews et al. (in review).

A handy way of visualising the results is to use the `plot()` method of the `column` object:
```
column.plot()
```
To access the calculation results directly:
- Temperature: `column.Temperature`
- Pressure: `column.P`
- Melt fraction of each lithology: `column.F`
- Aggregate melt fraction: `column.F_total`

To calculate crustal thickness (assuming passive decompression melting in a triangular spreading centre melting region) , the melt fractions must be integrated over the column:
```
tc = column.integrate_tri()
print(tc)
```
The method returns the crustal thickness in km, but the result is also stored in the column object. See documentation for further options for the integration calculation.

To estimate the crystallisation temperature of the melts, use the `MeltCrystallisationT()` method. The method returns two crystallisation temperature estimates, one for the shallowest melts, and one for the deepest melts. If no arguments are supplied to the method, the shallowest melts will be extracted from immediately below the calculated speading-centre crust, and melts will crystallise at the base of the calculated spreading-centre crust. Integration must have been performed beforehand. 
```
column.MeltCrystallisationT()
```
See documentation for how to apply this to method to other settings. 

## Online App
To see the module in action without installing the module or using any code, use our webapp at [pymeltapp.swmatthews.com](http://pymeltapp.swmatthews.com)
