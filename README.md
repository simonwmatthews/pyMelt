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
