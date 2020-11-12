
# pyMelt documentation
S. Matthews (University of Iceland) and K. Wong (University of Leeds)

## Introduction

pyMelt is a python package for calculating the melting behaviour of mantle comprising multiple lithologies. The module implements the melting equations developed by Phipps Morgan (2001) to calculate the melting behaviour of mantle comprising any chosen lithology. 

Currently supported calculations:

* Adiabatic decompression melting
* Isobaric melting

Parameters that can be calculated:

* The geotherm for decompressing mantle
* Melt fractions for each lithology
* Crustal thickness for passive-upwelling at a mid-ocean ridge
* Crystallisation temperatures (following the method in Matthews et al., 2016)

## Getting started
First the pyMelt module must be imported. This can be achieved either by including the .py file in the same working directory:


```python
import pyMelt as m
```

If operating on a Linux installation, the latest version can be downloaded from the [code repository](https://github.com/simonwmatthews?tab=repositories) and installed locally:

```shell
$ git clone https://github.com/simonwmatthews/pyMelt
```

In addition, the following python packages must be installed and imported for pyMelt to operate:

* ```numpy```
* ```scipy```
* ```matplotlib```
* ```pandas```

## Lithology objects
pyMelt offers a number of different lithologies that can be thermodynamically modelled either separately or in combination. pyMelt includes the new parameterisations for KLB-1, KG1, and silica-saturated pyroxenite of Matthews et al. (in review). The lithologies included in the module are:

* ```LithologyKLB1```: KLB-1 lherzolite (Matthews et al., in review)
* ```LithologyKG1```: KG1 silica-undersaturated pyroxenite (Matthews et al., in review)
* ```LithologyPx```: silica-saturated pyroxenite (Matthews et al., in review)
* ```LithologyShorttle```: KG1 silica-undersaturated pyroxenite (Shorttle et al. 2014)
* ```LithologyKatz```: lherzolite (Katz et al., 2003)
* ```LithologySimple```: G2 pyroxenite (Pertermann & Hirschmann, 2002)
* ```LithologyNonMelting```: non-melting harzburgite

Each lithology is treated as a python object, and an instance can be assigned to a variable:


```python
lz = m.LithologyKLB1()
px = m.LithologyKG1()
hz = m.LithologyNonMelting()
```

Each lithology object contains methods describing its thermodynamic properties. As each lithology differs in composition, the thermodynamic/melting behaviour of each lithology varies; the methods available to each lithology object in pyMelt therefore also varies. For a full list of available methods please refer to the comments present within the pyMelt code. 
Some methods that are common to several lithologies include: 
* ```TSolidus(self, P)```: temperature of the lithology solidus (&deg;C) at a given pressure in GPa
* ```TLiquidus(self, P)```: temperature of the lithology liquidus (&deg;C) at a given pressure in GPa
* ```F(self, P, T)```: melt fraction of the lithology at a given pressure (GPa) and temperature (&deg;C)
* ```dTdF(self, P, T)```: dT/dF of the lithology at a constant pressure
* ```dTdP(self, P, T)```: dT/dP of the lithology at a constant melt fraction


```python
lz_solidus = lz.TSolidus(2)
lz_liquidus = lz.TLiquidus(2)
lz_F = lz.F(2, 1500)
print(lz_solidus, lz_liquidus, lz_F)
```

    1397.5733298733503 1911.2678949281221 0.249495821395242
    

## Mantle objects

A ```mantle``` object is constructed from one or multiple lithologies in specified proportions, and comprises three arguments.

* Argument 1: a list of the defined lithology objects to be considered in the melting calculation.
* Argument 2: a list of floats (of equivalent length to the list of Argument 1) comprising the relative proportions of the lithologies listed in Argument 1. The floats do not have to be normalised.
* Argument 3: a list of strings (of equivalent length to the other lists) comprising the names by which the lithologies of Argument 1 will be labelled. These strings will be used in data outputs. If ```False```, default names will be chosen.


```python
mantle = m.mantle([lz,px,hz],[6,2,2],['Lz', 'Px', 'Hz'])
```

The following objects and methods can be called from the ```mantle``` class, with all temperatures in &deg;C and all pressures in GPa:

``` python 
bulk_properties(self, P=False, T=False)
```
returns bulk thermodynamic properties of the solid or partially molten mantle (alpha, CP, rho; and labelled as such). If ```P```, ```T``` are ```False``` then the properties of the solid mantle will be returned.

```python
solidus_intersection(self, Tp)
```
returns the pressures (GPa) at which the solidus for each lithology will be intersected, in the order passed when creating the mantle class, assuming the mantle follows the solid adiabat for a potential temperature ```Tp``` up until that point.

```python
solidus_intersection_isobaric(self, P)
``` 
returns the pressures at which the solidus for each lithology will be intersected, in the order passed when creating the mantle class, assuming the mantle is heated isobarically at a pressure ```P```.

```python
adiabat(self, P, Tp)
```
returns the temperature of the solid mantle at a given pressure ```P``` and potential temperature ```Tp```.

```python
F(self, P, T)
```
returns the melt fraction of each lithology, in the order passed when creating the mantle class, at a given pressure ```P``` and temperature ```T```. 

```python
dFdP(self, P, T)
```
returns the dF/dP for each lithology, in the order passed when creating the mantle class, at a given pressure ```P``` and temperature ```T```, following Eq(26) of Phipps Morgan (2001).

```python
adiababatic_gradient(self, P, T)
```
returns dT/dP if melting has gone to completion (or has not started) for the bulk mantle at a given pressure ```P``` and temperature ```T```.

```python
dTdP(self, P, T, dFdP)
```
returns dT/dP using Eq(28) of Phipps Morgan (2001) at a given pressure ```P```, temperature ```T```, and ```dFdP```. Selects the lithology to use by the one with the largest increase in melt fraction with decompression (although this choice should not matter). The ```dFdP``` function is not recalled in order to save re-calculation of the same numbers. This function should be primarily used to check the behaviour of the pyMelt code.


```python
example_F = mantle.F(2, 1500)
print(example_F)
```

    [0.24949582 0.57529601 0.        ]
    

## Adiabatic decompression melting
To calculate the consequences of adiabatic decompression melting for this ```mantle``` object, the method ```AdiabaticMelt_1D``` can be called which will return a new ```column``` object:

```python
AdiabaticMelt_1D(self, Tp, Pstart=8.0, Pend=0.01, steps=1001, ReportSSS=True)
```
This function simultaneous integration of dFdP and dTdP to obtain the thermal gradient through the melting region. F of each lithology is then calculated along the P-T path. Integration is performed using a fourth-order Runge-Kutta algorithm. The P-T path is allowed to overstep the solidus on the step prior to the start of melting. Input parameters are:

* ```Tp```: potential temperature in &deg;C at which to perform the calculation.
* ```Pstart```; ```Pend```: the pressures in GPa at which the begin and stop upwelling respectively. Default pressures are 8.0 and 0.01 GPa.
* ```steps```: the number of dP increments the melting region is divided into. Default is 1001.
* ```ReportSSS```: print to the console if the starting pressure is above the solidus of one of the melting lithologies. Either way the ode will calculate melt fraction at this pressure through conservation of entropy. Can be set to ```False``` if this action is deliberate.
```AdiabaticMelt_1D``` returns a ```MeltingColumn_1D``` object, on which further calculations may be performed, comprising ```P```, ```T```, ```F``` results of the Runge-Kutta integration.

An example:


```python
column = mantle.AdiabaticMelt_1D(1400.0)
```

The P-T relationship of the given lithologies is visualised using the ```PlotBoundaries``` method, which generates two plots on one figure:

* Left subfigure depicting the P-T relationship of lithology solidii and liquidii,
* Right subfigure depicting the melt fraction of each lithology at fixed temperature as a function of pressure (default of 1600 &deg;C).

```python
PlotBoundaries(self, Pmax=8.0, Pmin=0.0, steps=1000, T_F=1600, show=True)
``` 
Input parameters are:

* ```Pmax```; ```Pmin```: the maximum and minimum pressures to display; default pressures are 8.0 and 0.0 GPa respectively. 
* ```steps```: discretization to use; default of 1000 discrete points.
* ```T_F```: temperature at which to calculate melt fractions; default of 1600 &deg;C.

This method returns a matplotlib.figure object, which can be adjusted and saved.

An example:


```python
Figure1 = mantle.PlotBoundaries()
```

## 1D melting column
A ```MeltingColumn_1D``` object is the result of appling ```AdiabaticMelt_1D``` to a mantle object, and contains information of the P-T-melt fraction properties of the adiabatically decompressing mantle. The ```MeltingColumn_1D``` object is visualised by using its ```plot()``` method, which generates a figure with two subplots:
* Left subfigure depicting the thermal gradient of the melting region in relation to the melting lithology solidii,
* Right subfigure depicting the total melt fraction in relation to the melt fractions of the melting lithologies.

```python
plot(self, solidii=True, show=True)
```
This function generates a matplotlib.figure object given a calculated ```MeltingColumn_1D``` object. The ```solidii``` parameter can be set to ```False``` to hide the solidii of the chosen lithologies on the final plot; likewise the plot itself can be hidden by setting ```show``` to ```False```.


```python
Figure2 = column.plot()
```

The calculation results, presented as pandas dataframes, can be directly accessed using the following commands:

* Temperature: ```column.Temperature```
* Pressure: ```column.P```
* Melt fraction of each lithology: ```column.F```
* Aggregate melt fraction: ```column.F_total```


```python
print(column.Temperature.iloc[0], column.F_total.iloc[-1])
```

    1570.371818980012 0.3106134938317263
    

## Calculating crustal thickness
Crustal thickness can be calculated assuming passive decompression melting in a triangular spreading centre melting region similar to that of a mid-ocean ridge. If this is the case the melt fractions must be integrated over the column using the ```integrate_tri()``` function to return a crustal thickness in km.


```python
tc = column.integrate_tri()
print(tc)
```

    11.882728732371918
    

The method also adds the following attributes to the ```MeltingColumn_1D``` object class:

* ```dtcdP```: results from Eq(6) (White et al., 1992) for the total melt fraction.
* ```tc_int```: integrated crustal thickness as a function of pressure (up to 0 GPa).
* ```tc_P_int```: pressure exerted by the integrated crustal thickness as a function of pressure (up to 0 GPa).
* ```tc```: integrated crustal thickness at the point where the pressure it exerts is equal to the calculation pressure.
* ```P_base_of_crust```: pressure at the base of the crust, at the point where the pressure the generated crust exerts is equal to the calculation pressure.
* ```tc_lithology_contributions_int```: integrated proportion of generated crust derived from each lithology as a function of pressure.
* ```tc_lithology_contributions``` integrated proportion of generated crust derived from each lithology where the pressure the generated crust exerts is equal to the calculation pressure.


```python
P_crust = column.P_base_of_crust
tc_lith_cont = column.tc_lithology_contributions
print(P_crust) 
print(tc_lith_cont)
```

    0.37753999999999976
    [0.18696586 0.81303414 0.        ]
    

## Melt liquidus temperature
pyMelt can be used to estimate the liquidus of mantle-derived melts (crystallisation temperature). This is achieved using the ```MeltCrystallisationT()``` method. Triangular integration must have been performed beforehand to achieve a liquidus temperature; else an error will be returned.

```python
MeltCrystallisationT(self, ShallowMeltP=False, MeltStorageP=False, liqdTdP=39.16)
``` 
This function returns two crystallisation temperature estimates, the first for the melts at the top of the melting column (the shallowest melts), and the second for the melts at the bottom of the melting column (the deepest melts).

The inputs are as follows:

* ```ShallowMeltP```: pressure at which the shallowest melt should be extracted. If ```False``` (as is default) then this will be taken to be the base of the crust.
* ```MeltStorageP```: pressure at which crystallisation occurs. If ```False``` (as is default), then this will be taken to be the base of the crust.
* ```liqdTdP```: the Clapeyron slope of the melt liquidus/olivine saturation surface in &deg;C/GPa. The default value is 39.16, from Eq(15) of Putirka (2008).


```python
Tcrys_values = column.MeltCrystallisationT()
print(Tcrys_values)
```

    (1251.591564249882, 1346.9009090794675)
    

## pyMelt_MultiNest
pyMelt can be used in conjunction with the MultiNest algorithm (Feroz and Hobson, 2008; Feroz et al., 2009, 2013) via its python frontend, pyMultinest (Buchner et al., 2014). This permits the inversion of measured data (e.g. crystallisation temperature, crustal thickness) to obtain unknowns (e.g. potential temperature) via Bayesian inference. More details of the inversion methods are provided in Matthews et al. (in review).

For pyMelt_MultiNest to work, MultiNest and pyMultinest must be installed. The user is directed to the [pyMultinest installation instructions](https://johannesbuchner.github.io/PyMultiNest/) for further guidance.

## Citing pyMelt
If pyMelt enables or aids your research please cite our publication:

Matthews, S., Wong, K., Shorttle, O., Edmonds, M., & Maclennan, J. (in review). Do olivine crystallization temperatures faithfully record mantle temperature variability?
https://doi.org/10.31223/osf.io/hqbgy
