#########
Changelog
#########

The presently available stable version of pyMelt represents a largescale reorganisation of the
previous version of the module. The syntax has remained similar, but lithologies are now stored
in separate sub-modules.

v1.93
-----
* The `MeltingColumn` class has been condensed and all geological-setting specific functionality
  removed. It has also been renamed from `MeltingColumn_1D` to `MeltingColumn`
* There is a new `geosettings` module containing the geological-setting specific functionality
  removed from the `MeltingColumn` class.
* The default number of steps has been increased back to 1001 for the melting calculation.
* The calculation of crustal thickness in `geosettings.SpreadingCentre` now has a much smaller
  pressure decrement, so that small differences in crustal thickness can be determined more
  precisely.
