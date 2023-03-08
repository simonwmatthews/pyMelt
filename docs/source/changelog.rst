#########
Changelog
#########

The presently available stable version of pyMelt represents a largescale reorganisation of the
previous version of the module. The syntax has remained similar, but lithologies are now stored
in separate sub-modules.

v2.01
-----
* Fixed a bug when calculating the crustal thickness in spreadingCentre geosetting if the
  column has a weighting function applied to it. The calculated chemistry is unaffected.

v2.00
-----
* Bug fixes related to preventing magma freezing and H2O exsolution.

v1.95
-----
* Hydrous melting is now supported via the `hydrousLithology` class
* The `chemistry` module has been added for trace element calculations
* By default, adiabatic decompression calculations will now start at the solidus with a fixed
  pressure decrement.
* This is essentially a move to v2.0, which is reserved for the version of pyMelt that coincides
  with manuscript publication.

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
