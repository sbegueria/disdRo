# disdRo

An R package for reading and handling meteorological particle size and
velocity distribution (PSVD), or disdrometric matrix data, from raw data
files generated by **Thies LPM** or **Ott Parsivel** optical disdrometers.
Index (SPEI)** and the **Standardized Precipitation Index (SPE)**.


## Details

Disdrometers are devices used to measure the particle size distribution and
velocity of falling hydrometeors.
They have multiple uses in automatic meterological networks and in climatology.
The Thies LPM and the Ott Parsivel are among the most common commercial
disdrometers.
This package allows reading and working with raw data files generated by
these two disdrometers, within R.

It currently conatains the workhorse functions `dsd_read()`, `dsd_integrate()`
and `dsd_plot()`, but other functionalities may be added in the future.
Support for other disdrometer types may be also be incorporated, althought
there are no fixed plans for that.


## Installation

Install the latest stable development version from GitHub:

```r
library(devtools)
install_github('sbegueria', 'disdRo', build_vignettes=TRUE)
```


## References

You can cite this references if you use the SPEI library on your work:





## Version history

### Version 0.1, July 2017 (current on github).

First release of the disdRo package.



## To do list for next version (work in progress)

- [ ] Make the functions more flexible to allow for other telegram configurations



## Any problems?

Feel free to [write an issue](https://github.com/sbegueria/disdRo/issues)
if you have any questions or problems.


## Copyright and license

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.