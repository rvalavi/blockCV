## Test environments
* local Ubuntu install, R v4.1.0
* Windows Server 2008 R2 SP1, R-release, 32/64 bit
* Oracle Solaris 10, x86, 32 bit, R-release

## R CMD check results
There were no ERRORs or WARNINGs.

NOTES:

* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Roozbeh Valavi <valavi.r@gmail.com>’

* checking installed package size ... NOTE
  installed size is  8.7Mb
  sub-directories of 1Mb or more:
    extdata   7.8Mb

* checking Rd cross-references ... NOTE
Package unavailable to check Rd xrefs: ‘biomod2’

Resonse:
- This is the main email I use and it will always be available.
- the size of tallbal is less than 5MB
- 'biomod' package is not directly used in vignette or the package (the code related to this package is will not be run in vignette due to slow runtime).

## Downstream dependencies
The both downstream dependencies (forestecology and sdmApp packages) are checked and they are fine.
