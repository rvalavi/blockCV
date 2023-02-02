## Test environments
* MacOS on R v4.2.2
* Windows Server 2008 R2 SP1, R-release, 32/64 bit

## R CMD check results
There were no ERRORs or WARNINGs.

NOTES:
* checking installed package size ... NOTE
  installed size is  5.9Mb
  sub-directories of 1Mb or more:
    doc       3.6Mb
    extdata   1.9Mb

checking Rd cross-references ... NOTE Package unavailable to check Rd xrefs: ‘biomod2’

Response:
- the size of tallbal is less than 4.3MB
- 'biomod' package is not directly used in vignette or the package (the code related to this package is not executed in the vignette due to slow runtime).

## Downstream dependencies
All downstream dependencies are checked and they are fine.
