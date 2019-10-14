## Test environments

* Ubuntu 16.04 (on travis-ci), R 3.2.4
* Local Windows 2010, R 3.6.1 (devel and release)
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit (rhub)
* Fedora Linux, R-devel, clang, gfortran (rhub)


## R CMD check results

There were no ERRORs or WARNINGs. 

There were 2 NOTEs: 

* Checking CRAN incoming feasibility
  + New submission
  + Version contains large components (these are required to run the vignette)
* Checking Rd cross-references
  + Package unavailable to check Rd xrefs: ‘gemtc’ (gemtc package not in "Imports" in the DESCRIPTION, but it would only be needed for cross-referencing in the documentation of the package and so is not required for package functionality)


## Downstream dependencies

There are no downstream dependencies (yet!)
