## Test environments

* Ubuntu 16.04 on travis-ci, R 3.2.4
* Local (Windows 2010), R 4.0.2
* win-builder (release, dev)
* Ubuntu 20.04 on GitHub Actions (release, dev)
* Mac OS 10.15 on GitHub Actions (release)
* Windows Server 2019 on GitHub Actions (release)


## R CMD check results

This is a resubmission, having set further tests to run conditionally and reduced build time of vignette.

There was 1 NOTE: 

### NOTEs

* Checking Rd cross-references
  + Package unavailable to check Rd xrefs: ‘gemtc’ (gemtc package not in "Imports" in the DESCRIPTION, but it would only be needed for cross-referencing in the documentation of the package and so is not required for package functionality)


## Downstream dependencies

There are no downstream dependencies (yet!)
